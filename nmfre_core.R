# ============================================================
# NMF-RE (EU): Multiplier (Wild) Bootstrap + 1-step Newton
# - Function name: nmf.re
# - dfU cap: cap = dfU_cap_rate * (N*Q) unless dfU_cap (absolute) is supplied
# - Wild weights: exp_centered only  (w ~ Exp(1) - 1; mean 0, var 1)
#
# [ALIGN 2026-01-11]
# - X normalization now rescales BOTH Theta and U (identifiability-consistent)
# - MU updates use positive-part stabilization:
#     B_pos = pmax(ThetaA + U, 0) for X-step
#     Y_star_pos = pmax(pmax(Y,0) - XU, 0) for Theta-step
# - Bootstrap score weights centered: w = rexp(N) - 1
# - 1-step bootstrap replicates projected to Theta>=0
# - Removed obsolete wild_weight options; only "exp_centered" remains
# ============================================================

pmax_eps <- function(x, eps = 1e-10) pmax(x, eps)

safe_div <- function(num, den, eps = 1e-10) {
  res <- num / pmax_eps(den, eps)
  res[is.na(res) | is.infinite(res)] <- 0
  res
}

# ------------------------------------------------------------
# Normalize columns of X to sum 1, and rescale Theta AND U:
# X <- X D^{-1}, Theta <- D Theta, U <- D U
# ------------------------------------------------------------
normalize_X <- function(X, Theta, U = NULL, eps = 1e-12) {
  col_sums <- pmax(colSums(X), eps)
  Dinv <- diag(1 / col_sums, nrow = length(col_sums))
  D    <- diag(col_sums,     nrow = length(col_sums))

  Xn <- X %*% Dinv
  Th <- D %*% Theta

  if (!is.null(U)) {
    Un <- D %*% U
    return(list(X = Xn, Theta = Th, U = Un, scale = col_sums))
  }
  list(X = Xn, Theta = Th, scale = col_sums)
}

.nndsvdar <- function(Y, Q) {
  P <- nrow(Y)
  s <- svd(Y, nu = Q, nv = Q)
  W <- matrix(0, P, Q)
  W[, 1] <- sqrt(s$d[1]) * abs(s$u[, 1])

  if (Q > 1) {
    idx <- 2:Q
    U <- s$u[, idx, drop = FALSE]
    V <- s$v[, idx, drop = FALSE]
    D_sqrt <- sqrt(s$d[idx])

    U_p <- pmax(U, 0); U_n <- pmax(-U, 0)
    V_p <- pmax(V, 0); V_n <- pmax(-V, 0)

    norm_Up <- sqrt(colSums(U_p^2))
    norm_Un <- sqrt(colSums(U_n^2))
    norm_Vp <- sqrt(colSums(V_p^2))
    norm_Vn <- sqrt(colSums(V_n^2))

    Mp <- norm_Up * norm_Vp
    Mn <- norm_Un * norm_Vn

    eps <- .Machine$double.eps
    norm_Up_safe <- norm_Up; norm_Up_safe[norm_Up == 0] <- eps
    norm_Un_safe <- norm_Un; norm_Un_safe[norm_Un == 0] <- eps

    W_pos <- sweep(U_p, 2, (D_sqrt * Mp / norm_Up_safe), FUN = "*")
    W_neg <- sweep(U_n, 2, (D_sqrt * Mn / norm_Un_safe), FUN = "*")

    use_pos <- (Mp > Mn)
    W_combined <- W_pos
    W_combined[, !use_pos] <- W_neg[, !use_pos]
    W[, idx] <- W_combined
  }

  W[is.nan(W)] <- 0
  W[W < 0] <- 0

  avg_Y <- mean(Y)
  idx_zero <- which(W == 0)
  if (length(idx_zero) > 0) {
    W[idx_zero] <- stats::runif(length(idx_zero)) * avg_Y / 100
  }
  W
}

# ============================================================
# dfU helpers
# dfU(lambda) = N * sum_i d_i/(d_i + lambda), where d_i are eigenvalues of XtX
# ============================================================

.dfU_from_lambda <- function(d, N, lambda) {
  N * sum(d / (d + lambda))
}

.lambda_for_dfU_cap <- function(d, N, target,
                               lambda_min = 0,
                               lambda_max = 1e12,
                               tol = 1e-8) {
  if (.dfU_from_lambda(d, N, lambda_min) <= target) return(lambda_min)

  hi <- lambda_min
  df_hi <- Inf
  while (hi < lambda_max) {
    hi <- max(1e-12, if (hi == 0) 1e-6 else hi * 10)
    df_hi <- .dfU_from_lambda(d, N, hi)
    if (df_hi <= target) break
  }
  if (df_hi > target) return(lambda_max)

  f <- function(lam) .dfU_from_lambda(d, N, lam) - target
  out <- stats::uniroot(f, lower = lambda_min, upper = hi, tol = tol)
  out$root
}

# ============================================================
# Main function: nmf.re
# ============================================================
# ============================================================
# Main function: nmf.re
#  - Added: theta_p_side ("two.sided" or "one.sided") for boundary null p-values
#  - No auto-detection: if one.sided, apply to ALL Theta coefficients
# ============================================================
nmf.re <- function(
  Y, A = NULL, Q = 2,

  # ---- initialization ----
  X.ini = NULL, Theta.ini = NULL, U.ini = NULL,
  prefix = "Basis",

  # ---- variance handling ----
  sigma2 = 1,
  sigma2_update = TRUE,
  sigma2_update_start = 50,
  sigma2_update_every = 10,
  sigma2_update_rate = 0.05,
  sigma2_min = 1e-12,
  sigma2_max = 1e12,

  # tau2 moment matching (disabled when dfU cap is on)
  tau2 = 1,
  tau2_update = TRUE,
  tau2_update_start = 1,
  tau2_update_every = 1,
  tau2_update_rate = 0.2,
  tau2_min = 1e-12,
  tau2_max = 1e12,

  # ---- dfU control (NQ only) ----
  dfU_control = c("cap", "off"),
  dfU_cap = NULL,           # absolute cap (optional)
  dfU_cap_rate = 0.10,      # cap = rate * (N*Q) (e.g., 0.05)
  dfU_enforce_every = 1,
  dfU_lambda_max = 1e12,

  # ---- optimization ----
  maxit = 50000,
  epsilon = 1e-8,
  verbose = TRUE,
  seed = 1,

  # ---- inference options (Gaussian working) ----
  inference = TRUE,
  sigma2_hat = NULL,

  df_sigma = c("PN-df", "PN", "PN-QR"),
  theta_info_mode = c("IminusH", "plain"),
  sandwich = TRUE,
  sandwich_cr1 = TRUE,
  se_rule = c("sandwich", "naive", "max"),

  cov_ridge = 1e-8,

  # ---- p-value side for boundary null (Theta >= 0) ----
  theta_p_side = c("one.sided","two.sided"),

  # ---- multiplier (wild) bootstrap (executed when inference=TRUE) ----
  wild_B = 500,
  wild_level = 0.95,
  wild_weight = c("exp_centered"),
  wild_seed = 123,

  # ---- post-run message for dfU cap ----
  post_message = TRUE,

  # ---- objective trace (optional; for debugging/plots) ----
  save_obj_trace = TRUE
) {

  set.seed(seed)
  .eps <- 1e-10

  df_sigma <- match.arg(df_sigma)
  theta_info_mode <- match.arg(theta_info_mode)
  se_rule <- match.arg(se_rule)
  wild_weight <- match.arg(wild_weight)
  dfU_control <- match.arg(dfU_control)
  theta_p_side <- match.arg(theta_p_side)

  if (!is.numeric(dfU_cap_rate) || length(dfU_cap_rate) != 1 ||
      !is.finite(dfU_cap_rate) || dfU_cap_rate <= 0) {
    stop("'dfU_cap_rate' must be a positive finite number (e.g., 0.05).")
  }

  # ---- dimensions ----
  P <- nrow(Y); N <- ncol(Y)
  if (is.null(A)) A <- matrix(1, 1, N)
  if (!is.matrix(A)) A <- as.matrix(A)
  stopifnot(ncol(A) == N)
  K <- nrow(A)

  # ---- init ----
  if (is.null(X.ini)) {
    X <- .nndsvdar(pmax(Y, 0), Q)
  } else {
    X <- X.ini
    stopifnot(nrow(X) == P, ncol(X) == Q)
  }

  if (is.null(Theta.ini)) {
    Theta <- matrix(1, nrow = Q, ncol = K)
  } else {
    Theta <- Theta.ini
    stopifnot(nrow(Theta) == Q, ncol(Theta) == K)
  }

  if (is.null(U.ini)) {
    U <- matrix(0, Q, N)
  } else {
    U <- U.ini
    stopifnot(nrow(U) == Q, ncol(U) == N)
  }

  # ---- normalize X columns to sum 1 (rescale Theta AND U) ----
  normed <- normalize_X(pmax(X, .eps), pmax(Theta, .eps), U)
  X <- normed$X
  Theta <- normed$Theta
  U <- normed$U

  clip_val <- function(x, xmin, xmax) min(max(x, xmin), xmax)

  # ---- convergence bookkeeping ----
  obj <- NA_real_
  obj_prev <- Inf
  rel_change <- NA_real_

  converged <- FALSE
  stop_reason <- "maxit"
  iter_done <- 0

  obj_trace <- if (isTRUE(save_obj_trace)) rep(NA_real_, maxit) else NULL
  rss_trace <- if (isTRUE(save_obj_trace)) rep(NA_real_, maxit) else NULL

  # ---- dfU diagnostics ----
  dfU_last <- NA_real_
  dfU_cap_last <- NA_real_
  lambda_last <- NA_real_

  dfU_hit_cap <- FALSE
  dfU_hit_iter <- integer(0)

  # =========================================================
  # main loop
  # =========================================================
  for (iter in 1:maxit) {
    iter_done <- iter

    # (0) current components
    ThetaA <- Theta %*% A  # Q x N

    # (1) U-step: ridge BLUP with lambda = sigma2/tau2
    tau2 <- clip_val(tau2, tau2_min, tau2_max)
    sigma2 <- clip_val(sigma2, sigma2_min, sigma2_max)
    lambda <- sigma2 / tau2

    XtX <- crossprod(X)  # Q x Q

    # ---- dfU cap enforcement ----
    if (dfU_control == "cap" && (iter %% dfU_enforce_every == 0)) {

      dfU_cap_now <- if (!is.null(dfU_cap)) dfU_cap else dfU_cap_rate * (N * Q)
      dfU_cap_now <- max(dfU_cap_now, 1)

      d <- eigen(XtX, symmetric = TRUE, only.values = TRUE)$values
      d <- pmax(d, 0)

      dfU_now <- .dfU_from_lambda(d, N, lambda)

      if (dfU_now > dfU_cap_now) {
        lambda_new <- .lambda_for_dfU_cap(
          d, N, target = dfU_cap_now,
          lambda_min = pmax(lambda, 0),
          lambda_max = dfU_lambda_max
        )
        lambda <- lambda_new
        tau2 <- clip_val(sigma2 / pmax(lambda, 1e-12), tau2_min, tau2_max)
        dfU_now <- .dfU_from_lambda(d, N, lambda)

        dfU_hit_cap <- TRUE
        dfU_hit_iter <- c(dfU_hit_iter, iter)
      }

      dfU_last <- dfU_now
      dfU_cap_last <- dfU_cap_now
    }

    lambda_last <- lambda

    M <- XtX + diag(pmax(lambda, 1e-12), Q)
    cholM <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(cholM)) stop("Cholesky failed: XtX + lambda I not SPD. Increase lambda (or decrease tau2).")

    # residual without U
    fit_fixed <- X %*% ThetaA
    R0 <- Y - fit_fixed

    # solve for each column u_n
    for (n in 1:N) {
      rhs <- crossprod(X, R0[, n, drop = FALSE])  # Q x 1
      u <- backsolve(cholM, forwardsolve(t(cholM), rhs))
      U[, n] <- as.numeric(u)
    }

    # center U (identifiability)
    U <- sweep(U, 1, rowMeans(U), "-")

    # (1.5) tau2 update by moment matching (disabled when dfU cap is on)
    if (isTRUE(tau2_update) && dfU_control != "cap" &&
        iter >= tau2_update_start && (iter %% tau2_update_every == 0)) {

      term1 <- sum(U^2) / (N * Q)

      Minv <- tryCatch(chol2inv(cholM), error = function(e) NULL)
      if (!is.null(Minv)) {
        term2 <- sigma2 * sum(diag(Minv)) / Q
        tau2_new <- clip_val(term1 + term2, tau2_min, tau2_max)
        tau2 <- (1 - tau2_update_rate) * tau2 + tau2_update_rate * tau2_new
      }
    }

    # (2) X-step (EU MU) with positive-part stabilization
    Y_tilde <- pmax(Y, 0)
    B_sem <- ThetaA + U
    B_pos <- pmax(B_sem, 0)

    numX <- Y_tilde %*% t(B_pos)
    denX <- X %*% (B_pos %*% t(B_pos))
    X <- X * safe_div(numX, denX, eps = .eps)
    X <- pmax(X, .eps)

    # normalize columns sum 1; absorb scale into Theta AND U
    normed <- normalize_X(X, Theta, U)
    X <- pmax(normed$X, .eps)
    Theta <- pmax(normed$Theta, .eps)
    U <- normed$U

    # (3) Theta-step (EU MU) with positive-part stabilization
    Y_star <- Y_tilde - X %*% U
    Y_star_pos <- pmax(Y_star, 0)

    numT <- crossprod(X, Y_star_pos) %*% t(A)
    denT <- (crossprod(X, X) %*% Theta) %*% (A %*% t(A))
    Theta <- Theta * safe_div(numT, denT, eps = .eps)
    Theta <- pmax(Theta, .eps)

    # (3.5) sigma2 update from residual (optional; iteration-scale)
    if (isTRUE(sigma2_update) && iter >= sigma2_update_start && (iter %% sigma2_update_every == 0)) {
      ThetaA <- Theta %*% A
      fit <- X %*% (ThetaA + U)
      R <- Y - fit
      sigma2_new <- clip_val(mean(R^2), sigma2_min, sigma2_max)
      sigma2 <- (1 - sigma2_update_rate) * sigma2 + sigma2_update_rate * sigma2_new
    }

    # (4) objective & convergence (uses raw ThetaA+U, not truncated)
    ThetaA <- Theta %*% A
    fit <- X %*% (ThetaA + U)
    R <- Y - fit

    lambda_obj <- sigma2 / clip_val(tau2, tau2_min, tau2_max)
    obj <- sum(R^2) + lambda_obj * sum(U^2)
    rss <- sum(R^2)

    if (isTRUE(save_obj_trace)) {
      obj_trace[iter] <- obj
      rss_trace[iter] <- rss
    }

    if (is.finite(obj_prev)) {
      rel_change <- abs(obj_prev - obj) / (abs(obj_prev) + .eps)
    } else {
      rel_change <- NA_real_
    }

    if (!is.finite(obj)) {
      stop_reason <- "nonfinite_obj"
      converged <- FALSE
      break
    }

    if (is.finite(obj_prev) && rel_change < epsilon) {
      stop_reason <- "epsilon"
      converged <- TRUE
      break
    }

    obj_prev <- obj

    if (verbose && iter %% 100 == 0) {
      cat(sprintf("iter=%d  obj=%.6g  rel=%.3g  sigma2=%.4g  tau2=%.4g  lambda=%.4g",
                  iter, obj, rel_change, sigma2, tau2, sigma2 / tau2))
      if (dfU_control == "cap" && is.finite(dfU_last)) {
        cat(sprintf("  dfU=%.2f cap=%.2f lam_enf=%.4g", dfU_last, dfU_cap_last, lambda_last))
      }
      cat("\n")
    }
  }

  # ---- post loop: finalize convergence fields robustly ----
  obj_final <- obj
  rel_change_final <- rel_change

  if (iter_done < maxit && identical(stop_reason, "maxit")) {
    stop_reason <- "unknown_break"
    converged <- FALSE
  }

  # ---- reorder basis (optional heuristic) ----
  if (ncol(X) > 1) {
    w_ord <- matrix((1:P) / P, nrow = 1)
    score <- as.numeric(w_ord %*% X)
    index <- order(score)
    X <- X[, index, drop = FALSE]
    Theta <- Theta[index, , drop = FALSE]
    U <- U[index, , drop = FALSE]
  }

  # ---- final fitted matrices ----
  ThetaA <- Theta %*% A
  B_fixed <- ThetaA
  B_blup_raw <- ThetaA + U

  XB <- X %*% B_fixed
  XB.blup <- X %*% B_blup_raw

  # ---- final dfU recomputation (ensure consistency with final X, lambda) ----
  lambda_final <- sigma2 / clip_val(tau2, tau2_min, tau2_max)
  XtX_final <- crossprod(X)
  d_final <- eigen(XtX_final, symmetric = TRUE, only.values = TRUE)$values
  d_final <- pmax(d_final, 0)
  dfU_final <- .dfU_from_lambda(d_final, N, lambda_final)

  if (dfU_control == "cap") {
    dfU_cap_last <- if (!is.null(dfU_cap)) max(dfU_cap, 1) else max(dfU_cap_rate * (N * Q), 1)
  } else {
    dfU_cap_last <- NA_real_
  }
  dfU_last <- dfU_final
  lambda_last <- lambda_final

  # ---- names ----
  colnames(X) <- paste0(prefix, 1:ncol(X))
  rownames(Theta) <- colnames(X)
  colnames(Theta) <- rownames(A)
  rownames(U) <- rownames(Theta)
  colnames(U) <- colnames(Y)

  # =========================================================
  # Theta inference
  # =========================================================
  Theta.vec.cov <- NULL
  Theta.se <- NULL                # kept for backward compatibility (Hessian/Sandwich)
  Theta.se_hess <- NULL           # explicit Hessian/Sandwich SE
  Theta.se_boot <- NULL           # bootstrap SE
  coefficients <- NULL
  sigma2_used <- sigma2_hat

  Theta.ci.lower <- NULL
  Theta.ci.upper <- NULL
  Theta.boot.sd <- NULL

  if (isTRUE(inference)) {

    # Y_star for inference (conditioning on X,U)
    Y_star_inf <- Y - X %*% U
    R_theta <- Y_star_inf - X %*% (Theta %*% A)
    RSS_inf <- sum(R_theta^2)

    lambda_inf <- sigma2 / pmax(tau2, 1e-12)

    XtX_now <- crossprod(X)
    AAt <- A %*% t(A)

    M_inf <- XtX_now + diag(pmax(lambda_inf, 1e-12), Q)
    cholM_inf <- tryCatch(chol(M_inf), error = function(e) NULL)
    if (!is.null(cholM_inf)) {
      Minv <- chol2inv(cholM_inf)
    } else {
      Minv <- tryCatch(solve(M_inf), error = function(e) MASS::ginv(M_inf))
    }

    trH <- sum(diag(Minv %*% XtX_now))
    dfU_inf <- N * trH

    if (is.null(sigma2_used)) {
      if (df_sigma == "PN") {
        denom <- P * N
      } else if (df_sigma == "PN-QR") {
        denom <- max(P * N - Q * K, 1)
      } else { # "PN-df"
        dfTheta <- Q * K
        denom <- max(P * N - dfU_inf - dfTheta, 1)
      }
      sigma2_used <- RSS_inf / denom
    }

    if (theta_info_mode == "plain") {
      Info_core <- kronecker(AAt, XtX_now)
    } else {
      Xt_IH_X <- XtX_now - XtX_now %*% Minv %*% XtX_now
      Info_core <- kronecker(AAt, Xt_IH_X)
    }

    Info <- Info_core / max(sigma2_used, 1e-12)
    Info <- Info + diag(cov_ridge, nrow(Info))

    Hinv <- tryCatch(solve(Info), error = function(e) MASS::ginv(Info))
    V_naive <- Hinv

    V_sand <- NULL
    if (isTRUE(sandwich)) {
      Xt <- t(X)
      J <- matrix(0, Q * K, Q * K)

      for (n in 1:N) {
        a_n <- A[, n, drop = FALSE]
        g_n <- Xt %*% R_theta[, n, drop = FALSE]
        S_n <- -(g_n %*% t(a_n)) / max(sigma2_used, 1e-12)
        s_n <- as.vector(S_n)
        J <- J + tcrossprod(s_n)
      }

      if (isTRUE(sandwich_cr1) && N > 1) {
        J <- (N / (N - 1)) * J
      }
      V_sand <- Hinv %*% J %*% Hinv
    }

    if (se_rule == "naive") {
      Theta.vec.cov <- V_naive
    } else if (se_rule == "sandwich") {
      Theta.vec.cov <- if (is.null(V_sand)) V_naive else V_sand
    } else { # "max"
      if (is.null(V_sand)) {
        Theta.vec.cov <- V_naive
      } else {
        d1 <- pmax(diag(V_naive), 0)
        d2 <- pmax(diag(V_sand), 0)
        dmax <- pmax(d1, d2)
        Theta.vec.cov <- V_naive
        diag(Theta.vec.cov) <- dmax
      }
    }

    se_vec <- sqrt(pmax(diag(Theta.vec.cov), 0))
    Theta.se_hess <- matrix(se_vec, nrow = Q, ncol = K, byrow = FALSE)
    dimnames(Theta.se_hess) <- dimnames(Theta)

    # Backward compatible: Theta.se is the Hessian/Sandwich SE (used for z/p)
    Theta.se <- Theta.se_hess

    # ---- Multiplier (wild) bootstrap (centered Exp weights) ----
    set.seed(wild_seed)

    Xt <- t(X)
    score_mat <- matrix(0, Q * K, N)

    for (n in 1:N) {
      a_n <- A[, n, drop = FALSE]
      r_n <- R_theta[, n, drop = FALSE]
      g_n <- Xt %*% r_n
      G_n <- -(g_n %*% t(a_n)) / max(sigma2_used, 1e-12)
      score_mat[, n] <- as.vector(G_n)
    }

    theta_hat_vec <- as.vector(Theta)
    Theta_boot <- matrix(NA_real_, nrow = Q * K, ncol = wild_B)

    for (b in 1:wild_B) {
      # exp_centered: mean 0, var 1
      w <- stats::rexp(N, rate = 1) - 1
      grad_b <- as.vector(score_mat %*% w)

      theta_b <- theta_hat_vec - as.vector(Hinv %*% grad_b)

      # project to Theta>=0 (paper-consistent)
      theta_b <- pmax(theta_b, 0)

      Theta_boot[, b] <- theta_b
    }

    alpha <- 1 - wild_level
    lo_q <- alpha / 2
    hi_q <- 1 - alpha / 2

    lo <- apply(Theta_boot, 1, stats::quantile, probs = lo_q, na.rm = TRUE, names = FALSE)
    hi <- apply(Theta_boot, 1, stats::quantile, probs = hi_q, na.rm = TRUE, names = FALSE)

    Theta.ci.lower <- matrix(lo, nrow = Q, ncol = K, byrow = FALSE)
    Theta.ci.upper <- matrix(hi, nrow = Q, ncol = K, byrow = FALSE)
    dimnames(Theta.ci.lower) <- dimnames(Theta)
    dimnames(Theta.ci.upper) <- dimnames(Theta)

    sd_vec <- apply(Theta_boot, 1, stats::sd, na.rm = TRUE)
    Theta.boot.sd <- matrix(sd_vec, nrow = Q, ncol = K, byrow = FALSE)
    dimnames(Theta.boot.sd) <- dimnames(Theta)

    Theta.se_boot <- Theta.boot.sd

    # ---- coefficients (z/p computed with Hessian/Sandwich SE) ----
    Estimate <- as.vector(Theta)
    SE_hess_vec <- as.vector(Theta.se_hess)
    SE <- SE_hess_vec
    BSE_vec <- as.vector(Theta.boot.sd)

    # z-value
    z_value <- rep(NA_real_, length(Estimate))
    ok <- is.finite(Estimate) & is.finite(SE) & (SE > 0)
    z_value[ok] <- Estimate[ok] / SE[ok]

    # Wald (chi-square, df=1) for two-sided style (as in current implementation)
    Wald <- rep(NA_real_, length(Estimate))
    Wald[ok] <- z_value[ok]^2

    # two-sided p-value (current default)
    p_two <- rep(NA_real_, length(Estimate))
    p_two[ok] <- 1 - stats::pchisq(Wald[ok], df = 1)

    # one-sided p-value for boundary null: H0: Theta=0 vs H1: Theta>0
    p_one <- rep(NA_real_, length(Estimate))
    p_one[ok] <- stats::pnorm(z_value[ok], lower.tail = FALSE)

    # choose
    p_value <- if (theta_p_side == "one.sided") p_one else p_two

    coefficients <- data.frame(
      Variable = rep(colnames(Theta), each = Q),
      Basis    = rep(rownames(Theta), times = K),
      Estimate = Estimate,
      SE       = SE,
      BSE      = BSE_vec,
      z_value  = z_value,
      Wald     = Wald,
      p_value  = p_value,
      p_two_sided = p_two,
      p_one_sided = p_one,
      p_side_used = theta_p_side,
      row.names = NULL
    )

    if (!is.null(Theta.ci.lower) && !is.null(Theta.ci.upper)) {
      coefficients$CI_low  <- as.vector(Theta.ci.lower)
      coefficients$CI_high <- as.vector(Theta.ci.upper)
    }
  }

  # ---- Post-run message: dfU cap ----
  if (isTRUE(post_message) && dfU_control == "cap" &&
      is.finite(dfU_last) && is.finite(dfU_cap_last)) {
    dfU_cap_now <- dfU_cap_last
    if (isTRUE(dfU_hit_cap)) {
      message(sprintf(
        "[dfU-cap] Activated: dfU reached the cap during optimization.\nFinal dfU = %.2f (%.1f%% of NQ), cap = %.2f (%.1f%% of NQ).\nFirst activated at iter %d (hits=%d).",
        dfU_last, 100 * (dfU_last / (N * Q)),
        dfU_cap_now, 100 * (dfU_cap_now / (N * Q)),
        dfU_hit_iter[1], length(dfU_hit_iter)
      ))
    } else {
      message(sprintf(
        "[dfU-cap] Not activated: dfU never reached the cap.\nFinal dfU = %.2f (%.1f%% of NQ), cap = %.2f (%.1f%% of NQ).",
        dfU_last, 100 * (dfU_last / (N * Q)),
        dfU_cap_now, 100 * (dfU_cap_now / (N * Q))
      ))
    }
  }

  # ---- convenience probabilities ----
  colnorm_prob <- function(M, eps = 1e-12) {
    cs <- colSums(M)
    sweep(M, 2, pmax(cs, eps), "/")
  }

  out <- list(
    X = X,
    Theta = Theta,
    U = U,

    sigma2 = sigma2,
    tau2 = tau2,
    lambda = sigma2 / tau2,

    # convergence diagnostics
    converged = converged,
    stop_reason = stop_reason,
    iter = iter_done,
    maxit = maxit,
    epsilon = epsilon,
    obj_final = obj_final,
    rel_change_final = rel_change_final,
    obj_trace = if (isTRUE(save_obj_trace)) obj_trace[seq_len(iter_done)] else NULL,
    rss_trace = if (isTRUE(save_obj_trace)) rss_trace[seq_len(iter_done)] else NULL,

    # dfU diagnostics (final-consistent)
    dfU = dfU_last,
    dfU_cap = dfU_cap_last,
    dfU_cap_rate = dfU_cap_rate,
    lambda_enforced = lambda_last,
    dfU_hit_cap  = dfU_hit_cap,
    dfU_hit_iter = dfU_hit_iter,
    dfU_frac     = if (is.finite(dfU_last)) dfU_last / (N * Q) else NA_real_,
    dfU_cap_frac = if (is.finite(dfU_cap_last)) dfU_cap_last / (N * Q) else NA_real_,

    # fitted matrices (raw BLUP scores allowed to be real-valued)
    C = Theta,
    B = Theta %*% A,
    B.prob = colnorm_prob(pmax(Theta %*% A, 0)),
    B.blup = (Theta %*% A) + U,
    B.blup.pos = pmax((Theta %*% A) + U, 0),
    B.blup.prob = colnorm_prob(pmax((Theta %*% A) + U, 0)),
    XB = XB,
    XB.blup = XB.blup,

    # inference outputs
    sigma2_used = sigma2_used,
    Theta.vec.cov = Theta.vec.cov,
    Theta.se = Theta.se,                 # Hessian/Sandwich (compat)
    Theta.se_hess = Theta.se_hess,       # explicit
    Theta.se_boot = Theta.se_boot,       # explicit
    coefficients = coefficients,
    Theta.ci.lower = Theta.ci.lower,
    Theta.ci.upper = Theta.ci.upper,
    Theta.boot.sd  = Theta.boot.sd,

    # record inference option
    theta_p_side = theta_p_side
  )

  class(out) <- "nmfre"
  out
}



# ============================================================
# summary method
# - SE column shown as "HSE (BSE)" when BSE exists
# ============================================================
summary.nmfre <- function(object,
                          digits = 4,
                          coef_digits = 3,
                          p_digits = 4,
                          show_ci = FALSE,
                          max_rows = Inf,
                          ...) {

  x <- object

  cat("Summary: NMF-RE (EU)\n")
  cat("------------------------------------------------------------\n")

  # ---- convergence block (ALWAYS show) ----
  cat(sprintf("Converged: %s\n", ifelse(isTRUE(x$converged), "YES", "NO")))
  cat(sprintf("Stop reason: %s\n", if (!is.null(x$stop_reason)) x$stop_reason else "NA"))
  cat(sprintf("Iterations: %s / %s\n",
              if (!is.null(x$iter)) x$iter else "NA",
              if (!is.null(x$maxit)) x$maxit else "NA"))
  cat(sprintf("Objective (final): %s\n",
              if (!is.null(x$obj_final) && is.finite(x$obj_final)) format(x$obj_final, digits = 6) else "NA"))
  cat(sprintf("Relative change (final): %s  (epsilon = %s)\n",
              if (!is.null(x$rel_change_final) && is.finite(x$rel_change_final)) format(x$rel_change_final, digits = 3) else "NA",
              if (!is.null(x$epsilon) && is.finite(x$epsilon)) format(x$epsilon, digits = 3) else "NA"))
  cat("\n")

  # ---- size/variance ----
  cat(sprintf("Dimensions: X(%d,%d), Theta(%d,%d), U(%d,%d)\n",
              nrow(x$X), ncol(x$X),
              nrow(x$Theta), ncol(x$Theta),
              nrow(x$U), ncol(x$U)))
  cat(sprintf("sigma2 = %.6g, tau2 = %.6g, lambda = %.6g\n",
              x$sigma2, x$tau2, x$lambda))

  # ---- dfU cap diagnostics ----
  if (!is.null(x$dfU_cap) && is.finite(x$dfU_cap)) {
    cat(sprintf(
      "dfU(cap): dfU = %.2f (%.1f%% of NQ), cap = %.2f (%.1f%% of NQ), rate = %.3g\n",
      x$dfU, 100 * x$dfU_frac,
      x$dfU_cap, 100 * x$dfU_cap_frac,
      x$dfU_cap_rate
    ))
    cat(sprintf("dfU-cap activated: %s",
                ifelse(isTRUE(x$dfU_hit_cap), "YES", "NO")))
    if (isTRUE(x$dfU_hit_cap)) {
      it0 <- if (length(x$dfU_hit_iter)) x$dfU_hit_iter[1] else NA_integer_
      cat(sprintf(" (first at iter %s; hits=%d)\n", it0, length(x$dfU_hit_iter)))
    } else {
      cat("\n")
    }
  }
  cat("\n")

  # ---- Theta ----
  cat("Theta (rounded):\n")
  print(round(x$Theta, digits))
  cat("\n")

  # ---- coefficients ----
  if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
    cat("Theta regression coefficients:\n")

    tmp <- x$coefficients
    tmp$Estimate <- round(tmp$Estimate, coef_digits)

    # SE display: "HSE (BSE)" if BSE exists
    if ("BSE" %in% names(tmp)) {
      hse <- tmp$SE
      bse <- tmp$BSE
      tmp$SE <- ifelse(
        is.finite(bse),
        sprintf(paste0("%.", coef_digits, "f (%.", coef_digits, "f)"),
                round(hse, coef_digits), round(bse, coef_digits)),
        sprintf(paste0("%.", coef_digits, "f"),
                round(hse, coef_digits))
      )
    } else {
      tmp$SE <- round(tmp$SE, coef_digits)
    }

    tmp$z_value  <- round(tmp$z_value, 2)
    tmp$p_value  <- sprintf(paste0("%.", p_digits, "f"), tmp$p_value)

    keep_cols <- c("Variable", "Basis", "Estimate", "SE", "z_value", "p_value")

    if (isTRUE(show_ci) && all(c("CI_low", "CI_high") %in% names(tmp))) {
      tmp$CI_low  <- round(tmp$CI_low, coef_digits)
      tmp$CI_high <- round(tmp$CI_high, coef_digits)
      keep_cols <- c(keep_cols, "CI_low", "CI_high")
    }

    tmp <- tmp[, keep_cols, drop = FALSE]

    if (is.finite(max_rows) && nrow(tmp) > max_rows) {
      cat(sprintf("(showing first %d of %d rows)\n", max_rows, nrow(tmp)))
      tmp <- tmp[seq_len(max_rows), , drop = FALSE]
    }

    print(tmp, row.names = FALSE)
  }

  invisible(x)
}


# ============================================================
# Utility: scan dfU_cap_rate
# ============================================================
scan_dfU_cap_rate <- function(
  rates = (1:10) / 100,
  Y, A, Q,
  X.ini = NULL, Theta.ini = NULL, U.ini = NULL,
  verbose = FALSE,
  ...,

  ok_upper = 1.00,
  warn_upper = 0.80,

  digits = list(dfU = 2, ratio = 2, lambda = 4, tau2 = 4, sigma2 = 4)
) {
  N <- ncol(Y)
  NQ <- N * Q

  one_run <- function(rate) {

    r <- nmf.re(
      Y = Y, A = A, Q = Q,
      X.ini = X.ini, Theta.ini = Theta.ini, U.ini = U.ini,
      verbose = verbose,
      dfU_control = "cap",
      dfU_cap_rate = rate,
      ...
    )

    ratio <- if (is.finite(r$dfU) && is.finite(r$dfU_cap) && r$dfU_cap > 0) {
      r$dfU / r$dfU_cap
    } else {
      NA_real_
    }

    status <- NA_character_
    if (is.finite(ratio)) {
      if (ratio >= 0.99) status <- "NG"
      if (ratio >= 0.80 && ratio < 0.99) status <- "BETTER"
      if (ratio > 0 && ratio < 0.80) status <- "GOOD"
    }

    data.frame(
      rate = rate,
      dfU_cap = r$dfU_cap,
      dfU = round(r$dfU, 2),
      ratio = round(ratio, 2),
      status = status,
      converged = isTRUE(r$converged),
      tau2 = round(r$tau2, 3),
      sigma2 = round(r$sigma2, 3),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, lapply(rates, one_run))
  rownames(out) <- NULL
  out
}


.dfU_from_lambda <- function(d, N, lambda) {
  N * sum(d / (d + lambda))
}

.lambda_for_dfU_cap <- function(d, N, target,
                               lambda_min = 0,
                               lambda_max = 1e12,
                               tol = 1e-8) {
  if (.dfU_from_lambda(d, N, lambda_min) <= target) return(lambda_min)

  hi <- lambda_min
  df_hi <- Inf
  while (hi < lambda_max) {
    hi <- max(1e-12, if (hi == 0) 1e-6 else hi * 10)
    df_hi <- .dfU_from_lambda(d, N, hi)
    if (df_hi <= target) break
  }
  if (df_hi > target) return(lambda_max)

  f <- function(lam) .dfU_from_lambda(d, N, lam) - target
  out <- stats::uniroot(f, lower = lambda_min, upper = hi, tol = tol)
  out$root
}


  # ---- convenience probabilities ----
  colnorm_prob <- function(M, eps = 1e-12) {
    cs <- colSums(M)
    sweep(M, 2, pmax(cs, eps), "/")
  }
