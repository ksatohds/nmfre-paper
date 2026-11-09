# NMF-RE estimation core -- 2026-07-18.
# Supersedes nmfre_core20260627_fullEM.R.  Same outer-inner ECM, same EM-form variance
# updates, same semi-NMF X-step with the posterior-variance term.  TWO corrections:
#
#  (1) The row-centering of U is REMOVED.  The U/Theta indeterminacy of the model is
#      U -> U + Delta A , Theta -> Theta - Delta  (any Q x K Delta), so the condition that
#      removes it is  U A' = 0 , not the row-centering U 1 = 0.  The alternating scheme
#      attains U A' = 0 by itself at its fixed point; imposing the row-centering on top of
#      that DESTROYS it whenever 1' is not in the row space of A.  Consequence, measured on
#      the paper's examples:
#         centered  : ||U A'|| = 63.9 (HS) / 1691.6 (OD);  profiled score 1.14 / 27.2 ;
#                     Thetahat differs from the root of the profiled estimating equation by
#                     9.1% (HS) / 19.8% (OD)
#         uncentered: ||U A'|| = 0 ;  profiled score 8e-15 ;  Thetahat = that root exactly
#      Sections 4 and Appendix B derive the inference by profiling the UNCENTERED ridge
#      solution, so only the uncentered fit is the estimator the inference describes.
#      (Orthodont was unaffected either way: 1' IS in the row space of its A, so the fixed
#      point is already row-centered and ||U A'|| = 0.)
#
#  (2) The initial Theta is no longer clipped to be positive when theta_sign = "free".
#      The old line  normalize_X(pmax(X,.eps), pmax(Theta,.eps), U)  applied pmax to Theta
#      regardless of the sign convention, so any negative entry of a supplied Theta.ini was
#      reset to +eps.  Harmless for a non-negative NMF initialisation, but the C2 full-refit
#      bootstrap warm-starts from the fitted Theta, so under the sign-free default every
#      refit began by destroying the sign of the effect being studied.
#
# ============================================================
# NMF-RE: Semi-NMF with Random Effects
# Multiplier (Wild) Bootstrap + 1-step Newton
# - Function name: nmf.re
# - dfU cap: cap = dfU_cap_rate * (N*Q) unless dfU_cap (absolute) is supplied
# - Wild weights: exp_centered only  (w ~ Exp(1) - 1; mean 0, var 1)
#
# [PROTOTYPE 2026-06-27 — COMPLETE-EM VARIANT]  derived from
#   nmfre_core20260627_emsigma.R (which already uses the proper EM M-step for
#   sigma^2).  This variant ALSO makes the X-step a genuine EM M-step by adding
#   the posterior-variance term to the basis subproblem:
#       X-step:  min_{X>=0}  ||Y - X(Theta A + Uhat)||_F^2 + tr( X S X' ),
#       S = sum_n Cov(u_n | y_n) = N * sigma^2 * (X'X + lambda I)^{-1}  (E-step).
#   The gradient adds +2 X S, so the Ding-Li-Jordan multiplicative update simply
#   replaces B B' by G = B B' + S (S symmetric PD => monotone for the augmented
#   objective).  With this change the X, Theta, sigma^2, tau^2 steps are ALL exact
#   EM/ECM steps for the marginal Gaussian mixed model, i.e. a fully exact EM
#   (controlled by x_postvar = TRUE; set FALSE to recover the emsigma behaviour).
#
# [REVISION 2026-06-23]  (Psychometrika PSY-2026-0040, Major Revision)
# Revised from nmfre_core20260112.R. Key methodological changes:
#  - Theta is now SIGN-FREE by default (theta_sign = "free"): the essential
#    NMF constraint is X >= 0 (parts-based non-negative basis); the latent
#    scores B = Theta A + U are real-valued. This makes NMF-RE a *semi-NMF*
#    with random effects, and lets Theta be interpreted as covariate EFFECTS
#    (positive, zero, or negative). theta_sign = "nonneg" recovers the old
#    non-negative-Theta behaviour for compositional/intensity data.
#  - Theta-step (free): exact closed-form least squares
#      Theta = (X'X)^{-1} X' (Y - X U) A' (A A')^{-1},
#    solving the stated subproblem exactly (no positive-part surrogate).
#  - X-step: Ding-Li-Jordan (2010) semi-NMF multiplicative update, which is
#    monotone-decreasing for X >= 0 with a sign-free score matrix B, replacing
#    the previous positive-part B_pos = pmax(ThetaA + U, 0) stabilization.
#    x_update = "mu" recovers the old positive-part multiplicative update.
#  - Inference: for sign-free Theta the parameter is interior, so the default
#    p-value is TWO-sided and bootstrap replicates are NOT projected to >= 0.
#    For theta_sign = "nonneg" the boundary null (one-sided, project >= 0) is
#    retained, matching the submitted manuscript.
#
# [ALIGN 2026-01-11] (inherited)
# - X normalization rescales BOTH Theta and U (identifiability-consistent)
# - Bootstrap score weights centered: w = rexp(N) - 1
# - Only "exp_centered" wild weights remain
# ============================================================

pmax_eps <- function(x, eps = 1e-10) pmax(x, eps)

safe_div <- function(num, den, eps = 1e-10) {
  res <- num / pmax_eps(den, eps)
  res[is.na(res) | is.infinite(res)] <- 0
  res
}

# ------------------------------------------------------------
# Semi-NMF multiplicative update for X >= 0 with a sign-free
# score matrix B (Ding, Li & Jordan, 2010, IEEE TPAMI 32(1):45-55).
# Solves (monotonically) min_{X>=0} ||Y - X B||_F^2 with B real-valued:
#   X <- X * sqrt( [ (Y B')+ + X (B B')- ] / [ (Y B')- + X (B B')+ ] ).
# pos()/neg() are the elementwise positive/negative parts.
# ------------------------------------------------------------
.seminmf_update_X <- function(X, Y, B, eps = 1e-10) {
  YBt <- Y %*% t(B)          # P x Q
  BBt <- B %*% t(B)          # Q x Q
  pos <- function(M) pmax(M, 0)
  neg <- function(M) pmax(-M, 0)
  numer <- pos(YBt) + X %*% neg(BBt)
  denom <- neg(YBt) + X %*% pos(BBt)
  X * sqrt(safe_div(numer, denom, eps = eps))
}

# ------------------------------------------------------------
# COMPLETE-EM penalized semi-NMF update for X >= 0 (2026-06-27 prototype):
# Solves (monotonically) min_{X>=0} ||Y - X B||_F^2 + tr(X S X'), where S is a
# symmetric positive-definite matrix (here S = sum_n Cov(u_n|y_n), the E-step
# posterior-variance contribution).  The penalty gradient adds +2 X S, so the
# Ding-Li-Jordan update replaces B B' by G = B B' + S:
#   X <- X * sqrt( [ (Y B')+ + X G- ] / [ (Y B')- + X G+ ] ),  G = B B' + S.
# Monotonicity follows from the same auxiliary-function argument (the pos/neg
# split applies to any symmetric matrix in the quadratic form tr(X G X')).
# ------------------------------------------------------------
.seminmf_update_X_pen <- function(X, Y, B, S, eps = 1e-10) {
  YBt <- Y %*% t(B)              # P x Q
  G   <- B %*% t(B) + S          # Q x Q symmetric
  pos <- function(M) pmax(M, 0)
  neg <- function(M) pmax(-M, 0)
  numer <- pos(YBt) + X %*% neg(G)
  denom <- neg(YBt) + X %*% pos(G)
  X * sqrt(safe_div(numer, denom, eps = eps))
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
# dfU helper (diagnostic only)
# dfU(lambda) = N * sum_i d_i/(d_i + lambda), where d_i are eigenvalues of XtX
# ============================================================

.dfU_from_lambda <- function(d, N, lambda) {
  N * sum(d / (d + lambda))
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

  # ---- sign constraint on Theta and X-step variant ----
  # theta_sign = "free"   : real-valued Theta (semi-NMF; covariate effects, default)
  # theta_sign = "nonneg" : Theta >= 0 (compositional/intensity; submitted-paper behaviour)
  theta_sign = c("free", "nonneg"),
  # x_update = "seminmf"  : Ding et al. (2010) monotone semi-NMF update (default)
  # x_update = "mu"       : legacy positive-part multiplicative update
  x_update = c("seminmf", "mu"),
  # x_postvar = TRUE: COMPLETE-EM X-step (add posterior-variance term tr(X S X'),
  #   S = N*sigma2*(X'X+lambda I)^{-1}); FALSE recovers the conditional X-step.
  x_postvar = TRUE,

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

  # ---- dfU control (DEPRECATED: no cap is applied) ----
  # df_U is no longer capped: the variance components (sigma2, tau2) are
  # estimated by ECM and df_U/(NQ) is reported only as a complexity diagnostic.
  # These arguments are kept (inert) for backward compatibility with existing
  # call sites; "cap" is accepted but behaves identically to "off".
  dfU_control = c("off", "cap"),
  dfU_cap = NULL,           # (inert) retained for backward compatibility
  dfU_cap_rate = 0.10,      # (inert) retained for backward compatibility
  dfU_enforce_every = 1,    # (inert)
  dfU_lambda_max = 1e12,    # (inert)

  # ---- optimization ----
  # Outer-inner structure (revised 2026-06-26): the variance components
  # (sigma2, tau2) and the dfU-cap are updated only on the OUTER loop; between
  # outer updates lambda = sigma2/tau2 is held FIXED and the block-coordinate
  # updates (U, X, Theta) are iterated to convergence on the INNER loop. Each
  # inner phase therefore minimizes a single fixed-lambda objective and is
  # monotone nonincreasing. maxit is the TOTAL inner-iteration budget.
  maxit = 200000,            # TOTAL inner-iteration budget (backstop)
  epsilon = 1e-8,            # inner (fixed-lambda) relative-change tolerance
  inner_maxit = 10000,       # per-outer cap on inner iterations
  outer_maxit = 500,         # max number of variance/lambda updates (generous:
                             # the required count is data-dependent, so set it
                             # high enough that non-convergence means "needs more",
                             # not an artificial cutoff; e.g. HS1939 needs ~112)
  epsilon_outer = 1e-6,      # outer convergence: relative change in lambda
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

  # ---- p-value side ----
  # NULL (default): resolved by theta_sign -- "two.sided" for free Theta
  # (interior parameter), "one.sided" for nonneg (boundary null H0: Theta=0).
  theta_p_side = NULL,

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
  theta_sign <- match.arg(theta_sign)
  x_update <- match.arg(x_update)

  # Resolve p-value side from theta_sign unless the user set it explicitly.
  if (is.null(theta_p_side)) {
    theta_p_side <- if (theta_sign == "nonneg") "one.sided" else "two.sided"
  }
  theta_p_side <- match.arg(theta_p_side, c("one.sided", "two.sided"))

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
  ## 2026-07-18: clip the initial Theta only when it is constrained to be non-negative.
  Theta_ini0 <- if (theta_sign == "nonneg") pmax(Theta, .eps) else Theta
  normed <- normalize_X(pmax(X, .eps), Theta_ini0, U)
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
  nll_trace <- if (isTRUE(save_obj_trace)) rep(NA_real_, maxit) else NULL  # marginal negative log-likelihood (ECM objective)

  # ---- dfU diagnostics ----
  dfU_last <- NA_real_
  dfU_cap_last <- NA_real_
  lambda_last <- NA_real_

  dfU_hit_cap <- FALSE
  dfU_hit_iter <- integer(0)

  # =========================================================
  # main loop
  # =========================================================
  # Does lambda move at all across outer steps? (variance updates or a cap whose
  # floor tracks X). If not, a single outer pass = pure fixed-lambda descent.
  any_var_update <- isTRUE(sigma2_update) || isTRUE(tau2_update)

  total_iter <- 0L
  lambda <- sigma2 / tau2
  outer_break_iters <- integer(0)   # cumulative inner-iteration index at each outer boundary
  lambda_path <- numeric(0)         # lambda used in each outer (inner) phase

  for (outer in 1:outer_maxit) {

    # ----- fix lambda for this inner solve (variance held constant) -----
    tau2   <- clip_val(tau2,   tau2_min,   tau2_max)
    sigma2 <- clip_val(sigma2, sigma2_min, sigma2_max)
    lambda <- sigma2 / tau2

    lambda_used <- lambda
    lambda_last <- lambda

    # ----- INNER: fixed-lambda block-coordinate descent (monotone) -----
    obj_prev_inner <- Inf
    inner_done <- 0L
    for (inner in 1:inner_maxit) {
      total_iter <- total_iter + 1L
      iter_done  <- total_iter
      inner_done <- inner

      ThetaA <- Theta %*% A

      # (1) U-step: ridge / BLUP at the fixed lambda
      XtX   <- crossprod(X)
      M     <- XtX + diag(pmax(lambda, 1e-12), Q)
      cholM <- tryCatch(chol(M), error = function(e) NULL)
      if (is.null(cholM)) stop("Cholesky failed: XtX + lambda I not SPD. Increase lambda (or decrease tau2).")
      R0 <- Y - X %*% ThetaA
      for (n in 1:N) {
        rhs <- crossprod(X, R0[, n, drop = FALSE])
        U[, n] <- as.numeric(backsolve(cholM, forwardsolve(t(cholM), rhs)))
      }
      ## NOTE (2026-07-18): the row-centering of U was removed here; the fixed point
      ## satisfies U A' = 0, which is the condition that identifies U against Theta.
      ## See the header.

      # (2) X-step: X >= 0 with sign-free score matrix B = ThetaA + U
      B_sem <- ThetaA + U
      if (x_update == "seminmf") {
        if (isTRUE(x_postvar)) {
          ## COMPLETE-EM X-step: S = sum_n Cov(u_n|y_n) = N*sigma2*(X'X+lambda I)^{-1}.
          ## cholM = chol(X'X + lambda I) was just formed for the U-step (same X).
          S_pv <- N * clip_val(sigma2, sigma2_min, sigma2_max) * chol2inv(cholM)
          X <- .seminmf_update_X_pen(X, Y, B_sem, S_pv, eps = .eps)
        } else {
          X <- .seminmf_update_X(X, Y, B_sem, eps = .eps)
        }
      } else {
        Y_tilde <- pmax(Y, 0); B_pos <- pmax(B_sem, 0)
        numX <- Y_tilde %*% t(B_pos)
        denX <- X %*% (B_pos %*% t(B_pos))
        X <- X * safe_div(numX, denX, eps = .eps)
      }
      X <- pmax(X, .eps)
      normed <- normalize_X(X, Theta, U)
      X     <- pmax(normed$X, .eps)
      Theta <- if (theta_sign == "nonneg") pmax(normed$Theta, .eps) else normed$Theta
      U     <- normed$U

      # (3) Theta-step
      if (theta_sign == "free") {
        Y_star <- Y - X %*% U
        XtX_t  <- crossprod(X) + diag(1e-10, Q)
        AAt_t  <- tcrossprod(A) + diag(1e-10, K)
        rhs    <- crossprod(X, Y_star) %*% t(A)
        Theta  <- solve(XtX_t, rhs)
        Theta  <- t(solve(AAt_t, t(Theta)))
      } else {
        Y_tilde <- pmax(Y, 0); Y_star <- Y_tilde - X %*% U; Y_star_pos <- pmax(Y_star, 0)
        numT <- crossprod(X, Y_star_pos) %*% t(A)
        denT <- (crossprod(X, X) %*% Theta) %*% (A %*% t(A))
        Theta <- Theta * safe_div(numT, denT, eps = .eps)
        Theta <- pmax(Theta, .eps)
      }

      # (4) objective at the FIXED lambda (this is what the inner loop minimizes)
      ThetaA <- Theta %*% A
      R   <- Y - X %*% (ThetaA + U)
      obj <- sum(R^2) + lambda * sum(U^2)
      rss <- sum(R^2)
      if (isTRUE(save_obj_trace) && total_iter <= maxit) {
        obj_trace[total_iter] <- obj
        rss_trace[total_iter] <- rss
        # marginal negative log-likelihood: y_n ~ N(X Theta a_n, sigma2 I + tau2 X X^T)
        Rm   <- Y - X %*% ThetaA                      # P x N fixed-effects residual
        XtRm <- crossprod(X, Rm)                      # Q x N
        Mll  <- crossprod(X) + diag(pmax(lambda, 1e-12), Q)
        quad <- (sum(Rm^2) - sum(XtRm * solve(Mll, XtRm))) / clip_val(sigma2, sigma2_min, sigma2_max)
        logdetQ <- as.numeric(determinant(diag(Q) + (tau2 / clip_val(sigma2, sigma2_min, sigma2_max)) * crossprod(X),
                                          logarithm = TRUE)$modulus)
        nll_trace[total_iter] <- 0.5 * (N * (P * log(2 * pi) + P * log(clip_val(sigma2, sigma2_min, sigma2_max)) + logdetQ) + quad)
      }

      if (!is.finite(obj)) { stop_reason <- "nonfinite_obj"; converged <- FALSE; break }

      rel_change <- if (is.finite(obj_prev_inner))
        abs(obj_prev_inner - obj) / (abs(obj_prev_inner) + .eps) else NA_real_
      if (is.finite(obj_prev_inner) && rel_change < epsilon) break   # inner converged
      obj_prev_inner <- obj
      if (total_iter >= maxit) break                                  # total budget exhausted
    }

    outer_break_iters <- c(outer_break_iters, total_iter)
    lambda_path <- c(lambda_path, lambda_used)

    if (identical(stop_reason, "nonfinite_obj")) break

    # ----- OUTER: variance-component update from the converged inner fit -----
    ThetaA <- Theta %*% A
    R <- Y - X %*% (ThetaA + U)

    # ===== EM-SIGMA VARIANT (2026-06-27) =====================================
    # PROPER JOINT EM M-STEP for (sigma2, tau2) of the marginal mixed model
    # y_n ~ N(X Theta a_n, sigma2 I + tau2 X X^T), treating u_n as the latent
    # variable with posterior Cov(u_n | y_n) = sigma2 (X'X + lambda I)^{-1}.
    # The inner loop has converged at the FIXED lambda = sigma2_old/tau2_old, so
    # sigma2/tau2 here are the E-step values; both M-steps use that sigma2_old:
    #   sigma2_new = ||R||_F^2/(PN) + sigma2_old * dfU/(PN),  dfU = N*tr(H_lambda)
    #   tau2_new   = ||U||_F^2/(NQ) + (sigma2_old/Q)*tr((X'X+lambda I)^{-1})
    # with tr(H_lambda) = Q - lambda*tr((X'X+lambda I)^{-1}).  This is the ONLY
    # change from nmfre_core20260626.R: the original sigma2-step used the naive
    # residual mean square ||R||_F^2/(PN) (it OMITTED the +sigma2*dfU/(PN) term),
    # so it was a residual-moment update, not the EM M-step.  At the fixed point
    # this sigma2-step gives sigma2* = ||R||_F^2/(PN - dfU) (df-adjusted).
    sigma2_old <- sigma2
    Minv <- NULL
    if (isTRUE(sigma2_update) || isTRUE(tau2_update)) {
      Minv <- tryCatch(chol2inv(chol(crossprod(X) + diag(pmax(lambda, 1e-12), Q))),
                       error = function(e) NULL)
    }
    if (isTRUE(sigma2_update)) {
      if (!is.null(Minv)) {
        trH   <- Q - pmax(lambda, 1e-12) * sum(diag(Minv))   # tr(H_lambda) in [0,Q]
        dfU_e <- N * trH
        sigma2 <- clip_val(mean(R^2) + sigma2_old * dfU_e / (P * N), sigma2_min, sigma2_max)
      } else {
        sigma2 <- clip_val(mean(R^2), sigma2_min, sigma2_max)   # fallback: naive
      }
    }
    if (isTRUE(tau2_update)) {
      term1 <- sum(U^2) / (N * Q)
      if (!is.null(Minv)) {
        tau2 <- clip_val(term1 + sigma2_old * sum(diag(Minv)) / Q, tau2_min, tau2_max)
      }
    }
    # ===== end EM-SIGMA VARIANT ==============================================

    # lambda after the variance update (pure ECM; no cap)
    lambda_new <- clip_val(sigma2, sigma2_min, sigma2_max) / clip_val(tau2, tau2_min, tau2_max)

    if (verbose) {
      cat(sprintf("[outer %d] inner=%d obj=%.6g sigma2=%.4g tau2=%.4g lambda %.4g -> %.4g\n",
                  outer, inner_done, obj, sigma2, tau2, lambda_used, lambda_new))
    }

    # ----- outer convergence -----
    if (!any_var_update) { converged <- TRUE; stop_reason <- "fixed_lambda"; break }
    if (abs(lambda_new - lambda_used) / (abs(lambda_used) + .eps) < epsilon_outer) {
      converged <- TRUE; stop_reason <- "outer_lambda"; break
    }
    if (total_iter >= maxit) { stop_reason <- "maxit"; converged <- FALSE; break }
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

  dfU_cap_last <- NA_real_   # no cap is applied
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

      # project to Theta>=0 only under the non-negative (boundary) variant;
      # for sign-free Theta the parameter is interior, so no projection.
      if (theta_sign == "nonneg") theta_b <- pmax(theta_b, 0)

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
    nll_trace = if (isTRUE(save_obj_trace)) nll_trace[seq_len(iter_done)] else NULL,
    outer_break_iters = outer_break_iters,
    lambda_path = lambda_path,

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

    # record inference / method options
    theta_p_side = theta_p_side,
    theta_sign = theta_sign,
    x_update = x_update
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


# (scan_dfU_cap_rate and the .lambda_for_dfU_cap helper were removed together
#  with the dfU-cap procedure; df_U is now reported as a diagnostic only.)
