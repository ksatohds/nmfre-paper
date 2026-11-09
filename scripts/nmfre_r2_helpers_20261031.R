# ----------------------------------------------------------------------------
# Helpers for the R2 revision, built on the output of nmfkc develop 1.0.2 nmfre().
# (NEW file, 2026-10-31.)  The package itself is not modified; these functions take
# an nmfre() fit (a list with X, C, U, sigma2, tau2) and return transformed fits or
# summaries.  Candidates for the package once validated.
#
# Model: Y = X (C A + U) + E,  X >= 0 with 1'X = 1',  U ~ N(0, tau2 I_Q).
# Equal-likelihood family: (X, C, U) -> (X O, O'C, O'U) with O orthogonal, O 1 = 1
# and X O >= 0.  O = 11'/Q + V W V', V an orthonormal basis of 1-perp, W in
# SO(Q-1); (Q-1)(Q-2)/2 free parameters (Q = 2: none; Q = 3: one angle).
#
#   r2_O(par, Q)            rotation in the family
#   r2_rotate(fit, O)       apply it to a fit (X, C, U); Xi = XC and the marginal
#                           likelihood are unchanged
#   r2_orient(fit, rule)    canonical representative: "geomin" (simple structure on
#                           the row shares X / rowSums(X)) or "purity" (largest
#                           smallest-per-basis purity, i.e. closest to anchors)
#   r2_align(X, Xref)       column permutation of X closest to Xref
#   r2_ident_range(fit)     range of each entry of C over the family (Q = 3)
#   r2_wald_cols(fit.inf)   covariate-level joint Wald test H0: C[, k] = 0 (df = Q),
#                           invariant under the family
#   r2_full_refit_boot(...) full-refit wild bootstrap (C2), each replicate oriented
#                           by the same rule and aligned
# ----------------------------------------------------------------------------

r2_perp <- function(Q) {
  V <- qr.Q(qr(cbind(1, diag(Q))))[, 2:Q, drop = FALSE]
  V
}

r2_npar <- function(Q) if (Q <= 2) 0L else as.integer((Q - 1) * (Q - 2) / 2)

## W in SO(Q-1): Q = 3 by the angle; Q >= 4 by the Cayley transform of a skew matrix
r2_W <- function(par, Q) {
  m <- Q - 1
  if (m == 1) return(matrix(1, 1, 1))
  if (m == 2) return(matrix(c(cos(par), sin(par), -sin(par), cos(par)), 2, 2))
  S <- matrix(0, m, m); S[upper.tri(S)] <- par; S <- S - t(S)
  solve(diag(m) - S, diag(m) + S)
}

r2_O <- function(par, Q, V = r2_perp(Q)) {
  matrix(1 / Q, Q, Q) + V %*% r2_W(par, Q) %*% t(V)
}

r2_rotate <- function(fit, O) {
  out <- fit
  out$X <- fit$X %*% O
  out$C <- t(O) %*% fit$C
  if (!is.null(fit$U)) out$U <- t(O) %*% fit$U
  dimnames(out$X) <- dimnames(fit$X); dimnames(out$C) <- dimnames(fit$C)
  if (!is.null(fit$U)) dimnames(out$U) <- dimnames(fit$U)
  out
}

r2_permute <- function(fit, p) {
  out <- fit
  out$X <- fit$X[, p, drop = FALSE]; out$C <- fit$C[p, , drop = FALSE]
  if (!is.null(fit$U)) out$U <- fit$U[p, , drop = FALSE]
  colnames(out$X) <- colnames(fit$X); rownames(out$C) <- rownames(fit$C)
  if (!is.null(fit$U)) rownames(out$U) <- rownames(fit$U)
  out
}

## orientation criteria (smaller is better)
r2_crit <- function(X, rule = c("geomin", "purity"), eps = 0.01) {
  rule <- match.arg(rule)
  S <- X / rowSums(X)
  if (rule == "geomin") sum(apply(S^2 + eps, 1, prod)^(1 / ncol(X)))
  else -min(apply(S, 2, max))
}

r2_orient <- function(fit, rule = c("geomin", "purity"), eps = 0.01, tol = 1e-12,
                      ngrid = 7200, nstart = 10, seed = 1) {
  rule <- match.arg(rule)
  Q <- ncol(fit$X); np <- r2_npar(Q); V <- r2_perp(Q)
  if (np == 0L) return(list(fit = fit, O = diag(Q), par = numeric(0), crit = r2_crit(fit$X, rule, eps),
                            boundary = NA))
  f <- function(par) {
    Xo <- fit$X %*% r2_O(par, Q, V)
    if (min(Xo) < -tol) return(1e10 + sum(pmax(-Xo, 0)))
    r2_crit(Xo, rule, eps)
  }
  if (np == 1L) {
    ph <- seq(-pi, pi, length.out = ngrid + 1)[-1]
    v <- vapply(ph, f, 0)
    j <- which.min(v); h <- ph[2] - ph[1]
    o <- stats::optimize(f, c(ph[j] - h, ph[j] + h), tol = 1e-10)
    par <- if (o$objective <= v[j]) o$minimum else ph[j]
    feas <- v < 1e10
    ## boundary: does the optimum sit at an end of a feasible arc?
    jn <- c(if (j > 1) j - 1 else ngrid, if (j < ngrid) j + 1 else 1)
    boundary <- any(!feas[jn])
  } else {
    set.seed(seed)
    starts <- c(list(rep(0, np)), lapply(seq_len(nstart), function(i) stats::rnorm(np, 0, 0.3)))
    best <- NULL
    for (s0 in starts) {
      if (f(s0) >= 1e10) next
      o <- stats::optim(s0, f, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
      if (is.null(best) || o$value < best$value) best <- o
    }
    par <- best$par; boundary <- NA
  }
  O <- r2_O(par, Q, V)
  list(fit = r2_rotate(fit, O), O = O, par = par, crit = f(par), boundary = boundary)
}

r2_perms <- function(Q) {
  if (Q == 1) return(matrix(1L, 1, 1))
  p <- r2_perms(Q - 1); out <- NULL
  for (i in 1:Q) out <- rbind(out, cbind(i, ifelse(p >= i, p + 1L, p)))
  out
}

r2_align <- function(X, Xref) {
  P <- r2_perms(ncol(Xref)); best <- NULL; bd <- Inf
  for (r in seq_len(nrow(P))) {
    d <- sum((X[, P[r, ]] - Xref)^2)
    if (d < bd) { bd <- d; best <- P[r, ] }
  }
  best
}

## range of each entry of C over the admissible family (Q = 3: exact scan of the arc
## containing the fit)
r2_ident_range <- function(fit, ngrid = 72000, tol = 1e-12) {
  Q <- ncol(fit$X)
  if (Q == 2L) return(list(range = rbind(min = as.vector(fit$C), max = as.vector(fit$C)), arc = c(0, 0)))
  if (Q != 3L) stop("r2_ident_range: implemented for Q = 2, 3 only")
  V <- r2_perp(Q)
  ph <- seq(-pi, pi, length.out = ngrid + 1)
  ok <- vapply(ph, function(p) min(fit$X %*% r2_O(p, Q, V)) >= -tol, TRUE)
  i0 <- which.min(abs(ph)); lo <- i0; hi <- i0
  while (lo > 1 && ok[lo - 1]) lo <- lo - 1
  while (hi < length(ph) && ok[hi + 1]) hi <- hi + 1
  path <- t(vapply(ph[lo:hi], function(p) as.vector(t(r2_O(p, Q, V)) %*% fit$C), numeric(length(fit$C))))
  list(range = rbind(min = apply(path, 2, min), max = apply(path, 2, max)), arc = ph[c(lo, hi)])
}

## covariate-level joint Wald test from nmfre.inference(): vec(C) is column-major
r2_wald_cols <- function(fit.inf) {
  C <- fit.inf$C; V <- fit.inf$C.vec.cov; Q <- nrow(C); K <- ncol(C)
  do.call(rbind, lapply(seq_len(K), function(k) {
    id <- (k - 1) * Q + seq_len(Q); th <- C[, k]
    W <- as.numeric(t(th) %*% solve(V[id, id, drop = FALSE], th))
    data.frame(covariate = colnames(C)[k], W = W, df = Q,
               p_value = stats::pchisq(W, Q, lower.tail = FALSE))
  }))
}

## Marginal NLL (nmfkc internal; works with load_all() and with an installed package)
r2_nll <- function(Y, f, A) get(".nmfre.marginal.nll", envir = asNamespace("nmfkc"))(Y, f$X, f$C, A, f$sigma2, f$tau2)

## Fit from two starts -- nmfre()'s default and X.init = "spa" -- and keep the one with the
## smaller marginal NLL (added 2026-10-31 after C24: a default start was trapped in a
## degenerate solution whose NLL was 62-128 worse; the SPA start reached the better one).
r2_fit_best <- function(Y, A, Q, epsilon = 1e-7, maxit = 500000, seed = 1, ...) {
  fd <- tryCatch(suppressMessages(suppressWarnings(nmfre(Y, A, rank = Q, epsilon = epsilon, maxit = maxit,
          seed = seed, print.trace = FALSE, ...))), error = function(e) NULL)
  fs <- tryCatch(suppressMessages(suppressWarnings(nmfre(Y, A, rank = Q, epsilon = epsilon, maxit = maxit,
          seed = seed, X.init = "spa", print.trace = FALSE, ...))), error = function(e) NULL)
  nd <- if (is.null(fd)) Inf else r2_nll(Y, fd, A); ns <- if (is.null(fs)) Inf else r2_nll(Y, fs, A)
  f <- if (ns < nd) fs else fd
  if (!is.null(f)) f$start <- list(chosen = if (ns < nd) "spa" else "default", nll_default = nd, nll_spa = ns,
                                    deg_default = isTRUE(fd$degeneracy$flag), deg_spa = isTRUE(fs$degeneracy$flag))
  f
}

## C2: full-refit bootstrap, then each replicate is oriented by `rule` and its columns
## aligned to Xref.  Returns the replicate C's (B x QK) and Xi's (B x PK).
## scheme (added 2026-10-31 after C13-C15; the default keeps the validated behaviour):
##   "residual"   : as the paper's full_refit_boot(): Yb = X(CA+U) + w_n * r_n with U held
##                  at its BLUP and r_n = y_n - X(Ca_n+u_n), w ~ Exp(1) - 1.  For Q = 3 these
##                  data pull the refit to a degenerate solution (two bases merged, theta
##                  diverging) with a LOWER marginal NLL (C14-C15), so do not use it for Q > 1.
##   "marginal"   : wild multipliers on the marginal residuals, Yb = XCA + w_n (y_n - XCa_n);
##                  u_n and e_n are scaled together, keeping the shape tau2 XX' + sigma2 I.
##   "parametric" : Yb = X(CA + U*) + E*, U* ~ N(0, tau2 I), E* ~ N(0, sigma2 I).
## Variances are re-estimated on each refit, warm-started from the fit.
## restart.spa (added 2026-10-31): when a warm-started refit is flagged degenerate, refit from
## X.init = "spa" as well and keep the smaller marginal NLL on Yb.
r2_full_refit_boot <- function(fit, Y, A, B = 100, rule = c("geomin", "purity", "none"),
                               Xref = fit$X, seed = 1, epsilon = 1e-6, maxit = 50000,
                               scheme = c("residual", "marginal", "parametric"),
                               restart.spa = FALSE, ...) {
  rule <- match.arg(rule); scheme <- match.arg(scheme)
  Q <- ncol(fit$X); P <- nrow(Y); K <- nrow(A); N <- ncol(Y)
  fitted0 <- fit$X %*% (fit$C %*% A + fit$U); R <- Y - fitted0
  M0 <- fit$X %*% fit$C %*% A; Rm <- Y - M0
  Cb <- matrix(NA_real_, B, Q * K); Xib <- matrix(NA_real_, B, P * K); conv <- rep(NA, B)
  Xb <- matrix(NA_real_, B, P * Q); deg <- rep(NA, B); rst <- rep(FALSE, B)
  for (b in seq_len(B)) {
    set.seed(seed + b)
    Yb <- switch(scheme,
      residual   = fitted0 + sweep(R, 2, stats::rexp(N, 1) - 1, `*`),
      marginal   = M0 + sweep(Rm, 2, stats::rexp(N, 1) - 1, `*`),
      parametric = fit$X %*% (fit$C %*% A + matrix(stats::rnorm(Q * N, 0, sqrt(fit$tau2)), Q, N)) +
                   matrix(stats::rnorm(P * N, 0, sqrt(fit$sigma2)), P, N))
    fb <- tryCatch(suppressMessages(suppressWarnings(
      nmfre(Yb, A, rank = Q, X.init = fit$X, C.init = fit$C, U.init = fit$U,
            sigma2 = fit$sigma2, tau2 = fit$tau2, epsilon = epsilon, maxit = maxit,
            print.trace = FALSE, seed = seed + b, ...))), error = function(e) NULL)
    if (is.null(fb)) next
    deg[b] <- isTRUE(fb$degeneracy$flag)          # nmfkc >= 1.0.2.9000; NA before
    if (restart.spa && deg[b]) {
      fs <- tryCatch(suppressMessages(suppressWarnings(
        nmfre(Yb, A, rank = Q, X.init = "spa", epsilon = epsilon, maxit = maxit, seed = seed + b,
              print.trace = FALSE, ...))), error = function(e) NULL)
      if (!is.null(fs) && r2_nll(Yb, fs, A) < r2_nll(Yb, fb, A)) { fb <- fs; rst[b] <- TRUE }
    }
    if (rule != "none") fb <- r2_orient(fb, rule)$fit
    fb <- r2_permute(fb, r2_align(fb$X, Xref))
    Cb[b, ] <- as.vector(fb$C); Xib[b, ] <- as.vector(fb$X %*% fb$C); conv[b] <- isTRUE(fb$converged)
    Xb[b, ] <- as.vector(fb$X)
  }
  list(C = Cb, Xi = Xib, X = Xb, converged = conv, degenerate = deg, restarted = rst)
}
