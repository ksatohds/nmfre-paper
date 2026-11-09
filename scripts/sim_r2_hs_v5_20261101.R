# ----------------------------------------------------------------------------
# HS1939 simulation v5 (NEW file, 2026-11-01; work memo C79) = v4 (sim_r2_hs_v4_20261101.R, C72) with one change:
# the fit is the package default of nmfkc 1.0.3 develop a343c39 -- the start rule of v4 (default and SPA starts,
# default restarts while degenerate, at most 20) now inside nmfre() as X.init = "mix2", and the variances started
# from the residual covariance (probabilistic PCA) instead of sigma2 = tau2 = 1.  C75-C78: from sigma2 = tau2 = 1
# all data-driven starts fell into a degenerate fit in 509 of the 2800 v4 samples; from the new start none did.
# Author: "はい、その順で進めて、作業メモに記録して．".  Everything else (truths, seeds, inference, ESEM, anchor arm,
# outputs) as v4; the v4 description follows.
# HS1939 simulation v4 (NEW file, 2026-11-01; work memo C72).  Replaces v3 (C54), the anchor-free runs (C50, C60)
# and the interval comparison (C63) in one design.  Author: "数値実験の真のΘに0を置きませんか?" -> "exactを既定にして
# commit push して。" -> "それならワイルド区間と尤度プロファイルだけにしましょうか" -> "はい　その提案で良いです　3台の
# サーバーをうまく使って実施して下さい".
# Package: nmfkc 1.0.3, develop commit 041bf1c (X.update = "exact" by default), given as <lib>.
# Truth: the R2 HS truth (sim_r2_hs_truth_20261031.rds: P = 9 tests, N = 300, Q = 3, A = observed covariates with
# intercept) with two covariate effects set to 0 in the canonical representative: age.rev -> Basis1 (0.144 -> 0) and
# sex.2 -> Basis2 (0.050 -> 0).  For mixing e, X(e) = X* M(e) with M(e) = (1-e) I + e/(Q-1)(11' - I) (no anchor for
# e > 0); the canonical representative is the geomin member X_g(e) of the identification region of X(e) (it depends on
# X only), columns aligned to X*, and the data are generated from (X_g(e), Theta0), so the canonical Theta is Theta0 for
# every e while the population identification region widens with e.  e = 0: X_g = X* (anchor x7, region a point).
# Errors (the same seeds as v2-v3 and C63: U first, then a normal matrix Z): gaussian E = sigma Z; varhet E_pn =
# sigma sqrt(w_p) Z_pn, w_p proportional to the unique variances of a 3-factor ML factor analysis of the observed tests
# (mean 1, largest/smallest 4.0).
# Per replicate: fit as A25 (default and SPA starts, seed 1; default starts with seeds 2, 3, ... while degenerate, at most
# 20); geomin orientation (package); identification region; nmfre.inference with the one-step wild bootstrap (B = 1000)
# and the profile-likelihood intervals of the 9 covariate effects; columns aligned to the canonical representative;
# coverage of the canonical Theta and of the whole population region; rejection of theta = 0; Xi Wald intervals and
# ESEM-within-MIMIC (lavaan) for Xi; covariate-level Wald tests; recovery (congruence).  e = 0 with Gaussian errors only:
# the known-anchor arm (X.anchor = x7; wild bootstrap).  No full-refit bootstrap (C2).
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#        Rscript sim_r2_hs_v4_20261101.R <cores> <lib> <outdir> <err>:<eps>:<from>-<to> [...]
# Output: <outdir>/<err>_eps_<e>/rep_XXXX.rds (one per replicate; done replicates are skipped) + stdout log
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1]); lib <- args[2]; outroot <- args[3]; specs <- args[-(1:3)]
suppressMessages({library(nmfkc, lib.loc = lib); library(parallel); library(lavaan)})
nmfre <- nmfkc::nmfre; nmfre.inference <- nmfkc::nmfre.inference; orient_fit <- nmfkc:::.nmfre.orient.fit
source("nmfre_r2_helpers_20261031.R")
tr <- readRDS("sim_r2_hs_truth_20261031.rds")
Xs <- tr$X; Ths <- tr$C; A <- tr$A; s2 <- tr$sigma2; t2 <- tr$tau2
P <- nrow(Xs); Q <- ncol(Xs); N <- ncol(A); K <- nrow(A); z <- qnorm(0.975)
ys <- rownames(tr$Y_obs); cv <- setdiff(rownames(A), "intercept"); kcv <- match(cv, rownames(A))
cvi <- as.vector(outer(1:Q, (kcv - 1) * Q, "+"))
stopifnot(abs(Ths["Basis1", "age.rev"] - 0.1436) < 1e-3, abs(Ths["Basis2", "sex.2"] - 0.0505) < 1e-3)
Th0 <- Ths; Th0["Basis1", "age.rev"] <- 0; Th0["Basis2", "sex.2"] <- 0
fa <- factanal(t(tr$Y_obs), factors = 3); psi <- fa$uniquenesses * apply(tr$Y_obs, 1, var); w <- psi / mean(psi)
Mof <- function(e) (1 - e) * diag(Q) + e / (Q - 1) * (matrix(1, Q, Q) - diag(Q))
cong <- function(X, Y) colSums(X * Y) / sqrt(colSums(X^2) * colSums(Y^2))
pop <- function(X, C) orient_fit(list(X = X, C = C, U = matrix(0, Q, N), sigma2 = s2, tau2 = t2), orient = "none")
truth_of <- function(e) {
  og <- orient_fit(list(X = Xs %*% Mof(e), C = Ths, U = matrix(0, Q, N), sigma2 = s2, tau2 = t2), orient = "geomin")
  p <- r2_align(og$X, Xs); Xg <- og$X[, p]; dimnames(Xg) <- dimnames(Xs)
  re <- pop(Xg, Th0); again <- orient_fit(list(X = Xg, C = Th0, U = matrix(0, Q, N), sigma2 = s2, tau2 = t2), orient = "geomin")
  list(e = e, X = Xg, C = Th0, Xi = Xg %*% Th0, L = re$C.ident.lower, U = re$C.ident.upper,
       arc = re$identification$arc, moved_again = again$identification$orient$moved)
}
tasks <- do.call(rbind, lapply(specs, function(s) {
  pr <- strsplit(s, ":")[[1]]; rr <- as.integer(strsplit(pr[3], "-")[[1]])
  data.frame(err = pr[1], eps = as.numeric(pr[2]), r = rr[1]:rr[2], stringsAsFactors = FALSE)
}))
truths <- lapply(sort(unique(tasks$eps)), truth_of); names(truths) <- sprintf("%.2f", sort(unique(tasks$eps)))
cat(sprintf("nmfkc %s from %s | lavaan %s | R %s.%s | cores %d | OPENBLAS_NUM_THREADS %s | %s\n",
            as.character(packageVersion("nmfkc", lib.loc = lib)), lib, as.character(packageVersion("lavaan")),
            R.version$major, R.version$minor, cores, Sys.getenv("OPENBLAS_NUM_THREADS"), format(Sys.time())))
cat("Theta0 (canonical):\n"); print(round(Th0, 4))
cat("variance weights w_p:", sprintf("%.2f", w), "\n")
for (tt in truths) cat(sprintf("eps %.2f: region arc %.2f deg | width (covariates) max %.4f | geomin again moves %.1e | min congruence X_g vs X* %.4f | zeros in region: %s\n",
  tt$e, if (is.null(tt$arc)) 0 else diff(tt$arc) * 180 / pi, max((tt$U - tt$L)[, kcv]), tt$moved_again,
  min(cong(tt$X, Xs)), paste(sprintf("%s:%s [%.3f, %.3f]", c("age.rev", "sex.2"), c("B1", "B2"),
    c(tt$L["Basis1", "age.rev"], tt$L["Basis2", "sex.2"]), c(tt$U["Basis1", "age.rev"], tt$U["Basis2", "sex.2"])), collapse = "; ")))

lab <- outer(seq_len(P), 1:Q, function(i, j) paste0("l", i, j))
gam <- outer(1:Q, seq_along(cv), function(j, k) paste0("g", j, k))
m0 <- paste0(
  paste0('efa("efa1")*f', 1:Q, " =~ ", sapply(1:Q, function(j) paste0(lab[, j], "*", ys, collapse = " + ")), collapse = "\n"), "\n",
  paste0("f", 1:Q, " ~ ", sapply(1:Q, function(j) paste0(gam[j, ], "*", cv, collapse = " + ")), collapse = "\n"))
dat_obs <- data.frame(t(tr$Y_obs), t(A[cv, , drop = FALSE]))
pt <- parTable(sem(m0, data = dat_obs, rotation = "none", fixed.x = TRUE, meanstructure = TRUE))
fixed0 <- pt$label[pt$op == "=~" & pt$free == 0L & nzchar(pt$label)]
defs <- unlist(lapply(seq_len(P), function(i) lapply(seq_along(cv), function(k) {
  terms <- paste0(lab[i, ], "*", gam[, k]); keep <- !(lab[i, ] %in% fixed0)
  sprintf("xi_%d_%d := %s", i, k, paste(terms[keep], collapse = " + ")) })))
m_esem <- paste(c(m0, defs), collapse = "\n")

## v5: the package default of nmfkc a343c39 -- X.init = "mix2" (k-means and SPA starts, k-means restarts with seeds
## 2, 3, ... while the kept fit is degenerate, at most 20; the rule v4 ran outside the package) with the variances
## started from the residual covariance (probabilistic PCA; v4 started them at 1).  start keeps v4's fields.
fit_adaptive <- function(Y) {
  f <- tryCatch(suppressMessages(suppressWarnings(nmfre(Y, A, rank = Q, epsilon = 1e-7, maxit = 500000,
          seed = 1, print.trace = FALSE))), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- f$X.init.search
  nl <- setNames(s$nll, c("default", "spa", if (nrow(s) > 2) paste0("default", s$seed[-(1:2)])))
  f$start <- list(chosen = names(nl)[s$chosen], nll_default = s$nll[1], nll_spa = s$nll[2],
                  deg_default = isTRUE(s$degenerate[1]), deg_spa = isTRUE(s$degenerate[2]),
                  n_extra = nrow(s) - 2L, nll_all = nl, variance_start = f$variance.start)
  f
}
## coefficient table in the target's labelling: fit row p[q] <-> target row q
coef_table <- function(fo, fi, p, tt, with_profile = TRUE) {
  g <- function(M) if (is.null(M)) rep(NA_real_, length(cvi)) else as.vector(M[p, kcv, drop = FALSE])
  est <- g(fo$C); se <- g(fi$C.se); bl <- g(fi$C.ci.lower); bu <- g(fi$C.ci.upper)
  L <- g(fo$C.ident.lower); U <- g(fo$C.ident.upper)
  pl <- if (with_profile) g(fi$C.profile.lower) else NA; pu <- if (with_profile) g(fi$C.profile.upper) else NA
  tg <- as.vector(tt$C[, kcv]); pL <- as.vector(tt$L[, kcv]); pU <- as.vector(tt$U[, kcv])
  data.frame(coef = paste0(rep(cv, each = Q), ":B", rep(1:Q, length(cv))), target = tg, popL = pL, popU = pU,
             est = est, se = se, L = L, U = U, wild_lo = bl, wild_hi = bu, prof_lo = pl, prof_hi = pu,
             prof_reached = if (with_profile) as.vector(fi$C.profile.reached[p, kcv, drop = FALSE]) else NA,
             cov_wild = bl <= tg & tg <= bu, cov_prof = pl <= tg & tg <= pu,
             set_wild = bl <= pL & pU <= bu, set_prof = pl <= pL & pU <= pu,
             rej0_wild = bl > 0 | bu < 0, rej0_prof = pl > 0 | pu < 0)
}

one_rep <- function(err, e, r) {
  t0 <- proc.time()[3]; tt <- truths[[sprintf("%.2f", e)]]
  set.seed(1000000 + r)
  U <- matrix(rnorm(Q * N, 0, sqrt(t2)), Q, N); Z <- matrix(rnorm(P * N), P, N)
  E <- if (err == "gaussian") sqrt(s2) * Z else sqrt(s2 * w) * Z
  Y <- tt$X %*% (tt$C %*% A + U) + E; rownames(Y) <- ys
  f <- fit_adaptive(Y)
  if (is.null(f)) return(list(err = err, eps = e, r = r, ok = FALSE, msg = "nmfre failed"))
  t_fit <- proc.time()[3] - t0
  fo <- orient_fit(f, nonneg = FALSE, orient = "geomin"); fo$C.signed <- TRUE
  t1 <- proc.time()[3]
  fi <- suppressMessages(suppressWarnings(nmfre.inference(fo, Y, A, wild.bootstrap = TRUE, wild.B = 1000, wild.seed = r,
                                                          profile = TRUE, profile.which = cvi)))
  t_inf <- proc.time()[3] - t1
  p <- r2_align(fo$X, tt$X)
  ct <- coef_table(fo, fi, p, tt)
  V <- fi$C.vec.cov; Xa <- fo$X[, p]; Ca <- fo$C[p, , drop = FALSE]
  perm_idx <- as.vector(matrix(seq_len(Q * K), Q, K)[p, , drop = FALSE]); Va <- V[perm_idx, perm_idx]
  se_xi <- sapply(seq_len(K), function(k) { id <- (k - 1) * Q + 1:Q; sqrt(pmax(diag(Xa %*% Va[id, id] %*% t(Xa)), 0)) })
  xi <- Xa %*% Ca; xs_cv <- tt$Xi[, kcv, drop = FALSE]
  wald <- r2_wald_cols(fi)
  t2s <- proc.time()[3]
  dat <- data.frame(t(Y), t(A[cv, , drop = FALSE]))
  fe <- tryCatch(suppressWarnings(sem(m_esem, data = dat, rotation = "none", fixed.x = TRUE, meanstructure = TRUE)),
                 error = function(e) NULL)
  es <- if (!is.null(fe)) {
    pe <- parameterEstimates(fe); pe <- pe[pe$op == ":=", ]
    th_diag <- diag(lavInspect(fe, "est")$theta)[ys]
    list(ok = isTRUE(lavInspect(fe, "converged")), heywood = any(th_diag < 0),
         xi = matrix(pe$est, P, length(cv), byrow = TRUE), se = matrix(pe$se, P, length(cv), byrow = TRUE))
  } else list(ok = FALSE)
  t_esem <- proc.time()[3] - t2s
  ## known-anchor arm: e = 0, Gaussian errors
  anc <- NULL
  if (e == 0 && err == "gaussian") {
    t3 <- proc.time()[3]
    fa <- tryCatch(suppressMessages(suppressWarnings(nmfre(Y, A, rank = Q, X.anchor = 7, epsilon = 1e-7, maxit = 500000,
            seed = 1, print.trace = FALSE))), error = function(e) NULL)
    if (!is.null(fa)) {
      fa$C.signed <- TRUE
      fia <- suppressMessages(suppressWarnings(nmfre.inference(fa, Y, A, wild.bootstrap = TRUE, wild.B = 1000, wild.seed = r)))
      pa <- r2_align(fa$X, tt$X)
      anc <- list(cost = r2_nll(Y, fa, A) - r2_nll(Y, f, A), degenerate = isTRUE(fa$degeneracy$flag),
                  width = fa$identification$width, congruence = cong(fa$X[, pa], tt$X),
                  table = coef_table(fa, fia, pa, tt, with_profile = FALSE), time = proc.time()[3] - t3)
    }
  }
  list(err = err, eps = e, r = r, ok = TRUE, conv = isTRUE(f$converged), X.update = f$X.update,
       sigma2 = f$sigma2, tau2 = f$tau2, start = f$start, degenerate = isTRUE(f$degeneracy$flag),
       width = fo$identification$width, arc = fo$identification$arc, nll = r2_nll(Y, f, A),
       congruence = cong(Xa, tt$X), congruence_Xstar = cong(Xa, Xs), table = ct,
       cov_xi = abs(xi[, kcv] - xs_cv) <= z * se_xi[, kcv], err_xi = xi[, kcv] - xs_cv,
       esem = es, cov_xi_esem = if (isTRUE(es$ok)) abs(es$xi - xs_cv) <= z * es$se else matrix(NA, P, length(cv)),
       err_xi_esem = if (isTRUE(es$ok)) es$xi - xs_cv else matrix(NA, P, length(cv)),
       wald_p = setNames(wald$p_value, wald$covariate), anchor = anc,
       t_fit = t_fit, t_inf = t_inf, t_esem = t_esem, t_total = proc.time()[3] - t0)
}

rds_of <- function(err, e, r) file.path(outroot, sprintf("%s_eps_%.2f", err, e), sprintf("rep_%04d.rds", r))
for (i in seq_len(nrow(tasks))) dir.create(dirname(rds_of(tasks$err[i], tasks$eps[i], 1)), recursive = TRUE, showWarnings = FALSE)
todo <- tasks[!file.exists(mapply(rds_of, tasks$err, tasks$eps, tasks$r)), , drop = FALSE]
cat(sprintf("%d replicates to run (%d already done) | %s\n", nrow(todo), nrow(tasks) - nrow(todo), format(Sys.time())))
t_start <- proc.time()[3]
invisible(mclapply(seq_len(nrow(todo)), function(i) {
  x <- tryCatch(one_rep(todo$err[i], todo$eps[i], todo$r[i]),
                error = function(err) list(err = todo$err[i], eps = todo$eps[i], r = todo$r[i], ok = FALSE, msg = conditionMessage(err)))
  saveRDS(x, rds_of(todo$err[i], todo$eps[i], todo$r[i])); NULL
}, mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("wall %.0f s | %s\n", proc.time()[3] - t_start, format(Sys.time())))
for (g in unique(paste(tasks$err, tasks$eps))) {
  er <- sub(" .*", "", g); e <- as.numeric(sub(".* ", "", g)); rr <- tasks$r[tasks$err == er & tasks$eps == e]
  xs <- lapply(rr, function(r) { f <- rds_of(er, e, r); if (file.exists(f)) readRDS(f) })
  ok <- Filter(function(x) isTRUE(x$ok), xs)
  if (!length(ok)) { cat(sprintf("%s eps %.2f: no replicate ok\n", er, e)); next }
  D <- do.call(rbind, lapply(ok, `[[`, "table"))
  cat(sprintf("%s eps %.2f: %d ok of %d | degenerate %d | cov wild %.3f prof %.3f | zeros: size wild %.3f prof %.3f | time %.0f s\n",
              er, e, length(ok), length(rr), sum(sapply(ok, `[[`, "degenerate")), mean(D$cov_wild), mean(D$cov_prof, na.rm = TRUE),
              mean(D$rej0_wild[D$target == 0]), mean(D$rej0_prof[D$target == 0], na.rm = TRUE), mean(sapply(ok, `[[`, "t_total"))))
  bad <- Filter(function(x) !isTRUE(x$ok), xs); if (length(bad)) cat("  failed:", sapply(bad, `[[`, "r"), unique(sapply(bad, `[[`, "msg")), "\n")
}
cat("\nDONE\n")
