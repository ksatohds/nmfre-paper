# ----------------------------------------------------------------------------
# HS simulation v4 (C72): does the likelihood itself favour the degenerate fits?  (NEW file, 2026-11-01; work memo C75.)
# Author: "真のXから始める確認をサーバーで実行して".  C74: in the replicates whose kept fit stayed degenerate, all 22
# starts (default, SPA, 20 default restarts) reached the same NLL.  Here the same data are fitted once more, started
# at the truth (X = X_g(e), Theta = Theta0, U = the generated U, sigma2, tau2 true), with the package and settings of
# v4 (nmfkc 041bf1c, epsilon 1e-7, maxit 5e5, exact X-step).  Recorded: the NLL at the true parameters, the NLL and
# degeneracy measures of the truth-started fit, and v4's kept NLL (from the v4 replicate file).  Readings:
#   truth start ends degenerate at v4's NLL           -> the likelihood climbs from the truth into the degenerate path
#   truth start ends non-degenerate, NLL below v4's   -> v4's starts missed a better fit (a trap, not the likelihood)
#   truth start ends non-degenerate, NLL above v4's   -> a non-degenerate local optimum exists, the likelihood prefers
#                                                        the degenerate path
# Control: non-degenerate replicates (the first 40 of each design), where the truth start should reach v4's NLL.
# Data generation copied from sim_r2_hs_v4_20261101.R (same truths, seeds 1000000 + r, U first, then Z).
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#        Rscript check_v4_truth_start_20261101.R <cores> <lib> <tasks.csv> <from> <to> <outdir>
# Input: tasks.csv (err, eps, rep, degenerate, nll_v4) from make_tasks below (built locally from the v4 summary).
# Output: <outdir>/<err>_eps_<e>/rep_XXXX.rds (done ones skipped) + stdout log
# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1]); lib <- args[2]; tasks <- read.csv(args[3], stringsAsFactors = FALSE)
from <- as.integer(args[4]); to <- min(as.integer(args[5]), nrow(tasks)); outroot <- args[6]
suppressMessages({library(nmfkc, lib.loc = lib); library(parallel)})
nmfre <- nmfkc::nmfre; orient_fit <- nmfkc:::.nmfre.orient.fit
source("nmfre_r2_helpers_20261031.R")
tr <- readRDS("sim_r2_hs_truth_20261031.rds")
Xs <- tr$X; Ths <- tr$C; A <- tr$A; s2 <- tr$sigma2; t2 <- tr$tau2
P <- nrow(Xs); Q <- ncol(Xs); N <- ncol(A)
ys <- rownames(tr$Y_obs)
Th0 <- Ths; Th0["Basis1", "age.rev"] <- 0; Th0["Basis2", "sex.2"] <- 0
fa <- factanal(t(tr$Y_obs), factors = 3); psi <- fa$uniquenesses * apply(tr$Y_obs, 1, var); w <- psi / mean(psi)
Mof <- function(e) (1 - e) * diag(Q) + e / (Q - 1) * (matrix(1, Q, Q) - diag(Q))
truth_of <- function(e) {
  og <- orient_fit(list(X = Xs %*% Mof(e), C = Ths, U = matrix(0, Q, N), sigma2 = s2, tau2 = t2), orient = "geomin")
  p <- r2_align(og$X, Xs); Xg <- og$X[, p]; dimnames(Xg) <- dimnames(Xs)
  list(e = e, X = Xg, C = Th0)
}
todo <- tasks[from:to, , drop = FALSE]
truths <- lapply(sort(unique(todo$eps)), truth_of); names(truths) <- sprintf("%.2f", sort(unique(todo$eps)))
cat(sprintf("nmfkc %s from %s | R %s.%s | cores %d | OPENBLAS_NUM_THREADS %s | tasks %d-%d of %d | %s\n",
            as.character(packageVersion("nmfkc", lib.loc = lib)), lib, R.version$major, R.version$minor, cores,
            Sys.getenv("OPENBLAS_NUM_THREADS"), from, to, nrow(tasks), format(Sys.time())))
one <- function(err, e, r, degenerate, nll_v4) {
  t0 <- proc.time()[3]; tt <- truths[[sprintf("%.2f", e)]]
  set.seed(1000000 + r)
  U <- matrix(rnorm(Q * N, 0, sqrt(t2)), Q, N); Z <- matrix(rnorm(P * N), P, N)
  E <- if (err == "gaussian") sqrt(s2) * Z else sqrt(s2 * w) * Z
  Y <- tt$X %*% (tt$C %*% A + U) + E; rownames(Y) <- ys
  nll_true <- r2_nll(Y, list(X = tt$X, C = tt$C, sigma2 = s2, tau2 = t2), A)
  f <- suppressMessages(suppressWarnings(nmfre(Y, A, rank = Q, epsilon = 1e-7, maxit = 500000, seed = 1,
         X.init = tt$X, C.init = tt$C, U.init = U, sigma2 = s2, tau2 = t2, print.trace = FALSE)))
  dg <- f$degeneracy
  list(err = err, eps = e, r = r, ok = TRUE, degenerate_v4 = degenerate, nll_v4 = nll_v4, nll_true = nll_true,
       nll_ts = r2_nll(Y, f, A), conv = isTRUE(f$converged), iter = f$iter, degenerate_ts = isTRUE(dg$flag),
       max.cosine = dg$max.cosine, cancellation = dg$cancellation, kappa.X = dg$kappa.X,
       congruence = colSums(f$X * tt$X) / sqrt(colSums(f$X^2) * colSums(tt$X^2)),
       sigma2 = f$sigma2, tau2 = f$tau2, time = proc.time()[3] - t0)
}
rds_of <- function(err, e, r) file.path(outroot, sprintf("%s_eps_%.2f", err, e), sprintf("rep_%04d.rds", r))
for (i in seq_len(nrow(todo))) dir.create(dirname(rds_of(todo$err[i], todo$eps[i], 1)), recursive = TRUE, showWarnings = FALSE)
todo <- todo[!file.exists(mapply(rds_of, todo$err, todo$eps, todo$rep)), , drop = FALSE]
cat(sprintf("%d to run | %s\n", nrow(todo), format(Sys.time())))
t_start <- proc.time()[3]
invisible(mclapply(seq_len(nrow(todo)), function(i) {
  x <- tryCatch(one(todo$err[i], todo$eps[i], todo$rep[i], todo$degenerate[i], todo$nll_v4[i]),
                error = function(err) list(err = todo$err[i], eps = todo$eps[i], r = todo$rep[i], ok = FALSE,
                                           msg = conditionMessage(err)))
  saveRDS(x, rds_of(todo$err[i], todo$eps[i], todo$rep[i])); NULL
}, mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("wall %.0f s | %s\nDONE\n", proc.time()[3] - t_start, format(Sys.time())))
