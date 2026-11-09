# OD-HUB (3 hubs, Q = 3), sign-free C with an intercept, under the final nmfkc default (a343c39) and the interval
# procedure of the revision (A46/A51: Wald / one-step wild bootstrap and profile likelihood, no full-refit bootstrap).
# Replaces Part 2 of check_odhub_r2_20261102.R, which ran the full-refit bootstrap of C38 and was stopped (author:
# "waldと尤度プロファイルで十分になったはず").  (NEW file, 2026-11-02; work memo C83.)
# Fits: mix2 with seeds 1-3 (geomin orientation); the best one: degeneracy, identification-region width, C1
# (nmfre.inference, wild.B = 1000, seed 123) and profile-likelihood intervals of the 9 covariate effects (one process
# per coefficient).
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript check_odhub_signfree_r2_20261102.R <cores> <lib>
# Output: check_odhub_signfree_r2_20261102.{log,rds} (estimates and intervals only; no OD data)
args <- commandArgs(trailingOnly = TRUE); cores <- as.integer(args[1]); lib <- args[2]
suppressMessages({library(nmfkc, lib.loc = lib); library(readxl); library(parallel)})
nmfre <- nmfkc::nmfre; nmfre.inference <- nmfkc::nmfre.inference
nll_of <- nmfkc:::.nmfre.marginal.nll
cat(sprintf("nmfkc %s from %s | R %s.%s | cores %d | OPENBLAS_NUM_THREADS %s | %s\n",
            as.character(packageVersion("nmfkc", lib.loc = lib)), lib, R.version$major, R.version$minor, cores,
            Sys.getenv("OPENBLAS_NUM_THREADS"), format(Sys.time())))
dd <- suppressMessages(as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx", sheet = 1, skip = 2)))
dd <- dd[dd[, 1] == "1" & dd[, 3] == "02" & dd[, 5] != "100" & dd[, 7] != "100", ]; pref <- unique(dd[, 6])
M <- matrix(NA, 47, 47); rownames(M) <- colnames(M) <- pref
dd[, 5] <- as.numeric(dd[, 5]); dd[, 7] <- as.numeric(dd[, 7])
for (i in 1:47) for (j in 1:47) M[i, j] <- dd[which(dd[, 5] == i & dd[, 7] == j), 9]
M <- log(1 + M); hb <- c(13, 27, 40); Q <- 3
Y <- M[-hb, ]; A <- rbind(intercept = 1, M[hb, , drop = FALSE]); K <- nrow(A)
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
t0 <- proc.time()[3]
fs <- mclapply(1:3, function(s) {
  t1 <- proc.time()[3]
  f <- tryCatch(quiet(nmfre(Y, A, rank = Q, epsilon = 1e-7, maxit = 5e5, seed = s, orient = "geomin",
                            print.trace = FALSE)), error = function(e) conditionMessage(e))
  list(f = f, time = proc.time()[3] - t1)
}, mc.cores = min(3, cores))
for (s in 1:3) {
  f <- fs[[s]]$f
  if (is.character(f)) { cat(sprintf("seed %d FAILED: %s\n", s, f)); next }
  cat(sprintf("seed %d: NLL %.4f | degenerate %s | converged %s iter %s | %.0f s | starts:\n", s,
              nll_of(Y, f$X, f$C, A, f$sigma2, f$tau2), isTRUE(f$degeneracy$flag), format(f$converged),
              format(f$iter), fs[[s]]$time))
  print(f$X.init.search, row.names = FALSE, digits = 8)
}
ok <- which(!vapply(fs, function(z) is.character(z$f), TRUE))
nl <- sapply(fs[ok], function(z) nll_of(Y, z$f$X, z$f$C, A, z$f$sigma2, z$f$tau2))
fo <- fs[[ok[which.min(nl)]]]$f
cs <- fo$X / rep(sqrt(colSums(fo$X^2)), each = nrow(fo$X))
cat(sprintf("\nbest: NLL %.4f | sigma2 %.4f tau2 %.2f lambda %.5f dfU/(NQ) %.3f | max cosine of X columns %.3f | identification region max width %s\n",
            min(nl), fo$sigma2, fo$tau2, fo$lambda, fo$dfU.frac, max(crossprod(cs)[upper.tri(diag(Q))]),
            if (is.null(fo$C.ident.lower)) "NA" else format(signif(max(fo$C.ident.upper - fo$C.ident.lower), 3))))
cat("X (share of each prefecture by component), rows with share > 0.99 per column:\n")
Xp <- fo$X / rowSums(fo$X)
for (q in 1:Q) cat(sprintf("  Basis%d: %s\n", q, paste(rownames(Y)[Xp[, q] > 0.99], collapse = " ")))
## C1 (wild bootstrap) and profile likelihood, one coefficient per process
pw <- (Q + 1):(Q * K)
res <- mclapply(c(0L, pw), function(j) {
  t1 <- proc.time()[3]
  fi <- tryCatch(quiet(nmfre.inference(fo, Y, A, wild.bootstrap = TRUE, wild.B = 1000, wild.seed = 123,
                                       profile = j > 0, profile.which = if (j > 0) j else NULL)),
                 error = function(e) conditionMessage(e))
  list(j = j, fi = fi, time = proc.time()[3] - t1)
}, mc.cores = cores, mc.preschedule = FALSE)
fi <- res[[1]]$fi
if (is.character(fi)) stop("inference failed: ", fi)
cf <- fi$coefficients
cf$prof_lo <- NA_real_; cf$prof_hi <- NA_real_; cf$prof_reached <- NA; cf$prof_time <- NA_real_
for (z in res[-1]) {
  if (is.character(z$fi)) { cat(sprintf("profile %d FAILED: %s\n", z$j, z$fi)); next }
  cf$prof_lo[z$j] <- as.vector(z$fi$C.profile.lower)[z$j]; cf$prof_hi[z$j] <- as.vector(z$fi$C.profile.upper)[z$j]
  cf$prof_reached[z$j] <- as.vector(z$fi$C.profile.reached)[z$j]; cf$prof_time[z$j] <- z$time
}
cf$wald_lo <- cf$Estimate - qnorm(0.975) * cf$SE; cf$wald_hi <- cf$Estimate + qnorm(0.975) * cf$SE
options(width = 220)
cat("\ncoefficients (C1 = wild-bootstrap percentile interval; Wald = estimate +- 1.96 SE; profile likelihood):\n")
print(cf, digits = 4, row.names = FALSE)
saveRDS(list(fits = lapply(fs, function(z) if (is.character(z$f)) z$f else
                             z$f[intersect(names(z$f), c("X", "C", "sigma2", "tau2", "lambda", "dfU", "dfU.frac",
                                 "converged", "iter", "X.init.search", "variance.start", "degeneracy",
                                 "C.ident.lower", "C.ident.upper"))]),
             nll = nl, coefficients = cf, C.vec.cov = fi$C.vec.cov), "check_odhub_signfree_r2_20261102.rds")
cat(sprintf("\ntotal %.0f s\nDONE %s\n", proc.time()[3] - t0, format(Sys.time())))
