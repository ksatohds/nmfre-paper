# OD-HUB (3 hubs, Q = 3), non-negative variant (C >= 0): without vs with an intercept row, under the final nmfkc
# default (a343c39).  For the R2 draft 3.8 (why the appendix fit has no intercept).  This is Part 1 of
# check_odhub_r2_20261102.R, which was stopped in its Part 2 before it saved an RDS; the numbers here are the same
# computation, now with the RDS.  (NEW file, 2026-11-02; work memo C83.)
# Fits: mix2 with seeds 1-5 each; best fit: nmfre.inference (wild.B = 1000, seed 123, one-sided p).
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript check_odhub_intercept_r2_20261102.R <cores> <lib>
# Output: check_odhub_intercept_r2_20261102.{log,rds} (estimates and inference only; no OD data)
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
Y <- M[-hb, ]; A0 <- M[hb, , drop = FALSE]; A1 <- rbind(intercept = 1, A0)
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
keep <- c("X", "C", "sigma2", "tau2", "lambda", "dfU", "dfU.frac", "converged", "iter", "X.init.search",
          "variance.start", "degeneracy")
out <- list()
options(width = 200)
for (nm in c("noint", "int")) {
  A <- if (nm == "noint") A0 else A1
  fs <- mclapply(1:5, function(s) tryCatch(quiet(nmfre(Y, A, rank = Q, C.signed = FALSE, epsilon = 1e-7, maxit = 5e5,
          seed = s, print.trace = FALSE)), error = function(e) conditionMessage(e)), mc.cores = min(5, cores))
  nl <- sapply(fs, function(f) if (is.character(f)) NA else nll_of(Y, f$X, f$C, A, f$sigma2, f$tau2))
  dg <- sapply(fs, function(f) if (is.character(f)) NA else isTRUE(f$degeneracy$flag))
  b <- fs[[which.min(nl)]]
  fi <- quiet(nmfre.inference(b, Y, A, wild.bootstrap = TRUE, wild.B = 1000, wild.seed = 123))
  cat(sprintf("\n-- %s: NLL by seed %s | degenerate %s | sigma2 %.4f tau2 %.2f dfU/(NQ) %.3f\n", nm,
              paste(sprintf("%.4f", nl), collapse = " "), paste(dg, collapse = " "), b$sigma2, b$tau2, b$dfU.frac))
  print(fi$coefficients[, c("Basis", "Covariate", "Estimate", "SE", "BSE", "z_value", "p_value")], digits = 5,
        row.names = FALSE)
  out[[nm]] <- list(nll = nl, degenerate = dg, fit = b[intersect(names(b), keep)], coefficients = fi$coefficients,
                    C.vec.cov = fi$C.vec.cov)
}
d <- min(out$int$nll) - min(out$noint$nll)
cat(sprintf("\nintercept vs none: NLL %.4f vs %.4f, change %+.4f (2dNLL %+.3f; 3 parameters added)\n",
            min(out$int$nll), min(out$noint$nll), d, 2 * d))
out$d2nll <- 2 * d
saveRDS(out, "check_odhub_intercept_r2_20261102.rds")
cat(sprintf("DONE %s\n", format(Sys.time())))
