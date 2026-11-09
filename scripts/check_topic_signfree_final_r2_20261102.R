# Topic example (inaugural addresses, P = 67 words, N = 59, era dummies without an intercept, Q = 3): the non-negative
# variant (C >= 0) against the sign-free fit with the final nmfkc (develop 55ba7c9), for the main-text paragraph that
# introduces the non-negative variant (work memo A48, A83 (2); NEW file, 2026-11-02, C94).
# C65 (b46c55c) compared the sign-free fit (NLL 9879.71) with the old multiplicative non-negative fit (9916.24;
# 2 dNLL = 73), and C70/C81 showed that the exact non-negative fit reaches 9882.47.  Here, with one version:
# (1) both fits from the same starts (mix2, seeds 1-3; best by marginal NLL), 2 x NLL difference;
# (2) era effects of each fit (sign-free: geomin orientation, matched to the non-negative components by X congruence);
# (3) identification-region width of each fit;
# (4) document scores B = C A + U of the non-negative fit: how many entries are negative (U is sign-free) and the minimum.
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript check_topic_signfree_final_r2_20261102.R <cores> <lib>
# Output: check_topic_signfree_final_r2_20261102.{log,rds} (no data matrices in the rds; the corpus is public anyway)
args <- commandArgs(trailingOnly = TRUE); cores <- as.integer(args[1]); lib <- args[2]
suppressMessages({library(nmfkc, lib.loc = lib); library(parallel)})
nmfre <- nmfkc::nmfre; nll_of <- nmfkc:::.nmfre.marginal.nll
cat(sprintf("nmfkc %s from %s | R %s.%s | %s\n", as.character(packageVersion("nmfkc", lib.loc = lib)), lib,
            R.version$major, R.version$minor, format(Sys.time())))
tp <- readRDS("topic_data_20261031.rds"); Y <- tp$Y; A <- tp$A
cat(sprintf("Y %d x %d | A rows: %s\n", nrow(Y), ncol(Y), paste(rownames(A), collapse = ", ")))
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
options(width = 200)

tasks <- expand.grid(signed = c(FALSE, TRUE), seed = 1:3)
fits <- mclapply(seq_len(nrow(tasks)), function(i) {
  s <- tasks$signed[i]
  f <- tryCatch(quiet(nmfre(Y, A, rank = 3, C.signed = s, epsilon = 1e-7, maxit = 5e5, seed = tasks$seed[i],
                            orient = if (s) "geomin" else "none", print.trace = FALSE)), error = function(e) conditionMessage(e))
  list(signed = s, seed = tasks$seed[i], f = f)
}, mc.cores = min(cores, nrow(tasks)), mc.preschedule = FALSE)

pick <- function(s) {
  zz <- Filter(function(z) z$signed == s && !is.character(z$f), fits)
  nl <- sapply(zz, function(z) nll_of(Y, z$f$X, z$f$C, A, z$f$sigma2, z$f$tau2))
  cat(sprintf("%s: NLL by seed %s | degenerate %s\n", if (s) "sign-free" else "non-negative",
              paste(sprintf("%.4f", nl), collapse = " "), paste(sapply(zz, function(z) isTRUE(z$f$degeneracy$flag)), collapse = " ")))
  list(f = zz[[which.min(nl)]]$f, nll = min(nl))
}
nn <- pick(FALSE); sf <- pick(TRUE)
cat(sprintf("\n2 x (NLL non-negative - NLL sign-free) = %.3f\n", 2 * (nn$nll - sf$nll)))

## match the sign-free components to the non-negative ones by congruence of the X columns
cg <- crossprod(nn$f$X / rep(sqrt(colSums(nn$f$X^2)), each = nrow(Y)), sf$f$X / rep(sqrt(colSums(sf$f$X^2)), each = nrow(Y)))
perm <- apply(cg, 1, which.max)
stopifnot(length(unique(perm)) == 3)
cat("congruence of X columns (non-negative component -> sign-free component):",
    paste(sprintf("%d->%d %.3f", 1:3, perm, cg[cbind(1:3, perm)]), collapse = ", "), "\n")
Cn <- nn$f$C; Cs <- sf$f$C[perm, , drop = FALSE]
dimnames(Cn) <- dimnames(Cs) <- list(paste0("Topic", 1:3), rownames(A))
cat("\nera effects, non-negative variant:\n"); print(round(Cn, 2))
cat("era effects, sign-free (components matched):\n"); print(round(Cs, 2))
cat(sprintf("max |X_nn - X_sf| (matched) %.4f | max |Xi_nn - Xi_sf| %.3f\n", max(abs(nn$f$X - sf$f$X[, perm])),
            max(abs(nn$f$X %*% nn$f$C - sf$f$X %*% sf$f$C))))
cat(sprintf("identification width: non-negative %s, sign-free %s\n", format(signif(nn$f$identification$width, 3)),
            format(signif(sf$f$identification$width, 3))))
cat(sprintf("zeros in the non-negative C: %d (entries %s)\n", sum(Cn == 0),
            paste(which(Cn == 0, arr.ind = FALSE), collapse = ",")))

## document scores of the non-negative fit
B <- nn$f$C %*% A + nn$f$U
cat(sprintf("\nnon-negative fit: document scores B = C A + U: negative entries %d of %d (documents with any: %d of %d), min %.2f\n",
            sum(B < 0), length(B), sum(colSums(B < 0) > 0), ncol(B), min(B)))
cat(sprintf("variances: non-negative sigma2 %.4g tau2 %.5g | sign-free sigma2 %.4g tau2 %.5g\n",
            nn$f$sigma2, nn$f$tau2, sf$f$sigma2, sf$f$tau2))

saveRDS(list(nll = c(nonneg = nn$nll, signfree = sf$nll), d2nll = 2 * (nn$nll - sf$nll), C_nonneg = Cn, C_signfree = Cs,
             perm = perm, congruence = cg, width = c(nonneg = nn$f$identification$width, signfree = sf$f$identification$width),
             B_negative = c(n = sum(B < 0), of = length(B), docs = sum(colSums(B < 0) > 0), min = min(B)),
             version = as.character(packageVersion("nmfkc", lib.loc = lib)), lib = lib, time = format(Sys.time())),
        "check_topic_signfree_final_r2_20261102.rds")
cat(sprintf("\nsaved check_topic_signfree_final_r2_20261102.rds | %s\n", format(Sys.time())))
