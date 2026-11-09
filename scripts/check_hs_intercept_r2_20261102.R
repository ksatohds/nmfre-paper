# Holzinger-Swineford: with vs without an intercept row under the final nmfkc default (a343c39), for the sentence in
# Section 6.2 "twice the marginal log-likelihood is larger by 355 with the three intercepts" (that 355 was C10, nmfkc
# 1.0.2.9000).  mix2 fits with seeds 1-3 for each specification.  (NEW file, 2026-11-02; work memo C89.)
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript check_hs_intercept_r2_20261102.R <cores> <lib>
# Output: check_hs_intercept_r2_20261102.{log,rds}
args <- commandArgs(trailingOnly = TRUE); cores <- as.integer(args[1]); lib <- args[2]
suppressMessages({library(nmfkc, lib.loc = lib); library(parallel)})
nll_of <- nmfkc:::.nmfre.marginal.nll
cat(sprintf("nmfkc %s from %s | R %s.%s | %s\n", as.character(packageVersion("nmfkc", lib.loc = lib)), lib,
            R.version$major, R.version$minor, format(Sys.time())))
data(HolzingerSwineford1939, package = "lavaan")
d <- HolzingerSwineford1939; d <- d[complete.cases(d), ]
d$age.rev <- -(d$ageyr + d$agemo / 12); d$sex.2 <- ifelse(d$sex == 2, 1, 0)
d$school.GW <- ifelse(d$school == "Grant-White", 1, 0)
d[, c("id", "sex", "ageyr", "agemo", "school", "grade")] <- NULL; d <- nmfkc.normalize(d)
exh <- c("age.rev", "sex.2", "school.GW"); enh <- setdiff(colnames(d), exh)
Y <- t(as.matrix(d[, enh])); A0 <- t(as.matrix(d[, exh])); A1 <- rbind(intercept = 1, A0)
tasks <- expand.grid(spec = c("int", "noint"), seed = 1:3, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(tasks)), function(i) {
  A <- if (tasks$spec[i] == "int") A1 else A0
  f <- suppressMessages(suppressWarnings(nmfre(Y, A, rank = 3, epsilon = 1e-7, maxit = 5e5, seed = tasks$seed[i],
                                               print.trace = FALSE)))
  data.frame(spec = tasks$spec[i], seed = tasks$seed[i], nll = nll_of(Y, f$X, f$C, A, f$sigma2, f$tau2),
             degenerate = isTRUE(f$degeneracy$flag), converged = isTRUE(f$converged), sigma2 = f$sigma2, tau2 = f$tau2)
}, mc.cores = cores)
R <- do.call(rbind, res); print(R, row.names = FALSE, digits = 8)
b <- tapply(R$nll, R$spec, min)
cat(sprintf("\nbest NLL: intercept %.4f, none %.4f | 2dNLL = %.2f (3 parameters)\n", b[["int"]], b[["noint"]],
            2 * (b[["noint"]] - b[["int"]])))
saveRDS(list(fits = R, best = b), "check_hs_intercept_r2_20261102.rds")
cat(sprintf("DONE %s\n", format(Sys.time())))
