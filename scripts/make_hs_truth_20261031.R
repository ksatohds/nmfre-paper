# ----------------------------------------------------------------------------
# Truth for the HS1939-based simulation (NEW file, 2026-10-31; work memo A13).
# HS1939 with an intercept row in A, fitted by nmfkc 1.0.2.9000 nmfre() (seeds 1-5 and
# X.init = "spa" were run in check_r2_helpers_20261031.R; the best by marginal NLL is
# refitted here so that the truth comes from the current package version).
# Saves X*, Theta*, sigma2*, tau2*, the design A (P = 9, N = 300) and diagnostics of the
# truth (anchor rows, identification range of the truth = a single point?).
# ----------------------------------------------------------------------------
suppressMessages(pkgload::load_all("D:/My Documents/GitHub/nmfkc", quiet = TRUE,
                                   export_all = TRUE, helpers = FALSE))
setwd("D:/My Documents/Data/NMF-RE混合線形モデルr2/scripts_results")
source("nmfre_r2_helpers_20261031.R")
cat("nmfkc", as.character(packageVersion("nmfkc")), "\n")
data(HolzingerSwineford1939, package = "lavaan")
d <- HolzingerSwineford1939; d <- d[complete.cases(d), ]
d$age.rev <- -(d$ageyr + d$agemo / 12); d$sex.2 <- ifelse(d$sex == 2, 1, 0)
d$school.GW <- ifelse(d$school == "Grant-White", 1, 0)
d[, c("id", "sex", "ageyr", "agemo", "school", "grade")] <- NULL; d <- nmfkc.normalize(d)
ex <- c("age.rev", "sex.2", "school.GW"); en <- setdiff(colnames(d), ex)
Y <- t(d[, en]); A <- rbind(intercept = 1, t(d[, ex]))
fits <- lapply(1:5, function(s) suppressMessages(suppressWarnings(
  nmfre(Y, A, rank = 3, epsilon = 1e-7, maxit = 500000, seed = s, print.trace = FALSE))))
nl <- sapply(fits, function(f) .nmfre.marginal.nll(Y, f$X, f$C, A, f$sigma2, f$tau2))
cat("NLL by seed:", sprintf("%.4f", nl), "\n")
f <- fits[[which.min(nl)]]
o <- r2_orient(f, "geomin"); ft <- o$fit
ir <- r2_ident_range(ft)
cat(sprintf("converged %s | sigma2 %.5f tau2 %.5f lambda %.5f\n", isTRUE(f$converged), f$sigma2, f$tau2, f$lambda))
cat("identification range of the truth: arc", round(ir$arc * 180 / pi, 3), "deg; max width of Theta",
    signif(max(ir$range["max", ] - ir$range["min", ]), 3), "\n")
cat("orientation moved Theta by", signif(max(abs(ft$C - f$C)), 3), "\n")
S <- ft$X / rowSums(ft$X); rownames(S) <- rownames(Y)
cat("row shares of X*:\n"); print(round(S, 3))
cat("Theta*:\n"); print(round(ft$C, 4))
cat("degeneracy of the truth:\n"); str(.nmfre.degeneracy(ft$X, ft$C))
truth <- list(X = ft$X, C = ft$C, sigma2 = f$sigma2, tau2 = f$tau2, A = A,
              Y_obs = Y, nll = min(nl), ident_range = ir, version = as.character(packageVersion("nmfkc")))
saveRDS(truth, "sim_r2_hs_truth_20261031.rds")
cat("\nDONE\n")
