# Facts for the revision text "why X >= 0 for HS1939" (Reviewer 2, R2-3d).  (NEW file, 2026-10-31.)
# (1) Are all correlations among the nine tests positive (positive manifold)?
# (2) Unconstrained 3-factor EFA with oblique geomin rotation (as in ESEM): how many loadings are
#     negative, and how large are they?
# Output: check_hs_positive_manifold_20261031.{log,rds}
setwd("D:/My Documents/Data/NMF-RE混合線形モデルr2/scripts_results")
suppressMessages(library(lavaan))
d <- HolzingerSwineford1939[, paste0("x", 1:9)]
R <- cor(d)
cat("HolzingerSwineford1939 (lavaan", as.character(packageVersion("lavaan")), "), N =", nrow(d), "\n\n")
cat("correlations among x1-x9:\n"); print(round(R, 2))
off <- R[upper.tri(R)]
cat(sprintf("\noff-diagonal: min %.3f, max %.3f, negative %d of %d\n", min(off), max(off), sum(off < 0), length(off)))
fit <- efa(data = d, nfactors = 3, rotation = "geomin", rotation.args = list(geomin.epsilon = 0.01))[[1]]
L <- unclass(lavInspect(fit, "est")$lambda)
cat("\n3-factor EFA, ML, oblique geomin (epsilon 0.01): loadings\n"); print(round(L, 3))
cat(sprintf("negative loadings: %d of %d; most negative %.3f; below -0.10: %d\n",
            sum(L < 0), length(L), min(L), sum(L < -0.10)))
cat("factor correlations:\n"); print(round(lavInspect(fit, "est")$psi, 3))
saveRDS(list(R = R, L = L, psi = lavInspect(fit, "est")$psi), "check_hs_positive_manifold_20261031.rds")
