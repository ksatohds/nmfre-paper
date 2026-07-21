# 2026-07-19: copy of show_HS_X_cluster.R repointed at the CURRENT core.
# The original still sources nmfre_core20260627_fullEM.R, so the basis table it produced
# (main-text Table 4, tab:hs1939_basis) was left at the pre-2026-07-18 core while the
# surrounding tables were updated.  Existing scripts are not modified; this is a new file.
# Show the optimized basis X (9 tests x 3 factors) for the HS fit, plus the
# row-normalized X.prob and the hard cluster label X.cluster (= which factor
# each observed test loads on). nmf.re does NOT return X.cluster; this mirrors
# nmfkc()'s definition: X.prob = X/rowSums(X); X.cluster = apply(X.prob,1,which.max).
suppressPackageStartupMessages({library(lavaan); library(nmfkc)})
source("nmfre_core20260718.R")

data(HolzingerSwineford1939); d <- HolzingerSwineford1939; d <- d[complete.cases(d), ]
d$age.rev <- -(d$ageyr + d$agemo/12); d$sex.2 <- ifelse(d$sex==2,1,0)
d$school.GW <- ifelse(d$school=="Grant-White",1,0)
d[, c("id","sex","ageyr","agemo","school","grade")] <- NULL
d <- nmfkc::nmfkc.normalize(d)
exo <- c("age.rev","sex.2","school.GW"); endo <- setdiff(colnames(d), exo)
Y <- t(d[, endo]); A <- t(d[, exo]); Q0 <- 3
ini <- nmfkc::nmfkc(Y, A, Q = Q0, epsilon = 1e-8, method = "EU", seed = 1)
fit <- nmf.re(Y, A, Q = Q0, X.ini = ini$X, Theta.ini = ini$C, prefix = "Factor",
              theta_sign = "free", x_update = "seminmf", dfU_control = "off",
              sigma2 = 1, tau2 = 1, sigma2_update = TRUE, tau2_update = TRUE,
              inference = FALSE, verbose = FALSE, post_message = FALSE)
X <- fit$X; rownames(X) <- endo; colnames(X) <- paste0("Factor", 1:Q0)

cat("names(fit) contains 'X.cluster'? ->", "X.cluster" %in% names(fit), "\n\n")

cat("=== Optimized basis X (rows = tests, cols = factors; columns sum to 1) ===\n")
print(round(X, 3))
cat("colSums(X):", round(colSums(X), 3), "\n\n")

Xp <- X / (rowSums(X) + 1e-12)                 # row-normalized membership
cat("=== X.prob = X / rowSums(X)  (each row sums to 1: membership share) ===\n")
print(round(Xp, 3))

Xcl <- apply(Xp, 1, which.max)
cat("\n=== X.cluster = argmax over factors for each test ===\n")
print(data.frame(test = endo, cluster = paste0("Factor", Xcl),
                 dominant_share = round(apply(Xp, 1, max), 3)), row.names = FALSE)

cat("\n=== Factor -> tests grouping ===\n")
for (q in 1:Q0) cat(sprintf("  Factor%d: %s\n", q, paste(endo[Xcl == q], collapse = ", ")))
