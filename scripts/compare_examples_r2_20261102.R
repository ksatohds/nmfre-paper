# Compare the re-fit of the five examples with the nmfkc default of a343c39 (refit_examples_r2_20261102.rds, C81) with
# the earlier fits behind the tables: C46 HS, C40 college, C55 Orthodont and topic (draft section 3.7), C53 OD-HUB,
# and C69 topic/OD-HUB (exact non-negative C-step).  Components are matched by the Tucker congruence of X.
# (NEW file, 2026-11-02; work memo C81.)  Reads RDS only; no data.
# Usage: Rscript compare_examples_r2_20261102.R      Output: compare_examples_r2_20261102.{log,rds}
new <- readRDS("refit_examples_r2_20261102.rds")
hs <- readRDS("refit_hs_r2_20261031.rds"); co <- readRDS("refit_college_r2_20261031.rds")
e31 <- readRDS("refit_examples_r2_20261031.rds"); od <- readRDS("refit_odhub_appendix_20261031.rds")
ex69 <- readRDS("refit_nonneg_examples_exact_20261101.rds")
old <- list(
  hs        = list(list(lab = "C46", fit = hs$stage1$free_int$best, cf = hs$inference$free_int)),
  college   = list(list(lab = "C40", fit = co$stage1$free_int$best, cf = co$inference$free_int)),
  orthodont = list(list(lab = "C55", fit = e31$orthodont$best, cf = e31$orthodont$coefficients)),
  topic     = list(list(lab = "C55 (mu, draft 3.7)", fit = e31$topic$best, cf = e31$topic$coefficients),
                   list(lab = "C69 (exact C)", fit = ex69$topic$best, cf = ex69$topic$coefficients)),
  odhub     = list(list(lab = "C53 (mu, R1 table)", fit = od$best, cf = od$coefficients),
                   list(lab = "C69 (exact C)", fit = ex69$odhub$best, cf = ex69$odhub$coefficients)))
perms <- function(n) if (n == 1) matrix(1L) else do.call(rbind, lapply(seq_len(n), function(i)
  cbind(i, perms(n - 1)[, , drop = FALSE] + (perms(n - 1) >= i))))
congr <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
## p[q] = column of the new fit matched to column q of the old fit
match_cols <- function(Xn, Xo) {
  P <- perms(ncol(Xo)); s <- apply(P, 1, function(p) sum(sapply(seq_along(p), function(q) congr(Xn[, p[q]], Xo[, q]))))
  P[which.max(s), ]
}
options(width = 200); out <- list()
for (nm in names(old)) {
  z <- new[[nm]]; cat(sprintf("\n==== %s ====\n", nm))
  if (!isTRUE(z$ok)) { cat("new fit FAILED\n"); next }
  f <- z$fit; cn <- z$inf$coefficients; Q <- ncol(f$X); K <- ncol(f$C)
  for (o in old[[nm]]) {
    p <- match_cols(f$X, o$fit$X)
    cg <- sapply(seq_len(Q), function(q) congr(f$X[, p[q]], o$fit$X[, q]))
    ## coefficient tables are in vec(C) order (Basis fastest); rows of the new table reordered to the old labels
    idx <- as.vector(matrix(seq_len(Q * K), Q, K)[p, , drop = FALSE])
    a <- cn[idx, ]
    t <- data.frame(Basis_old = o$cf$Basis, Covariate = o$cf$Covariate, est_old = o$cf$Estimate, est_new = a$Estimate,
                    se_old = o$cf$SE, se_new = a$SE, bse_old = o$cf$BSE, bse_new = a$BSE,
                    p_old = o$cf$p_value, p_new = a$p_value)
    if (!is.null(a$prof_lo)) { t$prof_lo <- a$prof_lo; t$prof_hi <- a$prof_hi }
    if (is.null(a$prof_lo) && !is.null(z$inf$C.profile.lower)) {
      t$prof_lo <- as.vector(z$inf$C.profile.lower)[idx]; t$prof_hi <- as.vector(z$inf$C.profile.upper)[idx] }
    sig_change <- sum((o$cf$p_value < 0.05) != (a$p_value < 0.05))
    sz <- function(g) c(sigma2 = g$sigma2, tau2 = g$tau2, lambda = g$sigma2 / g$tau2)
    cat(sprintf("-- vs %s | X congruence %s | max |X diff| %.2e | max |est diff| %.3g (max |diff|/SE %.3f) | max |SE diff| %.3g | significance (5%%) changes %d\n",
                o$lab, paste(sprintf("%.6f", cg), collapse = " "), max(abs(f$X[, p] - o$fit$X)),
                max(abs(t$est_new - t$est_old)), max(abs(t$est_new - t$est_old) / t$se_old),
                max(abs(t$se_new - t$se_old)), sig_change))
    cat("   size old:", paste(sprintf("%s %.6g", names(sz(o$fit)), sz(o$fit)), collapse = ", "),
        "| new:", paste(sprintf("%s %.6g", names(sz(f)), sz(f)), collapse = ", "), "\n")
    print(t, row.names = FALSE, digits = 5)
    out[[paste(nm, o$lab)]] <- list(perm = p, congruence = cg, table = t, sig_change = sig_change)
  }
}
saveRDS(out, "compare_examples_r2_20261102.rds")
cat(sprintf("\nDONE %s\n", format(Sys.time())))
