# LaTeX rows and Figure 4 from the re-fit of the five examples with the nmfkc default of a343c39 (C81):
# appendix tables tab:topic_theta and tab:odhub_theta (components labelled as in R1, matched by the X congruence of
# compare_examples_r2_20261102.R), the rows of the summary table tab:example_summary, and Figure 4 drawn as in
# refit_examples_r2_20261031.R.  (NEW file, 2026-11-02; work memo C82.)  Reads RDS only (topic A from
# topic_data_20261031.rds); no OD data.
# Usage: Rscript emit_tables_examples_r2_20261102.R
# Output: emit_tables_examples_r2_20261102.{log,rds}, Figure4_r2_20261102.pdf (the manuscript's Figure4.pdf is not touched)
new <- readRDS("refit_examples_r2_20261102.rds"); cmp <- readRDS("compare_examples_r2_20261102.rds")
pv <- function(p) ifelse(p < 0.001, "$<0.001$", sprintf("%.4f", p))
## rows of the new coefficient table in the R1 order (vec order, R1 labels): new row index for each R1 entry
rows_r1 <- function(nm, key) {
  z <- new[[nm]]; p <- cmp[[key]]$perm; Q <- ncol(z$fit$X); K <- ncol(z$fit$C)
  z$inf$coefficients[as.vector(matrix(seq_len(Q * K), Q, K)[p, , drop = FALSE]), ]
}
fmt <- function(cf, lab_cov, lab_basis, w) {
  sprintf("%-*s & %s & %7.3f & %6.3f & %6.3f & %5.2f & %s \\\\", w, lab_cov, lab_basis, cf$Estimate, cf$SE, cf$BSE,
          cf$z_value, pv(cf$p_value))
}
cat("%% tab:topic_theta (body)\n")
tc <- rows_r1("topic", "topic C55 (mu, draft 3.7)")
era <- rep(c("Early (1775--1865)", "Industrial (1865--1945)", "Postwar (1945--)"), each = 3)
cat(fmt(tc, era, rep(sprintf("Topic %d", 1:3), 3), 23), sep = "\n")
cat("\n%% tab:odhub_theta (body)\n")
oc <- rows_r1("odhub", "odhub C53 (mu, R1 table)")
hub <- rep(c("Tokyo", "Osaka", "Fukuoka"), each = 3)
cat(fmt(oc, hub, rep(sprintf("Region%d", 1:3), 3), 8), sep = "\n")
cat("\n%% tab:example_summary (rows; R1 column order N, P, Q, NQ, sigma2, tau2, lambda, df_U, r)\n")
lab <- c(hs = "HS1939", college = "College", orthodont = "Orthodont", topic = "Topic model", odhub = "OD--HUB")
for (nm in names(lab)) {
  z <- new[[nm]]; f <- z$fit
  s2 <- if (f$sigma2 < 0.1) sprintf("%.5f", f$sigma2) else sprintf("%.3f", f$sigma2)
  t2 <- if (f$tau2 < 1) sprintf("%.4f", f$tau2) else if (f$tau2 < 1000) sprintf("%.2f", f$tau2) else sprintf("%.1f", f$tau2)
  cat(sprintf("%-11s & %-3d & %-2d & %d & %-4d & %s & %s & %.4f & %.2f & %.3f \\\\\n", lab[[nm]], z$N, z$P, z$Q,
              z$N * z$Q, s2, t2, f$lambda, f$dfU, f$dfU.frac))
}
## Figure 4 (topic proportions of the fixed effects), as refit_examples_r2_20261031.R, components in R1 order
tp <- readRDS("topic_data_20261031.rds"); At <- tp$A
p <- cmp[["topic C55 (mu, draft 3.7)"]]$perm; Cn <- new$topic$fit$C[p, , drop = FALSE]
stopifnot(ncol(Cn) == nrow(At))
Bf <- pmax(Cn %*% At, 0); cs <- colSums(Bf); cs[cs == 0] <- 1; Bp4 <- sweep(Bf, 2, cs, "/")
lab4 <- c("Topic1 (institutional)", "Topic2 (union/foreign)", "Topic3 (civic)"); mycol4 <- c("white", "gray80", "gray90")
pdf("Figure4_r2_20261102.pdf", height = 7, width = 10)
par(mar = c(8,4,1,1) + 0.1, cex = 1.0)
barplot(Bp4, col = mycol4, legend = FALSE, las = 3, border = "gray70", cex.names = 0.8, ylab = "Topic proportion")
legend("topright", legend = lab4, bg = "white", fill = mycol4)
dev.off()
Xp <- new$topic$fit$X[, p] / rowSums(new$topic$fit$X)
rownames(Xp) <- rownames(tp$Y)
for (q in 1:3) cat(sprintf("Topic%d words with share > 0.5: %s\n", q, paste(rownames(Xp)[Xp[, q] > 0.5], collapse = ", ")))
cat("era proportions (Topic1, Topic2, Topic3):\n")
print(round(t(sapply(split(seq_len(ncol(At)), apply(At, 2, which.max)), function(j) rowMeans(Bp4[, j, drop = FALSE]))), 3))
cat("WROTE Figure4_r2_20261102.pdf\n")
## (added 2026-11-02, author: keep script, log and RDS together) the emitted rows and the Figure 4 proportions
saveRDS(list(topic_rows = tc, odhub_rows = oc, topic_perm = p, Bp4 = Bp4, era_props =
               t(sapply(split(seq_len(ncol(At)), apply(At, 2, which.max)), function(j) rowMeans(Bp4[, j, drop = FALSE])))),
        "emit_tables_examples_r2_20261102.rds")
