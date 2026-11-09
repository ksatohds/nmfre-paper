# Table bodies of the topic model and OD-HUB (non-negative variant) with the boundary coefficients marked instead of
# given the conventional one-sided p = 0.5 and z = 0, which are not the outcome of a calibrated test (work memo A109,
# R2改定案.md §5.5).  NEW file, 2026-11-02; emit_tab_examples_r2_20261102.R is unchanged.  Reads its RDS only
# (emit_tab_examples_r2_20261102.rds: the rows of the C81 fit in the R1 component order); no fitting, no OD data.
# The active rows are written exactly as before; a row with estimate 0 gets "0" for the estimate, its SE and BSE as
# before, and "---" / "boundary" for z and p.
# Usage: Rscript emit_tab_examples_boundary_r2_20261102.R
# Output: emit_tab_examples_boundary_r2_20261102.{log,rds}, tab_topic_theta_v2.tex, tab_odhub_theta_v2.tex
old <- readRDS("emit_tab_examples_r2_20261102.rds")
pv <- function(p) ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p))
row7 <- function(cf, lab1, lab2) {
  act <- sprintf("%s & %s & %7.3f & %6.3f & %6.3f & %6.2f & %s \\\\", lab1, lab2, cf$Estimate, cf$SE, cf$BSE,
                 cf$z_value, pv(cf$p_value))
  bnd <- sprintf("%s & %s & %7s & %6.3f & %6.3f & %6s & %s \\\\", lab1, lab2, "0", cf$SE, cf$BSE, "---", "boundary")
  ifelse(cf$Estimate == 0, bnd, act)
}
sink("emit_tab_examples_boundary_r2_20261102.log", split = TRUE)
cat("emit_tab_examples_boundary_r2_20261102 |", format(Sys.time()), "\n")
tc <- old$topic_theta
era <- rep(c("Early (1775--1865)     ", "Industrial (1865--1945)", "Postwar (1945--)       "), each = 3)
tx <- row7(tc, era, rep(sprintf("Topic %d", 1:3), 3))
writeLines(tx, "tab_topic_theta_v2.tex")
cat("\n-- tab_topic_theta_v2.tex\n"); cat(tx, sep = "\n")
cat(sprintf("topic: %d boundary rows, %d active rows\n", sum(tc$Estimate == 0), sum(tc$Estimate != 0)))
oc <- old$odhub_theta
hub <- rep(c("Tokyo  ", "Osaka  ", "Fukuoka"), each = 3)
tx2 <- row7(oc, hub, rep(sprintf("Region%d", 1:3), 3))
writeLines(tx2, "tab_odhub_theta_v2.tex")
cat("\n-- tab_odhub_theta_v2.tex\n"); cat(tx2, sep = "\n")
cat(sprintf("odhub: %d boundary rows, %d active rows\n", sum(oc$Estimate == 0), sum(oc$Estimate != 0)))
## the active rows must be identical to the earlier table bodies
same <- function(a, b) { x <- readLines(a); y <- readLines(b); k <- !grepl("boundary", y); identical(x[k], y[k]) }
cat(sprintf("active rows identical to tab_topic_theta.tex: %s; to tab_odhub_theta.tex: %s\n",
            same("tab_topic_theta.tex", "tab_topic_theta_v2.tex"), same("tab_odhub_theta.tex", "tab_odhub_theta_v2.tex")))
cat("DONE\n")
sink()
saveRDS(list(topic = tx, odhub = tx2), "emit_tab_examples_boundary_r2_20261102.rds")
