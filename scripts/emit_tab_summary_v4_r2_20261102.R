# Rows of the summary table tab:example_summary after the topic model was moved back to the main text (A110): the two
# main-text examples (HS, topic model) first, then the two of Appendix F.  NEW file, 2026-11-02
# (emit_tab_summary_v3_r2_20261102.R and its tab_ex_summary_v3.tex are not modified).  Same formatting; only the rule
# moves.  Reads refit_examples_r2_20261102.rds (C81) only.
# Usage: Rscript emit_tab_summary_v4_r2_20261102.R   Output: emit_tab_summary_v4_r2_20261102.{log,rds}, tab_ex_summary_v4.tex
new <- readRDS("refit_examples_r2_20261102.rds")
lab <- c(hs = "HS1939", topic = "Topic model", orthodont = "Orthodont", odhub = "OD--HUB")
fs <- function(x) if (x < 0.1) sprintf("%.4f", x) else if (x < 10) sprintf("%.3f", x) else if (x < 1000) sprintf("%.2f", x) else sprintf("%.1f", x)
tx <- sapply(names(lab), function(nm) {
  z <- new[[nm]]; f <- z$fit
  sprintf("%-11s & %d & %d & %d & %s & %s & %.4f & %.3f \\\\", lab[[nm]], z$N, z$P, z$Q, fs(f$sigma2), fs(f$tau2),
          f$lambda, f$dfU.frac)
})
tx <- append(tx, "\\midrule", after = 2)
writeLines(tx, "tab_ex_summary_v4.tex")
sink("emit_tab_summary_v4_r2_20261102.log", split = TRUE)
cat("emit_tab_summary_v4_r2_20261102 |", format(Sys.time()), "\n-- tab_ex_summary_v4.tex\n"); cat(tx, sep = "\n")
v3 <- readLines("tab_ex_summary_v3.tex")
cat(sprintf("data rows identical to tab_ex_summary_v3.tex: %s\n", identical(v3[v3 != "\\midrule"], unname(tx[tx != "\\midrule"]))))
cat("DONE\n")
sink()
saveRDS(tx, "emit_tab_summary_v4_r2_20261102.rds")
