# Orthodont table of Appendix F without the intercept row (A99): the intercept is reported as an estimate only, in
# the caption, as the component intercepts of the main text.  NEW file, 2026-11-02 (emit_tab_orthodont_r2_20261102.R
# and its tab_orthodont_theta.tex are not modified).  Reads refit_examples_r2_20261102.rds (C81) only.
# Usage: Rscript emit_tab_orthodont_v2_r2_20261102.R
# Output: emit_tab_orthodont_v2_r2_20261102.{log,rds}, tab_orthodont_theta_v2.tex
z <- readRDS("refit_examples_r2_20261102.rds")$orthodont
cf <- z$inf$coefficients
pv <- function(p) ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p))
ic <- cf[cf$Covariate == "intercept", ]; mc <- cf[cf$Covariate == "male", ]
tx <- sprintf("Male & %.3f & %.3f & %.3f & %.2f & %s \\\\", mc$Estimate, mc$SE, mc$BSE, mc$z_value, pv(mc$p_value))
writeLines(tx, "tab_orthodont_theta_v2.tex")
sink("emit_tab_orthodont_v2_r2_20261102.log", split = TRUE)
cat("emit_tab_orthodont_v2_r2_20261102 |", format(Sys.time()), "\n-- tab_orthodont_theta_v2.tex\n"); cat(tx, sep = "\n")
cat(sprintf("Orthodont intercept (estimate only): %.3f\n", ic$Estimate))
cat("DONE\n")
sink()
saveRDS(list(row = tx, intercept = ic$Estimate), "emit_tab_orthodont_v2_r2_20261102.rds")
