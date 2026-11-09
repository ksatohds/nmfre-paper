# Table body and basis values for the Orthodont example of Appendix F (A94).  NEW file, 2026-11-02.
# Reads refit_examples_r2_20261102.rds (C81) only.  The intercept is reported as an estimate only (no z, p), as the
# component intercepts in the main text (A24): its conditional interval ignores the estimation error of X.
# Usage: Rscript emit_tab_orthodont_r2_20261102.R
# Output: emit_tab_orthodont_r2_20261102.{log,rds}, tab_orthodont_theta.tex
z <- readRDS("refit_examples_r2_20261102.rds")$orthodont
cf <- z$inf$coefficients
pv <- function(p) ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p))
ic <- cf[cf$Covariate == "intercept", ]; mc <- cf[cf$Covariate == "male", ]
tx <- c(sprintf("Intercept & %.3f & %.3f & %.3f & -- & -- \\\\", ic$Estimate, ic$SE, ic$BSE),
        sprintf("Male      & %6.3f & %.3f & %.3f & %5.2f & %s \\\\", mc$Estimate, mc$SE, mc$BSE, mc$z_value, pv(mc$p_value)))
writeLines(tx, "tab_orthodont_theta.tex")
X <- as.vector(z$fit$X)
sink("emit_tab_orthodont_r2_20261102.log", split = TRUE)
cat("emit_tab_orthodont_r2_20261102 |", format(Sys.time()), "\n-- tab_orthodont_theta.tex\n"); cat(tx, sep = "\n")
cat(sprintf("Orthodont basis X: %s | sum %.6f\n", paste(sprintf("%.4f", X), collapse = ", "), sum(X)))
cat(sprintf("Orthodont dfU/(NQ) %.3f\n", z$fit$dfU.frac))
cat("DONE\n")
sink()
saveRDS(list(rows = tx, X = X), "emit_tab_orthodont_r2_20261102.rds")
