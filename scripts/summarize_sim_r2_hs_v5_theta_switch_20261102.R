# v5 simulation (C79): recovery of Theta (bias and RMSE of the oriented estimate against the canonical target, by
# coefficient) and the operating characteristics of the switching rule of Section 2.5 (wild bootstrap when the
# estimated identification region is a point, profile likelihood otherwise).  Reviewer 2 asked for the recovery of
# X and Theta; the summaries so far reported coverage only (work memo A98, A100).  NEW file, 2026-11-02.
# Reads the per-replicate files out_hs_v5/<design>/rep_*.rds (no recomputation).  A region is a point when its
# width is <= 1e-8, as in summarize_sim_r2_hs_v5_20261102.R.
# Usage: Rscript summarize_sim_r2_hs_v5_theta_switch_20261102.R
# Output: summarize_sim_r2_hs_v5_theta_switch_20261102.{log,rds}, tab_app_sim_theta.tex
fl <- list.files("out_hs_v5", pattern = "^rep_.*rds$", recursive = TRUE, full.names = TRUE)
rows <- do.call(rbind, lapply(fl, function(f) {
  x <- readRDS(f); if (!isTRUE(x$ok)) return(NULL); tb <- x$table
  data.frame(err = x$err, eps = x$eps, rep = x$r, point = x$width <= 1e-8, coef = tb$coef, target = tb$target,
             est = tb$est, cov_wild = tb$cov_wild, cov_prof = tb$cov_prof, rej_wild = tb$rej0_wild, rej_prof = tb$rej0_prof)
}))
rows$cov_switch <- ifelse(rows$point, rows$cov_wild, rows$cov_prof)
rows$rej_switch <- ifelse(rows$point, rows$rej_wild, rows$rej_prof)
rows$zero <- abs(rows$target) < 1e-12
sink("summarize_sim_r2_hs_v5_theta_switch_20261102.log", split = TRUE)
cat(sprintf("replicates %d, coefficient rows %d | %s\n\n", length(unique(paste(rows$err, rows$eps, rows$rep))), nrow(rows), format(Sys.time())))
cat("== switching rule (wild if the estimated region is a point, else profile): pooled over the nine effects ==\n")
sw <- do.call(rbind, lapply(split(rows, list(rows$err, rows$eps), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], point = mean(d$point), cov_wild = mean(d$cov_wild), cov_prof = mean(d$cov_prof),
  cov_switch = mean(d$cov_switch), size0_switch = mean(d$rej_switch[d$zero]), power_switch = mean(d$rej_switch[!d$zero]))))
sw <- sw[order(match(sw$err, c("gaussian", "varhet")), sw$eps), ]
print(sw, row.names = FALSE, digits = 3)
cat("\n== recovery of Theta: bias and RMSE of the estimate against the canonical target, by coefficient ==\n")
rc <- do.call(rbind, lapply(split(rows, list(rows$err, rows$eps, rows$coef), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], coef = d$coef[1], target = d$target[1], bias = mean(d$est - d$target),
  rmse = sqrt(mean((d$est - d$target)^2)))))
rc <- rc[order(match(rc$err, c("gaussian", "varhet")), rc$eps, rc$coef), ]
print(rc, row.names = FALSE, digits = 3)
cat("\n== pooled over the nine effects ==\n")
pl <- do.call(rbind, lapply(split(rc, list(rc$err, rc$eps), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], max_abs_bias = max(abs(d$bias)), rmse_min = min(d$rmse), rmse_max = max(d$rmse))))
pl <- pl[order(match(pl$err, c("gaussian", "varhet")), pl$eps), ]
print(pl, row.names = FALSE, digits = 3)
## table for Appendix E: equal error variances, eps = 0 and 0.2
pretty <- function(cf) { p <- strsplit(cf, ":B"); sapply(p, function(z) sprintf("\\texttt{%s} $\\to$ %s", z[1], z[2])) }
a <- rc[rc$err == "gaussian" & rc$eps == 0, ]; b <- rc[rc$err == "gaussian" & rc$eps == 0.2, ]
a <- a[order(a$coef), ]; b <- b[order(b$coef), ]; stopifnot(all(a$coef == b$coef))
f3 <- function(x) formatC(x, format = "f", digits = 3)
body <- sprintf("%s & %s & %s & %s & %s & %s \\\\", pretty(a$coef), f3(a$target), f3(a$bias), f3(a$rmse), f3(b$bias), f3(b$rmse))
writeLines(c("\\begin{tabular}{lrrrrr}", "\\toprule",
  "& & \\multicolumn{2}{c}{$\\varepsilon=0$} & \\multicolumn{2}{c}{$\\varepsilon=0.2$} \\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  "Effect & True value & Bias & RMSE & Bias & RMSE \\\\", "\\midrule", body, "\\bottomrule", "\\end{tabular}"),
  "tab_app_sim_theta.tex")
cat("\n-- tab_app_sim_theta.tex\n"); cat(readLines("tab_app_sim_theta.tex"), sep = "\n")
cat("DONE\n")
sink()
saveRDS(list(switch = sw, recovery = rc, pooled = pl), "summarize_sim_r2_hs_v5_theta_switch_20261102.rds")
