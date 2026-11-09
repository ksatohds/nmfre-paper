# HS simulation v5 (C79): further summaries for the rewrite of Section 5 (A60) that summarize_sim_r2_hs_v5_20261102 did
# not print: (1) the width of the estimated identification region when it is not a point (median, 90th percentile,
# against the mean standard error); (2) calibration of the conditional standard errors: Monte Carlo SD and MAD of the
# estimates over the mean SE, per coefficient; (3) bias of the variance estimates.  Reads
# summarize_sim_r2_hs_v5_20261102.rds only (no recomputation).  (NEW file, 2026-11-02; work memo C88.)
# Usage: Rscript summarize_sim_r2_hs_v5_extra_20261102.R    Output: summarize_sim_r2_hs_v5_extra_20261102.{log,rds}
S <- readRDS("summarize_sim_r2_hs_v5_20261102.rds"); D <- S$D; I <- S$I
s2_true <- 0.01527; t2_true <- 0.26263            # C20 (sim_r2_hs_truth_20261031.rds)
options(width = 200)
cat(sprintf("v5: %d replicates, %d coefficient rows | %s\n", nrow(I), nrow(D), format(Sys.time())))
## (1) widths of the estimated identification region
w <- do.call(rbind, lapply(split(I, list(I$err, I$eps), drop = TRUE), function(d) {
  dd <- D[D$err == d$err[1] & D$eps == d$eps[1], ]
  nz <- d$width > 1e-8
  data.frame(err = d$err[1], eps = d$eps[1], n = nrow(d), point = mean(!nz),
             width_med_nonpoint = if (any(nz)) median(d$width[nz]) else NA,
             width_p90_nonpoint = if (any(nz)) quantile(d$width[nz], 0.9) else NA,
             coefwidth_med_nonpoint = { cw <- dd$U - dd$L; cw <- cw[cw > 1e-8]; if (length(cw)) median(cw) else NA },
             mean_se = mean(dd$se), sigma2_mean = mean(d$sigma2), tau2_mean = mean(d$tau2),
             sigma2_relbias = mean(d$sigma2) / s2_true - 1, tau2_relbias = mean(d$tau2) / t2_true - 1)
}))
w <- w[order(w$err, w$eps), ]
cat("\n== (1) estimated identification region and (3) variances ==\n"); print(w, row.names = FALSE, digits = 4)
## (2) SE calibration per coefficient (non-degenerate replicates)
ok <- merge(D, I[, c("err", "eps", "rep", "degenerate")], by = c("err", "eps", "rep"))
ok <- ok[!ok$degenerate, ]
cal <- do.call(rbind, lapply(split(ok, list(ok$err, ok$eps, ok$coef), drop = TRUE), function(d)
  data.frame(err = d$err[1], eps = d$eps[1], coef = d$coef[1], n = nrow(d), mean_se = mean(d$se),
             sd_est = sd(d$est), mad_est = mad(d$est), sd_over_se = sd(d$est) / mean(d$se),
             mad_over_se = mad(d$est) / mean(d$se))))
cal <- cal[order(cal$err, cal$eps, cal$coef), ]
rng <- do.call(rbind, lapply(split(cal, list(cal$err, cal$eps), drop = TRUE), function(d)
  data.frame(err = d$err[1], eps = d$eps[1], sd_over_se_min = min(d$sd_over_se), sd_over_se_max = max(d$sd_over_se),
             mad_over_se_min = min(d$mad_over_se), mad_over_se_max = max(d$mad_over_se))))
rng <- rng[order(rng$err, rng$eps), ]
cat("\n== (2) conditional SE against the spread of the estimates (9 coefficients; range over coefficients) ==\n")
print(rng, row.names = FALSE, digits = 3)
cat("\nper coefficient, gaussian e = 0:\n"); print(cal[cal$err == "gaussian" & cal$eps == 0, ], row.names = FALSE, digits = 3)
saveRDS(list(width = w, se_calibration = cal, se_range = rng), "summarize_sim_r2_hs_v5_extra_20261102.rds")
cat(sprintf("DONE %s\n", format(Sys.time())))
