# Summary of the HS simulation v5 (NEW file, 2026-11-02; work memo C79; a copy of summarize_sim_r2_hs_v4_20261101.R
# with the input and output names changed).  Reads out_hs_v5/<err>_eps_<e>/rep_XXXX.rds
# (collected from the three servers) -- no recomputation; light, runs locally.
# Tables by error design x e: fits (degenerate, extra starts, estimated region width vs population), recovery
# (congruence), the two intervals (one-step wild bootstrap, profile likelihood): coverage of the canonical Theta and of
# the whole population region, size at the two true zeros, power for the seven non-zero effects, rejection of 0 where
# 0 lies in the population region, lengths; Xi (NMF-RE vs ESEM: RMSE, coverage, Heywood); covariate-level Wald tests;
# the known-anchor arm (e = 0, Gaussian); time.
# Usage: Rscript summarize_sim_r2_hs_v5_20261102.R      Output: summarize_sim_r2_hs_v5_20261102.{log,rds}
root <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else "out_hs_v5"
dirs <- list.dirs(root, recursive = FALSE)
res <- unlist(lapply(dirs, function(d) lapply(list.files(d, pattern = "^rep_.*rds$", full.names = TRUE), readRDS)),
              recursive = FALSE)
ok <- Filter(function(x) isTRUE(x$ok), res)
cat(sprintf("replicates read %d, ok %d | %s\n", length(res), length(ok), format(Sys.time())))
bad <- Filter(function(x) !isTRUE(x$ok), res)
if (length(bad)) for (x in bad) cat("  failed:", x$err, x$eps, x$r, x$msg, "\n")
D <- do.call(rbind, lapply(ok, function(x) cbind(err = x$err, eps = x$eps, rep = x$r, x$table)))
I <- do.call(rbind, lapply(ok, function(x) data.frame(err = x$err, eps = x$eps, rep = x$r, degenerate = x$degenerate,
  extra = x$start$n_extra, width = x$width, cong_min = min(x$congruence), nll = x$nll, sigma2 = x$sigma2, tau2 = x$tau2,
  xi_cov = mean(x$cov_xi), xi_rmse = sqrt(mean(x$err_xi^2)), esem_ok = isTRUE(x$esem$ok), heywood = isTRUE(x$esem$heywood),
  xi_cov_esem = mean(x$cov_xi_esem), xi_rmse_esem = sqrt(mean(x$err_xi_esem^2)),
  t_total = x$t_total, t_fit = x$t_fit, t_inf = x$t_inf)))
W <- do.call(rbind, lapply(ok, function(x) data.frame(err = x$err, eps = x$eps, rep = x$r, t(x$wald_p))))
options(width = 220)
grp <- function(df) split(df, list(df$err, df$eps), drop = TRUE)

cat("\n== fits ==\n")
print(do.call(rbind, lapply(grp(I), function(d) data.frame(err = d$err[1], eps = d$eps[1], n = nrow(d),
  degenerate = sum(d$degenerate), extra_starts = sum(d$extra > 0), point_est_region = mean(d$width <= 1e-8),
  width_med = median(d$width), cong_med = median(d$cong_min), cong_p05 = quantile(d$cong_min, 0.05),
  sigma2 = mean(d$sigma2), tau2 = mean(d$tau2), t_med = median(d$t_total), t_max = max(d$t_total)))),
  digits = 4, row.names = FALSE)

cat("\n== intervals for the covariate effects (9 pooled; zero = the two true zeros; nonzero = the other seven) ==\n")
tab <- do.call(rbind, lapply(grp(D), function(d) {
  z0 <- d$target == 0; nz <- !z0; inreg <- d$popL <= 0 & 0 <= d$popU & nz
  data.frame(err = d$err[1], eps = d$eps[1], n = nrow(d) / 9,
             cov_wild = mean(d$cov_wild), cov_prof = mean(d$cov_prof, na.rm = TRUE),
             set_wild = mean(d$set_wild), set_prof = mean(d$set_prof, na.rm = TRUE),
             size0_wild = mean(d$rej0_wild[z0]), size0_prof = mean(d$rej0_prof[z0], na.rm = TRUE),
             power_wild = mean(d$rej0_wild[nz]), power_prof = mean(d$rej0_prof[nz], na.rm = TRUE),
             rej0_inregion_wild = if (any(inreg)) mean(d$rej0_wild[inreg]) else NA,
             rej0_inregion_prof = if (any(inreg)) mean(d$rej0_prof[inreg], na.rm = TRUE) else NA,
             len_wild = mean(d$wild_hi - d$wild_lo), len_prof = mean(d$prof_hi - d$prof_lo, na.rm = TRUE),
             prof_reached = mean(d$prof_reached, na.rm = TRUE))
}))
print(tab, digits = 3, row.names = FALSE)
cat("\n== by coefficient (coverage wild / profile; rejection of 0), e = 0 and e = 0.2 ==\n")
bc <- do.call(rbind, lapply(split(D[D$eps %in% c(0, 0.2), ], list(D$err[D$eps %in% c(0, 0.2)], D$eps[D$eps %in% c(0, 0.2)],
  D$coef[D$eps %in% c(0, 0.2)]), drop = TRUE), function(d) data.frame(err = d$err[1], eps = d$eps[1], coef = d$coef[1],
  target = d$target[1], popL = d$popL[1], popU = d$popU[1], snr = abs(d$target[1]) / mean(d$se),
  cov_wild = mean(d$cov_wild), cov_prof = mean(d$cov_prof, na.rm = TRUE),
  rej0_wild = mean(d$rej0_wild), rej0_prof = mean(d$rej0_prof, na.rm = TRUE))))
print(bc[order(bc$err, bc$eps, bc$snr), ], digits = 3, row.names = FALSE)

cat("\n== Xi = X Theta (covariate columns): NMF-RE (C1 Wald) vs ESEM ==\n")
print(do.call(rbind, lapply(grp(I), function(d) data.frame(err = d$err[1], eps = d$eps[1],
  rmse_nmfre = mean(d$xi_rmse), cov_nmfre = mean(d$xi_cov), rmse_esem = mean(d$xi_rmse_esem, na.rm = TRUE),
  cov_esem = mean(d$xi_cov_esem, na.rm = TRUE), esem_ok = mean(d$esem_ok), heywood = mean(d$heywood)))),
  digits = 4, row.names = FALSE)

cat("\n== covariate-level Wald tests (all three covariates are non-null): rejection at 0.05 ==\n")
print(do.call(rbind, lapply(grp(W), function(d) data.frame(err = d$err[1], eps = d$eps[1],
  t(colMeans(d[, setdiff(names(d), c("err", "eps", "rep")), drop = FALSE] < 0.05))))), digits = 3, row.names = FALSE)

cat("\n== known-anchor arm (e = 0, Gaussian) ==\n")
an <- Filter(function(x) !is.null(x$anchor), ok)
if (length(an)) {
  AD <- do.call(rbind, lapply(an, function(x) x$anchor$table)); cost <- sapply(an, function(x) x$anchor$cost)
  cat(sprintf("n %d | cost (marginal NLL, anchored - free): median %.3f, 90th %.3f, max %.3f | degenerate %d | width 0 in %d | cong min median %.4f\n",
              length(an), median(cost), quantile(cost, 0.9), max(cost), sum(sapply(an, function(x) x$anchor$degenerate)),
              sum(sapply(an, function(x) isTRUE(x$anchor$width <= 1e-8))), median(sapply(an, function(x) min(x$anchor$congruence)))))
  cat(sprintf("wild coverage %.3f | size at the zeros %.3f | power (non-zero) %.3f | time median %.0f s\n",
              mean(AD$cov_wild), mean(AD$rej0_wild[AD$target == 0]), mean(AD$rej0_wild[AD$target != 0]),
              median(sapply(an, function(x) x$anchor$time))))
}
saveRDS(list(D = D, I = I, W = W, tab = tab, by_coef = bc), "summarize_sim_r2_hs_v5_20261102.rds")
cat(sprintf("\nDONE %s\n", format(Sys.time())))
