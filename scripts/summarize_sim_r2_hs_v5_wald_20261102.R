# Wald intervals for the HS simulation v5 (NEW file, 2026-11-02; work memo C79, A81).
# Author: "RDSが残っているなら，waldは再計算できるのでは？".  v5 computed the one-step wild bootstrap and the profile
# likelihood intervals for Theta, but not the Wald interval.  Each replicate's table (out_hs_v5/<err>_eps_<e>/rep_XXXX.rds,
# $table) keeps the estimate (est), its standard error from nmfre.inference (se = C.se) and the target, so the Wald
# interval est +- z * se is computed here without refitting.  The wild and profile columns are recomputed from the same
# tables as a check against summarize_sim_r2_hs_v5_20261102.log.
# Light; runs locally.  Usage: Rscript summarize_sim_r2_hs_v5_wald_20261102.R
# Output: summarize_sim_r2_hs_v5_wald_20261102.{log,rds}
root <- "out_hs_v5"
z <- qnorm(0.975)
res <- unlist(lapply(list.dirs(root, recursive = FALSE), function(d)
  lapply(list.files(d, pattern = "^rep_.*rds$", full.names = TRUE), readRDS)), recursive = FALSE)
ok <- Filter(function(x) isTRUE(x$ok), res)
cat(sprintf("replicates read %d, ok %d | %s\n", length(res), length(ok), format(Sys.time())))
D <- do.call(rbind, lapply(ok, function(x) cbind(err = x$err, eps = x$eps, rep = x$r, degenerate = isTRUE(x$degenerate), x$table)))
D$wald_lo <- D$est - z * D$se
D$wald_hi <- D$est + z * D$se
D$cov_wald <- D$wald_lo <= D$target & D$target <= D$wald_hi
D$set_wald <- D$wald_lo <= D$popL & D$popU <= D$wald_hi
D$rej0_wald <- D$wald_lo > 0 | D$wald_hi < 0
D$len_wald <- D$wald_hi - D$wald_lo
D$len_wild <- D$wild_hi - D$wild_lo
D$len_prof <- D$prof_hi - D$prof_lo

options(width = 200)
by_design <- do.call(rbind, lapply(split(D, list(D$err, D$eps), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], n = nrow(d),
  cov_wald = mean(d$cov_wald), cov_wild = mean(d$cov_wild), cov_prof = mean(d$cov_prof, na.rm = TRUE),
  size0_wald = mean(d$rej0_wald[d$target == 0]), size0_wild = mean(d$rej0_wild[d$target == 0]),
  size0_prof = mean(d$rej0_prof[d$target == 0], na.rm = TRUE),
  power_wald = mean(d$rej0_wald[d$target != 0]), power_wild = mean(d$rej0_wild[d$target != 0]),
  power_prof = mean(d$rej0_prof[d$target != 0], na.rm = TRUE),
  set_wald = mean(d$set_wald), set_wild = mean(d$set_wild), set_prof = mean(d$set_prof, na.rm = TRUE),
  len_wald = mean(d$len_wald), len_wild = mean(d$len_wild), len_prof = mean(d$len_prof, na.rm = TRUE))))
by_design <- by_design[order(by_design$err, by_design$eps), ]
cat("\n== by design: coverage of the canonical Theta, size at the true zeros, power, coverage of the whole population region, mean length ==\n")
print(by_design, row.names = FALSE, digits = 3)

by_coef <- do.call(rbind, lapply(split(D, list(D$err, D$eps, D$coef), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], coef = d$coef[1], target = d$target[1],
  cov_wald = mean(d$cov_wald), cov_wild = mean(d$cov_wild), cov_prof = mean(d$cov_prof, na.rm = TRUE))))
by_coef <- by_coef[order(by_coef$err, by_coef$eps, by_coef$coef), ]
cat("\n== by coefficient, e = 0 and e = 0.2 ==\n")
print(by_coef[by_coef$eps %in% c(0, 0.2), ], row.names = FALSE, digits = 3)
cat("\nrange over coefficients of |cov_wald - cov_wild|:", format(range(abs(by_coef$cov_wald - by_coef$cov_wild)), digits = 3), "\n")

# known-anchor arm (e = 0, Gaussian): same tables, wild bootstrap only
A <- do.call(rbind, lapply(ok, function(x) if (!is.null(x$anchor) && !is.null(x$anchor$table))
  cbind(err = x$err, eps = x$eps, rep = x$r, x$anchor$table)))
if (!is.null(A) && nrow(A)) {
  A$cov_wald <- abs(A$est - A$target) <= z * A$se
  A$rej0_wald <- abs(A$est) > z * A$se
  cat(sprintf("\n== known-anchor arm (n = %d rows): coverage Wald %.3f wild %.3f | size at zeros Wald %.3f wild %.3f | power Wald %.3f wild %.3f ==\n",
              nrow(A), mean(A$cov_wald), mean(A$cov_wild), mean(A$rej0_wald[A$target == 0]), mean(A$rej0_wild[A$target == 0]),
              mean(A$rej0_wald[A$target != 0]), mean(A$rej0_wild[A$target != 0])))
}

saveRDS(list(by_design = by_design, by_coef = by_coef, z = z, n_ok = length(ok),
             note = "Wald = est +- z * se (se = nmfre.inference C.se); wild and profile recomputed from the same tables"),
        "summarize_sim_r2_hs_v5_wald_20261102.rds")
cat("\nsaved summarize_sim_r2_hs_v5_wald_20261102.rds\n")
