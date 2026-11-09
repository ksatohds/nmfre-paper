# Width of the population identification region of each covariate effect in the v5 designs (C72: largest
# width 0, 0.084, 0.187, 0.475 for eps = 0, 0.05, 0.1, 0.2), recorded from the saved per-replicate tables so
# that the manuscript can cite a log (A91).  NEW file, 2026-11-02.  popL/popU are the bounds of the
# population region (the same in every replicate of a design); reads the first replicate of each design.
# Usage: Rscript check_v5_population_width_20261102.R   Output: check_v5_population_width_20261102.{log,rds}
des <- expand.grid(err = c("gaussian", "varhet"), eps = c("0.00", "0.05", "0.10", "0.20"), stringsAsFactors = FALSE)
out <- do.call(rbind, lapply(seq_len(nrow(des)), function(i) {
  d <- sprintf("out_hs_v5/%s_eps_%s", des$err[i], des$eps[i])
  x <- readRDS(list.files(d, pattern = "^rep_.*rds$", full.names = TRUE)[1])
  tb <- x$table; w <- tb$popU - tb$popL
  data.frame(err = des$err[i], eps = as.numeric(des$eps[i]), max_width = max(w), coef_max = tb$coef[which.max(w)],
             n_positive = sum(w > 1e-8))
}))
sink("check_v5_population_width_20261102.log", split = TRUE)
cat("v5 population identification region of the nine covariate effects |", format(Sys.time()), "\n")
print(out, row.names = FALSE, digits = 4)
cat("DONE\n")
sink()
saveRDS(out, "check_v5_population_width_20261102.rds")
