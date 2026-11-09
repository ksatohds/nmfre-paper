# Same computation as check_v5_single_start_20261102.R, with an RDS of the result added (author: "RDSも出すスクリプトを作って再集計して";
# keep script, log and RDS together).  The original script is not modified.  (NEW file, 2026-11-02.)
# HS simulation v5 (C79): what would one start (k-means only) have lost against mix2 (k-means and SPA)?  From the
# recorded starts of each replicate (start$nll_default, start$nll_spa, deg_*).  (NEW file, 2026-11-02; work memo C80.)
# Author: "初期値は１個にした方が早いし，ＫＫＴに到達してない場合には乱数を増やせば十分では？"
# Usage: Rscript check_v5_single_start_rds_20261102.R     Output: check_v5_single_start_rds_20261102.log
D <- do.call(rbind, lapply(list.files("out_hs_v5", pattern = "^rep_.*rds$", recursive = TRUE, full.names = TRUE),
  function(f) { x <- readRDS(f); s <- x$start
    data.frame(err = x$err, eps = x$eps, rep = x$r, nll_km = s$nll_default, nll_spa = s$nll_spa,
               deg_km = s$deg_default, deg_spa = s$deg_spa, chosen = s$chosen, n_extra = s$n_extra) }))
D$loss <- 2 * (D$nll_km - pmin(D$nll_km, D$nll_spa))     # what k-means alone would have lost (before restarts)
cat(sprintf("v5 %d replicates | %s\n", nrow(D), format(Sys.time())))
options(width = 200)
print(do.call(rbind, lapply(split(D, list(D$err, D$eps), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], n = nrow(d), km_degenerate = sum(d$deg_km), spa_degenerate = sum(d$deg_spa),
  spa_chosen = sum(d$chosen == "spa"), loss_gt_0.01 = sum(d$loss > 0.01), loss_gt_1 = sum(d$loss > 1),
  loss_gt_1_km_not_degenerate = sum(d$loss > 1 & !d$deg_km), max_loss = max(d$loss)))), row.names = FALSE, digits = 4)
saveRDS(list(D = D), "check_v5_single_start_rds_20261102.rds")   # per replicate: NLLs of the starts, flags, loss
cat(sprintf("all: k-means alone loses > 1 in %d of %d, of which k-means was not flagged degenerate in %d\nDONE %s\n",
            sum(D$loss > 1), nrow(D), sum(D$loss > 1 & !D$deg_km), format(Sys.time())))
