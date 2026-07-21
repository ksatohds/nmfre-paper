# ----------------------------------------------------------------------------
# Aggregate the simulation output into the three LaTeX table fragments of the paper:
#   tab_recovery.tex  (Appendix Table D1: bias/RMSE, all four methods)
#   tab_coverage.tex  (Appendix Table D2: C1 vs C2 coverage, NMF-RE arms)
#   tab_signed.tex    (Table 1: negative-effect head-to-head, free vs nonneg)
#
# Input: the two raw per-replicate CSVs written by
#   NMF-RE_SimRecoveryCovECM_20260718.R
# (one per error distribution), i.e.
#   sim_recovery_CovECM_20260718_raw_orth_gaussian.csv
#   sim_recovery_CovECM_20260718_raw_orth_exp_centered.csv
#
# RMSE of the composite covariate effect (RMSE(Xi)) is a POOLED Monte Carlo RMSE,
#   sqrt( mean over replicates and P elements of squared error ),
# which equals sqrt( mean_i Cmale_rmse_i^2 ) because Cmale_rmse_i is the per-replicate
# RMSE over the P elements.  Bias(theta), RMSE(theta) and the coverages are plain means.
#
# Run from the repository root after the simulation has produced the two raw CSVs:
#   Rscript make_tables.R
# ----------------------------------------------------------------------------

gfile <- "sim_recovery_CovECM_20260718_raw_orth_gaussian.csv"
efile <- "sim_recovery_CovECM_20260718_raw_orth_exp_centered.csv"
if (!file.exists(gfile) || !file.exists(efile))
  stop("Run NMF-RE_SimRecoveryCovECM_20260718.R first to produce the raw CSVs.")
raw <- rbind(read.csv(gfile, stringsAsFactors = FALSE),
             read.csv(efile, stringsAsFactors = FALSE))

agg <- do.call(rbind, lapply(
  split(raw, list(raw$method, raw$N, raw$scenario, raw$err_dist), drop = TRUE),
  function(d) {
    th <- d$theta_true[1]
    data.frame(
      method = d$method[1], N = d$N[1], scenario = d$scenario[1], err_dist = d$err_dist[1],
      Bias_theta  = if (all(is.na(d$est))) NA else mean(d$est - th, na.rm = TRUE),
      RMSE_theta  = if (all(is.na(d$est))) NA else sqrt(mean((d$est - th)^2, na.rm = TRUE)),
      RMSE_Cmale  = sqrt(mean(d$Cmale_rmse^2, na.rm = TRUE)),          # pooled MC RMSE
      Bias_Cmale  = mean(d$Cmale_bias, na.rm = TRUE),
      Cover_C1_bse = mean(d$cover_C1_bse, na.rm = TRUE),
      Cover_C2_se  = mean(d$cover_C2_se,  na.rm = TRUE),
      dfU = mean(d$dfU_frac, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
rownames(agg) <- NULL

fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
meth_lab <- c(NMFRE_free = "NMF-RE (free)", NMFRE_nonneg = "NMF-RE (nonneg)",
              MIMIC = "MIMIC", OLS = "OLS")
scen_lab <- c(A_null = "Null", B_alt = "Positive", C_neg = "Negative")
err_lab  <- c(gaussian = "Gaussian", exp_centered = "Exp.")
agg$scen_l <- scen_lab[agg$scenario]; agg$err_l <- err_lab[agg$err_dist]
agg$err_ord <- match(agg$err_dist, c("gaussian", "exp_centered"))

blk <- function(d, cols, header, file, comment) {
  d <- d[order(d$scen_l, d$N, d$method, d$err_ord), ]
  con <- file(file, "w", encoding = "UTF-8")
  writeLines(c(comment, "\\begin{tabular}{lrllrrr}", "\\toprule", header, "\\midrule"), con)
  prev <- ""
  for (i in seq_len(nrow(d))) {
    key <- paste(d$scen_l[i], d$N[i]); nb <- key != prev
    if (i > 1 && nb) writeLines("\\midrule", con)
    sc <- if (nb) d$scen_l[i] else ""; nn <- if (nb) as.character(d$N[i]) else ""
    prev <- key; r <- d[i, ]
    writeLines(sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", sc, nn,
      meth_lab[r$method], r$err_l, fmt(cols[[1]](r)), fmt(cols[[2]](r)), fmt(cols[[3]](r))), con)
  }
  writeLines(c("\\bottomrule", "\\end{tabular}"), con); close(con)
}

## Appendix Table D1: recovery (all four methods) --- RMSE(Xi) pooled
blk(agg, list(function(r) r$Bias_theta, function(r) r$RMSE_theta, function(r) r$RMSE_Cmale),
    "Scenario & $N$ & Method & Error & Bias($\\theta$) & RMSE($\\theta$) & RMSE($\\Xi$) \\\\",
    "tab_recovery.tex",
    "% recovery; variances by ECM; RMSE(Xi) pooled over replicates and elements")

## Appendix Table D2: coverage (NMF-RE arms; C1 vs C2)
blk(subset(agg, method %in% c("NMFRE_free", "NMFRE_nonneg")),
    list(function(r) r$Cover_C1_bse, function(r) r$Cover_C2_se, function(r) r$dfU),
    "Scenario & $N$ & Method & Error & C1 (cond.) & C2 (refit) & $\\mathrm{df}_U/(NQ)$ \\\\",
    "tab_coverage.tex",
    "% C1 (bootstrap-SD Wald) vs C2 (full-refit) coverage, nominal 0.95, full ECM")

## Table 1: negative scenario, Gaussian, free vs nonneg
sg <- subset(agg, scenario == "C_neg" & err_dist == "gaussian" &
               method %in% c("NMFRE_free", "NMFRE_nonneg"))
sg <- sg[order(sg$N, sg$method), ]
con <- file("tab_signed.tex", "w", encoding = "UTF-8")
writeLines(c("% Table 1: C_neg (true theta_male = -9.43), free vs nonneg; RMSE(Xi) pooled",
             "\\begin{tabular}{rlrrrr}", "\\toprule",
             "$N$ & Method & Bias($\\theta$) & RMSE($\\Xi$) & C1 cover & C2 cover \\\\",
             "\\midrule"), con)
prev <- -1
for (i in seq_len(nrow(sg))) {
  r <- sg[i, ]
  if (i > 1 && r$N != prev) writeLines("\\midrule", con)
  nn <- if (r$N != prev) as.character(r$N) else ""; prev <- r$N
  writeLines(sprintf("%s & %s & %s & %s & %s & %s \\\\", nn, meth_lab[r$method],
    fmt(r$Bias_theta), fmt(r$RMSE_Cmale), fmt(r$Cover_C1_bse), fmt(r$Cover_C2_se)), con)
}
writeLines(c("\\bottomrule", "\\end{tabular}"), con); close(con)

cat("Wrote tab_recovery.tex, tab_coverage.tex, tab_signed.tex\n")
