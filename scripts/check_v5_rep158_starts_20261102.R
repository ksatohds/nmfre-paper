# The one v5 replicate flagged as degenerate (gaussian, eps = 0.2, rep 158): record what its saved result
# says about the starts and the likelihood, and compare with the truth-started fit of C75, so that the
# manuscript sentence on it has a log to cite (work memo C88; A91).  NEW file, 2026-11-02.
# Reads out_hs_v5/gaussian_eps_0.20/rep_0158.rds (2.9 KB) and the matching file of out_v4_truthstart; no fitting.
# Usage: Rscript check_v5_rep158_starts_20261102.R   Output: check_v5_rep158_starts_20261102.{log,rds}
f <- "out_hs_v5/gaussian_eps_0.20/rep_0158.rds"
x <- readRDS(f)
ft <- list.files("out_v4_truthstart", pattern = "^rep_.*rds$", recursive = TRUE, full.names = TRUE)
ts <- NULL
for (g in ft) { y <- readRDS(g); if (y$err == "gaussian" && y$eps == 0.2 && y$r == 158) { ts <- y; break } }
sink("check_v5_rep158_starts_20261102.log", split = TRUE)
cat("file", f, "|", format(Sys.time()), "\n")
cat(sprintf("design %s eps %.2f rep %d | degenerate %s | starts: default (k-means) degenerate %s, spa degenerate %s, extra k-means starts %d\n",
            x$err, x$eps, x$r, x$degenerate, x$start$deg_default, x$start$deg_spa, x$start$n_extra))
nl <- x$start$nll_all
cat(sprintf("starts recorded %d | marginal NLL of the starts: min %.4f max %.4f range %.2e | kept fit NLL %.4f\n",
            length(nl), min(nl), max(nl), max(nl) - min(nl), x$nll))
print(round(nl, 4))
if (!is.null(ts)) {
  cat(sprintf("truth-started fit (C75): NLL %.4f | 2*(NLL_v5 - NLL_truthstart) = %.2e | NLL at the true parameters %.4f\n",
              ts$nll_ts, 2 * (x$nll - ts$nll_ts), ts$nll_true))
} else cat("truth-started fit: not found\n")
cat(sprintf("smallest congruence with the population basis %.3f | estimated region width %.3f\n",
            min(x$congruence), x$width))
cat("DONE\n")
sink()
saveRDS(list(file = f, start = x$start, nll = x$nll, degenerate = x$degenerate, congruence = x$congruence,
             width = x$width, truthstart = if (is.null(ts)) NULL else ts[c("nll_ts", "nll_true")]),
        "check_v5_rep158_starts_20261102.rds")
