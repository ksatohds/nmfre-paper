# Figure B1 redrawn with the iteration count fixed.  NEW file, 2026-11-09 (work memo A211).
# plot_nll_monotone_4ex_r2_20261102.R labels each panel "over z$n iterations", but z$n is the number of evaluations of the
# marginal NLL including the starting value, so every panel label was one too large (4901/973/23/624 vs 4900/972/22/623;
# the text's 6,517 = inner + outer steps is correct).  This script draws the same figure from the same input and labels
# the panels with inner + outer steps (= z$n - 1), the count used in the text.  The earlier script is left unchanged.
# The trace is drawn as steps (type = "s"): joining iterations by straight lines (type = "l") shows a drop between two
# iterations as a diagonal (author, 2026-07-16: "斜めの線があります．階段ではないです").
# Reads check_monotone_r2_20261102.rds (C87) only; no fitting.
# Usage (in scripts_results/): Rscript plot_nll_monotone_4ex_r2_20261109.R
# Output: plot_nll_monotone_4ex_r2_20261109.{log,rds}, Figure_nll_monotone_4ex_r2_20261109.{pdf,png}
res <- readRDS("check_monotone_r2_20261102.rds")
labs <- c(hs = "(a) Holzinger--Swineford (Q=3, sign-free)", topic = "(b) Topic model (Q=3, non-negative)",
          orthodont = "(c) Orthodont (Q=1, sign-free)", odhub = "(d) OD-HUB (Q=3, non-negative)")
cnt <- sapply(names(labs), function(nm) { z <- res[[nm]]
  brk <- z$ends[z$ends < z$n]
  c(evaluations = z$n, inner_steps = z$n - z$outer, outer_steps = z$outer - 1, iterations = z$n - 1,
    ticks = length(brk), inc_inner = z$n_inc_nll_inner, inc_outer = z$n_inc_nll_outer) })
stopifnot(all(cnt["iterations", ] == cnt["inner_steps", ] + cnt["outer_steps", ]))
draw <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 5.2, 2.6, 1.0), mgp = c(2.9, 0.8, 0), cex.lab = 1.05)
  for (nm in names(labs)) {
    z <- res[[nm]]; it <- seq_len(z$n) - 1          # 0 = starting value
    # z$ends = 1-based index of the last evaluation of each inner loop; the outer step that follows gives evaluation
    # ends+1, i.e. iteration ends (0-based), so the tick at ends lies on the drop drawn by type = "s".
    plot(it, z$nll, type = "s", lwd = 2.2, col = "black", xlab = "Iteration",
         ylab = "Marginal negative log-likelihood", main = labs[[nm]])
    brk <- z$ends[z$ends < z$n]
    if (length(brk)) rug(brk, side = 1, ticksize = 0.035, lwd = 0.8, col = "gray35")
    legend("topright", bty = "n", cex = 0.75, text.col = "gray25",
           legend = c(sprintf("%.1f to %.1f over %d iterations", z$nll[1], z$nll[z$n], z$n - 1),
                      sprintf("%d outer variance steps (ticks)", length(brk))))
  }
}
pdf("Figure_nll_monotone_4ex_r2_20261109.pdf", width = 9, height = 7.2); draw(); dev.off()
png("Figure_nll_monotone_4ex_r2_20261109.png", width = 1350, height = 1080, res = 150); draw(); dev.off()
sink("plot_nll_monotone_4ex_r2_20261109.log", split = TRUE)
cat("plot_nll_monotone_4ex_r2_20261109 |", format(Sys.time()), "\n"); print(cnt)
cat(sprintf("four examples: %d iterations in all (inner %d + outer %d), increases of the marginal NLL %d\n",
            sum(cnt["iterations", ]), sum(cnt["inner_steps", ]), sum(cnt["outer_steps", ]),
            sum(cnt["inc_inner", ] + cnt["inc_outer", ])))
cat("panel labels (iterations):", paste(names(labs), cnt["iterations", ], collapse = ", "), "\n")
cat("WROTE Figure_nll_monotone_4ex_r2_20261109.pdf/.png\nDONE\n")
sink()
saveRDS(cnt, "plot_nll_monotone_4ex_r2_20261109.rds")
