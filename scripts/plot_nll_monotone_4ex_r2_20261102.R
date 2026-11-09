# Figure B1 without the college data (removed from the paper, A102): the marginal NLL over the whole run for the four
# examples of the paper, drawn as in check_monotone_r2_20261102.R (monochrome; ticks = outer variance steps), and the
# iteration counts quoted in the text.  NEW file, 2026-11-02 (A107).  Reads check_monotone_r2_20261102.rds (C87) only.
# Usage: Rscript plot_nll_monotone_4ex_r2_20261102.R
# Output: plot_nll_monotone_4ex_r2_20261102.{log,rds}, Figure_nll_monotone_4ex_r2_20261102.{pdf,png}
res <- readRDS("check_monotone_r2_20261102.rds")
labs <- c(hs = "(a) Holzinger--Swineford (Q=3, sign-free)", topic = "(b) Topic model (Q=3, non-negative)",
          orthodont = "(c) Orthodont (Q=1, sign-free)", odhub = "(d) OD-HUB (Q=3, non-negative)")
draw <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 5.2, 2.6, 1.0), mgp = c(2.9, 0.8, 0), cex.lab = 1.05)
  for (nm in names(labs)) {
    z <- res[[nm]]; it <- seq_len(z$n)
    plot(it, z$nll, type = "l", lwd = 2.2, col = "black", xlab = "Update (iteration)",
         ylab = "Marginal negative log-likelihood", main = labs[[nm]])
    brk <- z$ends[z$ends < z$n]
    if (length(brk)) rug(brk, side = 1, ticksize = 0.035, lwd = 0.8, col = "gray35")
    legend("topright", bty = "n", cex = 0.75, text.col = "gray25",
           legend = c(sprintf("%.1f to %.1f over %d iterations", z$nll[1], z$nll[z$n], z$n),
                      sprintf("%d outer variance steps (ticks)", length(brk))))
  }
}
pdf("Figure_nll_monotone_4ex_r2_20261102.pdf", width = 9, height = 7.2); draw(); dev.off()
png("Figure_nll_monotone_4ex_r2_20261102.png", width = 1350, height = 1080, res = 150); draw(); dev.off()
cnt <- sapply(names(labs), function(nm) { z <- res[[nm]]
  c(inner_steps = z$n - z$outer, outer_steps = z$outer - 1, inc_inner = z$n_inc_nll_inner, inc_outer = z$n_inc_nll_outer,
    fixed_inc = z$fixed$n_inc_nll) })
sink("plot_nll_monotone_4ex_r2_20261102.log", split = TRUE)
cat("plot_nll_monotone_4ex_r2_20261102 |", format(Sys.time()), "\n"); print(cnt)
cat(sprintf("four examples: %d iterations in all (inner %d + outer %d), increases of the marginal NLL %d\n",
            sum(cnt["inner_steps", ] + cnt["outer_steps", ]), sum(cnt["inner_steps", ]), sum(cnt["outer_steps", ]),
            sum(cnt["inc_inner", ] + cnt["inc_outer", ])))
cat("WROTE Figure_nll_monotone_4ex_r2_20261102.pdf/.png\nDONE\n")
sink()
saveRDS(cnt, "plot_nll_monotone_4ex_r2_20261102.rds")
