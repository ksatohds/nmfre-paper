# LaTeX tables for the Q = 3 simulation (HS design, v5; work memo C79, C88, C92, A88).  NEW file, 2026-11-02.
# Author: "シミュレーションの表も作って".  No recomputation: reads the saved summaries
#   summarize_sim_r2_hs_v5_20261102.rds (wild, profile, fits, Xi vs ESEM),
#   summarize_sim_r2_hs_v5_extra_20261102.rds (region width, variances),
#   summarize_sim_r2_hs_v5_wald_20261102.rds (Wald),
# and, for the known-anchor arm only, the per-replicate files out_hs_v5/gaussian_eps_0.00/rep_*.rds.
# Writes tabular bodies to be \input inside table floats (captions are written in the manuscript):
#   tab_sim_q3.tex        main text: share of point regions; coverage, size at the true zeros, power (Wald, wild, profile)
#   tab_app_sim_fits.tex  Appendix E: fits (degenerate, point share, widths, congruence, variance bias)
#   tab_app_sim_coef.tex  Appendix E: coverage by coefficient (equal and unequal error variances, e = 0 and 0.2)
#   tab_app_sim_xi.tex    Appendix E: composite effects Xi, NMF-RE vs ESEM (RMSE, coverage, Heywood), known-anchor arm
# Usage: Rscript emit_tables_sim_v5_r2_20261102.R      Output: emit_tables_sim_v5_r2_20261102.{log,rds} and the four .tex
s5 <- readRDS("summarize_sim_r2_hs_v5_20261102.rds")
ex <- readRDS("summarize_sim_r2_hs_v5_extra_20261102.rds")
wd <- readRDS("summarize_sim_r2_hs_v5_wald_20261102.rds")
f3 <- function(x) formatC(x, format = "f", digits = 3)
f2 <- function(x) formatC(x, format = "f", digits = 2)
lab_err <- c(gaussian = "Equal", varhet = "Unequal")
ord <- function(df) df[order(match(df$err, c("gaussian", "varhet")), df$eps), ]

## ---- main table ----
B <- ord(wd$by_design)
P <- ord(ex$width)
stopifnot(all(B$err == P$err), all(B$eps == P$eps))
main <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                lab_err[B$err], f2(B$eps), f3(P$point),
                f3(B$cov_wald), f3(B$cov_wild), f3(B$cov_prof),
                f3(B$size0_wald), f3(B$size0_wild), f3(B$size0_prof),
                f3(B$power_wald), f3(B$power_wild), f3(B$power_prof))
main <- append(main, "\\midrule", after = 4)
head_main <- c("\\begin{tabular}{llrrrrrrrrrr}", "\\toprule",
  "& & Point & \\multicolumn{3}{c}{Coverage} & \\multicolumn{3}{c}{Rejection at true zeros} & \\multicolumn{3}{c}{Power} \\\\",
  "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}\\cmidrule(lr){10-12}",
  "Errors & $\\varepsilon$ & region & Wald & Wild & Profile & Wald & Wild & Profile & Wald & Wild & Profile \\\\", "\\midrule")
writeLines(c(head_main, main, "\\bottomrule", "\\end{tabular}"), "tab_sim_q3.tex")

## ---- fits (Appendix E) ----
I <- s5$I
fits <- do.call(rbind, lapply(split(I, list(I$err, I$eps), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], n = nrow(d), degenerate = sum(d$degenerate), extra = sum(d$extra > 0),
  cong_med = median(d$cong_min), cong_p05 = unname(quantile(d$cong_min, 0.05)))))
fits <- ord(fits); Wx <- ord(ex$width)
stopifnot(all(fits$err == Wx$err), all(fits$eps == Wx$eps))
fr <- sprintf("%s & %s & %d & %d & %s & %s & %s & %s & %s & %s & %s \\\\",
              lab_err[fits$err], f2(fits$eps), fits$n, fits$degenerate, f3(Wx$point),
              f3(Wx$coefwidth_med_nonpoint), f3(Wx$mean_se), f3(fits$cong_med), f3(fits$cong_p05),
              ifelse(fits$err == "varhet", "--", f3(Wx$sigma2_relbias)), f3(Wx$tau2_relbias))
# (unequal error variances: there is no single true sigma^2, so its relative bias is not defined)
fr <- append(fr, "\\midrule", after = 4)
writeLines(c("\\begin{tabular}{llrrrrrrrrr}", "\\toprule",
  "& & & & Point & Median & Mean & \\multicolumn{2}{c}{Congruence of $\\widehat X$} & \\multicolumn{2}{c}{Relative bias} \\\\",
  "\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}",
  "Errors & $\\varepsilon$ & $R$ & Degenerate & region & width & SE & median & 5\\% & $\\widehat\\sigma^2$ & $\\widehat\\tau^2$ \\\\", "\\midrule",
  fr, "\\bottomrule", "\\end{tabular}"), "tab_app_sim_fits.tex")

## ---- coverage by coefficient (Appendix E) ----
C <- wd$by_coef
pretty <- function(cf) { p <- strsplit(cf, ":B"); sapply(p, function(z) sprintf("\\texttt{%s} $\\to$ %s", z[1], z[2])) }
coef_rows <- function(err) {
  a <- C[C$err == err & C$eps == 0, ]; b <- C[C$err == err & C$eps == 0.2, ]
  a <- a[order(a$coef), ]; b <- b[order(b$coef), ]; stopifnot(all(a$coef == b$coef))
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", pretty(a$coef), f3(a$target),
          f3(a$cov_wald), f3(a$cov_wild), f3(a$cov_prof), f3(b$cov_wald), f3(b$cov_wild), f3(b$cov_prof))
}
writeLines(c("\\begin{tabular}{lrrrrrrr}", "\\toprule",
  "& & \\multicolumn{3}{c}{$\\varepsilon=0$} & \\multicolumn{3}{c}{$\\varepsilon=0.2$} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
  "Effect & True value & Wald & Wild & Profile & Wald & Wild & Profile \\\\", "\\midrule",
  "\\multicolumn{8}{l}{\\emph{Equal error variances}} \\\\", coef_rows("gaussian"), "\\midrule",
  "\\multicolumn{8}{l}{\\emph{Unequal error variances}} \\\\", coef_rows("varhet"),
  "\\bottomrule", "\\end{tabular}"), "tab_app_sim_coef.tex")

## ---- Xi vs ESEM, and the known-anchor arm (Appendix E) ----
xi <- do.call(rbind, lapply(split(I, list(I$err, I$eps), drop = TRUE), function(d) data.frame(
  err = d$err[1], eps = d$eps[1], rmse = mean(d$xi_rmse), cov = mean(d$xi_cov),
  rmse_esem = mean(d$xi_rmse_esem[d$esem_ok]), cov_esem = mean(d$xi_cov_esem[d$esem_ok]),
  esem_ok = mean(d$esem_ok), heywood = mean(d$heywood))))
xi <- ord(xi)
xr <- sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", lab_err[xi$err], f2(xi$eps), f3(xi$rmse), f3(xi$rmse_esem),
              f3(xi$cov), f3(xi$cov_esem), f3(xi$heywood))
xr <- append(xr, "\\midrule", after = 4)
writeLines(c("\\begin{tabular}{llrrrrr}", "\\toprule",
  "& & \\multicolumn{2}{c}{RMSE of $\\widehat\\Xi$} & \\multicolumn{2}{c}{Coverage of $\\Xi$} & ESEM \\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  "Errors & $\\varepsilon$ & NMF-RE & ESEM & NMF-RE & ESEM & Heywood \\\\", "\\midrule",
  xr, "\\bottomrule", "\\end{tabular}"), "tab_app_sim_xi.tex")

## known-anchor arm (e = 0, equal variances): Wald and wild from the per-replicate tables
z <- qnorm(0.975)
fl <- list.files("out_hs_v5/gaussian_eps_0.00", pattern = "^rep_.*rds$", full.names = TRUE)
An <- do.call(rbind, lapply(fl, function(f) { x <- readRDS(f); if (isTRUE(x$ok) && !is.null(x$anchor$table)) x$anchor$table }))
an <- data.frame(n_rows = nrow(An), cov_wald = mean(abs(An$est - An$target) <= z * An$se), cov_wild = mean(An$cov_wild),
                 size0_wald = mean((abs(An$est) > z * An$se)[An$target == 0]), size0_wild = mean(An$rej0_wild[An$target == 0]),
                 power_wald = mean((abs(An$est) > z * An$se)[An$target != 0]), power_wild = mean(An$rej0_wild[An$target != 0]))

options(width = 200)
cat("== main table (tab_sim_q3.tex) ==\n"); print(cbind(B[, c("err", "eps")], point = P$point, B[, c("cov_wald", "cov_wild", "cov_prof",
  "size0_wald", "size0_wild", "size0_prof", "power_wald", "power_wild", "power_prof")]), row.names = FALSE, digits = 3)
cat("\n== fits (tab_app_sim_fits.tex) ==\n"); print(cbind(fits, Wx[, c("point", "coefwidth_med_nonpoint", "mean_se", "sigma2_relbias", "tau2_relbias")]), row.names = FALSE, digits = 3)
cat("\n== Xi vs ESEM (tab_app_sim_xi.tex) ==\n"); print(xi, row.names = FALSE, digits = 3)
cat("\n== known-anchor arm (e = 0, equal variances) ==\n"); print(an, row.names = FALSE, digits = 3)
for (f in c("tab_sim_q3.tex", "tab_app_sim_fits.tex", "tab_app_sim_coef.tex", "tab_app_sim_xi.tex")) {
  cat("\n----", f, "----\n"); cat(readLines(f), sep = "\n")
}
saveRDS(list(main = cbind(B, point = P$point), fits = cbind(fits, Wx), coef = C, xi = xi, anchor = an),
        "emit_tables_sim_v5_r2_20261102.rds")
cat("\nsaved emit_tables_sim_v5_r2_20261102.rds\n")
