# Table 1 of the main text (Q = 3 simulation) without the "Point region" column.  NEW file, 2026-11-08 (work memo A165).
# Author: "Point regionの列は本当に必要ですか？ CoverageとRejection at true zeros Power まで示していますよ？"
# The share of point regions stays in Appendix E (tab_app_sim_fits.tex, same values) and in the text of Section 5 (3).
# No recomputation: reads emit_tables_sim_v5_r2_20261102.rds ($main), written by emit_tables_sim_v5_r2_20261102.R,
# whose tab_sim_q3.tex is left untouched.  The script checks that the remaining nine columns match tab_sim_q3.tex.
# Usage: Rscript emit_tab_sim_q3_nopoint_20261108.R     Output: tab_sim_q3_nopoint_20261108.tex (+ .log via sink)
sink("emit_tab_sim_q3_nopoint_20261108.log", split = TRUE)
M <- readRDS("emit_tables_sim_v5_r2_20261102.rds")$main
f3 <- function(x) formatC(x, format = "f", digits = 3)
f2 <- function(x) formatC(x, format = "f", digits = 2)
lab_err <- c(gaussian = "Equal", varhet = "Unequal")
B <- M[order(match(M$err, c("gaussian", "varhet")), M$eps), ]
rows <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                lab_err[B$err], f2(B$eps),
                f3(B$cov_wald), f3(B$cov_wild), f3(B$cov_prof),
                f3(B$size0_wald), f3(B$size0_wild), f3(B$size0_prof),
                f3(B$power_wald), f3(B$power_wild), f3(B$power_prof))
rows <- append(rows, "\\midrule", after = 4)
head <- c("\\begin{tabular}{llrrrrrrrrr}", "\\toprule",
  "& & \\multicolumn{3}{c}{Coverage} & \\multicolumn{3}{c}{Rejection at true zeros} & \\multicolumn{3}{c}{Power} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
  "Errors & $\\varepsilon$ & Wald & Wild & Profile & Wald & Wild & Profile & Wald & Wild & Profile \\\\", "\\midrule")
out <- c(head, rows, "\\bottomrule", "\\end{tabular}")
writeLines(out, "tab_sim_q3_nopoint_20261108.tex")
# check against the existing table: drop its third column ("Point region") and compare the data rows
old <- readLines("tab_sim_q3.tex")
dat_old <- grep("^(Equal|Unequal) ", old, value = TRUE)
dat_new <- grep("^(Equal|Unequal) ", out, value = TRUE)
drop3 <- function(s) { p <- strsplit(sub(" \\\\\\\\$", "", s), " & ", fixed = TRUE)[[1]]; paste(p[-3], collapse = " & ") }
strip <- function(s) sub(" \\\\\\\\$", "", s)
same <- identical(unname(vapply(dat_old, drop3, "")), unname(vapply(dat_new, strip, "")))
cat(sprintf("rows %d | identical to tab_sim_q3.tex without the Point region column: %s\n", length(dat_new), same))
cat(out, sep = "\n")
sink()
