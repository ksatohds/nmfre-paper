# LaTeX table bodies for the empirical examples of the R2 manuscript, written to .tex files to be \input
# (author: tables are not typed by hand; work memo A45, A91).  NEW file, 2026-11-02.
# Reads refit_examples_r2_20261102.rds (C81: the five examples refitted with the nmfkc default of a343c39) and,
# for the component labels of the topic model and OD-HUB as in R1, compare_examples_r2_20261102.rds (C82).
# Reads RDS only; no fitting, no OD data.  The rows agree with the drafts of 改訂文案20261031.md §3.7, §3.12, §3.16.
# Usage: Rscript emit_tab_examples_r2_20261102.R
# Output: emit_tab_examples_r2_20261102.{log,rds} and tab_ex_summary.tex, tab_hs_theta.tex, tab_hs_basis.tex,
#         tab_college_theta.tex, tab_topic_theta.tex, tab_odhub_theta.tex
new <- readRDS("refit_examples_r2_20261102.rds"); cmp <- readRDS("compare_examples_r2_20261102.rds")
pv <- function(p) ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p))
row7 <- function(cf, lab1, lab2) {
  sprintf("%s & %s & %7.3f & %6.3f & %6.3f & %6.2f & %s \\\\", lab1, lab2, cf$Estimate, cf$SE, cf$BSE,
          cf$z_value, pv(cf$p_value))
}
coefs <- function(nm) { cf <- new[[nm]]$inf$coefficients; cf[cf$Covariate != "intercept", ] }
pick <- function(cf, cov, basis) cf[match(paste(cov, basis), paste(cf$Covariate, cf$Basis)), ]
out <- list()
sink("emit_tab_examples_r2_20261102.log", split = TRUE)
cat("emit_tab_examples_r2_20261102 |", format(Sys.time()), "\n")

## HS (Q = 3; Basis1 = visual, Basis2 = textual, Basis3 = speed, as in C81/C89)
cf <- coefs("hs"); cv <- rep(c("age.rev", "sex.2", "school.GW"), each = 3); bs <- rep(sprintf("Basis%d", 1:3), 3)
r <- pick(cf, cv, bs); stopifnot(!anyNA(r$Estimate))
tx <- row7(r, sprintf("%-9s", cv), rep(sprintf("Factor%d", 1:3), 3))
writeLines(tx, "tab_hs_theta.tex"); out$hs_theta <- r
cat("\n-- tab_hs_theta.tex\n"); cat(tx, sep = "\n")
ic <- new$hs$inf$coefficients; ic <- ic[ic$Covariate == "intercept", ]
cat(sprintf("HS intercepts: %s\n", paste(sprintf("%.3f", ic$Estimate), collapse = ", ")))
X <- new$hs$fit$X; stopifnot(nrow(X) == 9, ncol(X) == 3)
desc <- c("Visual perception", "Cubes", "Lozenges", "Paragraph comprehension", "Sentence completion", "Word meaning",
          "Speeded addition", "Speeded counting", "Speeded discrimination")
tx <- sapply(1:9, function(i) {
  v <- sprintf("%.3f", X[i, ]); j <- which.max(X[i, ]); v[j] <- sprintf("\\textbf{%s}", v[j])
  sprintf("$x_%d$ & %-23s & %s \\\\", i, desc[i], paste(v, collapse = " & "))
})
writeLines(tx, "tab_hs_basis.tex"); out$hs_basis <- X
cat("\n-- tab_hs_basis.tex\n"); cat(tx, sep = "\n")
cat(sprintf("HS exact zeros of X (row, column): %s\n",
            paste(apply(which(X == 0, arr.ind = TRUE), 1, function(z) sprintf("(x%d, %d)", z[1], z[2])), collapse = " ")))

## College (Q = 2; Basis1 = faculty qualifications, Basis2 = selectivity and resources, as in C40/C81)
cf <- coefs("college")
cvn <- c("Apps", "Accept", "Enroll", "F.Undergrad", "Room.Board", "S.F.Ratio")
cvl <- c("Apps", "Accept", "Enroll", "F.Undergrad", "Room.Board", "S.F.Ratio (reversed)")
r <- pick(cf, rep(cvn, each = 2), rep(c("Basis1", "Basis2"), 6)); stopifnot(!anyNA(r$Estimate))
tx <- row7(r, sprintf("%-20s", rep(cvl, each = 2)), rep(c("Faculty    ", "Selectivity"), 6))
writeLines(tx, "tab_college_theta.tex"); out$college_theta <- r
cat("\n-- tab_college_theta.tex\n"); cat(tx, sep = "\n")
ic <- new$college$inf$coefficients; ic <- ic[ic$Covariate == "intercept", ]
cat(sprintf("College intercepts: %s\n", paste(sprintf("%.3f", ic$Estimate), collapse = ", ")))
Xc <- new$college$fit$X; sh <- Xc / rowSums(Xc)
cat("College row shares of X (Basis1, Basis2):\n"); print(round(sh, 3))

## topic model and OD-HUB (non-negative variant; components in the R1 order, as emit_tables_examples_r2_20261102.R)
rows_r1 <- function(nm, key) {
  z <- new[[nm]]; p <- cmp[[key]]$perm; Q <- ncol(z$fit$X); K <- ncol(z$fit$C)
  z$inf$coefficients[as.vector(matrix(seq_len(Q * K), Q, K)[p, , drop = FALSE]), ]
}
tc <- rows_r1("topic", "topic C55 (mu, draft 3.7)")
era <- rep(c("Early (1775--1865)     ", "Industrial (1865--1945)", "Postwar (1945--)       "), each = 3)
tx <- row7(tc, era, rep(sprintf("Topic %d", 1:3), 3))
writeLines(tx, "tab_topic_theta.tex"); out$topic_theta <- tc
cat("\n-- tab_topic_theta.tex (one-sided p)\n"); cat(tx, sep = "\n")
oc <- rows_r1("odhub", "odhub C53 (mu, R1 table)")
hub <- rep(c("Tokyo  ", "Osaka  ", "Fukuoka"), each = 3)
tx <- row7(oc, hub, rep(sprintf("Region%d", 1:3), 3))
writeLines(tx, "tab_odhub_theta.tex"); out$odhub_theta <- oc
cat("\n-- tab_odhub_theta.tex (one-sided p)\n"); cat(tx, sep = "\n")

## summary table: N, P, Q, sigma2, tau2, lambda, r = dfU/(NQ); main-text examples first, then the appendix ones
lab <- c(hs = "HS1939", college = "College", topic = "Topic model", orthodont = "Orthodont", odhub = "OD--HUB")
fs <- function(x) if (x < 0.1) sprintf("%.4f", x) else if (x < 10) sprintf("%.3f", x) else if (x < 1000) sprintf("%.2f", x) else sprintf("%.1f", x)
tx <- sapply(names(lab), function(nm) {
  z <- new[[nm]]; f <- z$fit
  sprintf("%-11s & %d & %d & %d & %s & %s & %.4f & %.3f \\\\", lab[[nm]], z$N, z$P, z$Q, fs(f$sigma2), fs(f$tau2),
          f$lambda, f$dfU.frac)
})
tx <- append(tx, "\\midrule", after = 3)
writeLines(tx, "tab_ex_summary.tex"); out$summary <- tx
cat("\n-- tab_ex_summary.tex\n"); cat(tx, sep = "\n")
cat("\nDONE\n")
sink()
saveRDS(out, "emit_tab_examples_r2_20261102.rds")
