# Holzinger-Swineford (Section 6.2) for R2 from the final fit (refit_examples_r2_20261102.rds, C81: nmfkc a343c39,
# intercept, sign-free, Q = 3, geomin orientation; identification region a point): LaTeX rows of the basis table
# (tab:hs1939_basis) and of the covariate-effect table (tab:hs1939_theta, R1 format), the intercepts (estimates only,
# A24), the profile-likelihood intervals, and Figure 2 (signed path diagram, drawn as regen_Figure2_20260718.R:
# A -> X edges for |Theta| >= 0.05, negative dashed; the intercept is not drawn).  (NEW file, 2026-11-02; work memo C89.)
# Reads RDS only (no refit).  Usage: Rscript make_hs_r2_tables_figure_20261102.R
# Output: make_hs_r2_tables_figure_20261102.{log,rds}, Figure2_r2_20261102.{pdf,png} (the manuscript's Figure2.pdf is
# not touched)
suppressPackageStartupMessages({library(DiagrammeR); library(DiagrammeRsvg); library(rsvg)})
z <- readRDS("refit_examples_r2_20261102.rds")$hs
X <- z$fit$X; C <- z$fit$C; cf <- z$inf$coefficients
cat(sprintf("HS fit: NLL %.4f | identification region max width %g | sigma2 %.5f tau2 %.4f lambda %.4f dfU %.2f r %.3f\n",
            z$nll, max(z$fit$C.ident.upper - z$fit$C.ident.lower), z$fit$sigma2, z$fit$tau2, z$fit$lambda, z$fit$dfU,
            z$fit$dfU.frac))
tests <- paste0("x", 1:9); rownames(X) <- tests
desc <- c("Visual perception", "Cubes", "Lozenges", "Paragraph comprehension", "Sentence completion", "Word meaning",
          "Speeded addition", "Speeded counting", "Speeded discrimination")
dom <- apply(X, 1, which.max)
cat("dominant component per test:", paste(sprintf("%s=%d", tests, dom), collapse = " "), "\n")
stopifnot(identical(unname(dom), c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)))     # visual / verbal / speed
fac <- paste0("Factor", 1:3); colnames(X) <- fac; rownames(C) <- fac
cat("\n%% tab:hs1939_basis (body)\n")
rows_basis <- sapply(1:9, function(i) {
  v <- sprintf("%.3f", X[i, ]); v[dom[i]] <- sprintf("\\textbf{%s}", v[dom[i]])
  sprintf("$x_%d$ & %-23s & %s \\\\", i, desc[i], paste(v, collapse = " & "))
})
cat(rows_basis, sep = "\n")
cat(sprintf("exact zeros in X: %s\n", paste(which(X < 1e-8, arr.ind = FALSE), collapse = " ")))
print(round(X, 4))
## covariate effects (R1 format: Covariate, Basis, Estimate, SE, BSE, z, p two-sided)
cf$Basis <- sub("^Basis", "Factor", cf$Basis)
cf$prof_lo <- as.vector(z$inf$C.profile.lower); cf$prof_hi <- as.vector(z$inf$C.profile.upper)   # vec order = rows of cf
cv <- cf[cf$Covariate != "intercept", ]
pv <- function(p) ifelse(p < 0.001, "$<0.001$", ifelse(p < 0.1, sprintf("%.4f", p), sprintf("%.3f", p)))
rows_theta <- sprintf("%-9s & %s & %6.3f & %.3f & %.3f & %5.2f & %s \\\\", cv$Covariate, cv$Basis, cv$Estimate, cv$SE,
                      cv$BSE, cv$z_value, pv(cv$p_value))
cat("\n%% tab:hs1939_theta (body)\n"); cat(rows_theta, sep = "\n")
cat("\nintercepts (estimates only, A24):", paste(sprintf("%s %.3f", cf$Basis[cf$Covariate == "intercept"],
                                                     cf$Estimate[cf$Covariate == "intercept"]), collapse = ", "), "\n")
iv <- data.frame(coef = paste(cv$Covariate, cv$Basis), est = cv$Estimate, wild_lo = cv$CI_low, wild_hi = cv$CI_high,
                 prof_lo = cv$prof_lo, prof_hi = cv$prof_hi)
cat("\nintervals (wild-bootstrap percentile; profile likelihood):\n"); print(iv, row.names = FALSE, digits = 3)

## ---- Figure 2: signed DOT (as regen_Figure2_20260718.R, helpers from the nmfkc namespace) ----
ns <- asNamespace("nmfkc")
for (h in c(".nmfkc_dot_cluster_nodes", ".nmfkc_dot_digits_from_threshold", ".nmfkc_dot_format_coef",
            ".nmfkc_dot_header", ".nmfkc_dot_penwidth", ".nmfkc_dot_sanitize_id")) assign(h, get(h, envir = ns))
nmfkc.DOT.signed <- function(X, C, threshold = 0.05, rankdir = "TB", fill = FALSE, weight_scale = 5,
                             Y.title = "Observation(Y)", X.title = "Factor(X)", A.title = "Covariate(A)",
                             neg.color = "black", neg.style = "dashed") {
  P <- nrow(X); Q <- ncol(X); K <- ncol(C)
  Y_ids <- .nmfkc_dot_sanitize_id(paste0("Y_", seq_len(P))); X_ids <- .nmfkc_dot_sanitize_id(paste0("X_", seq_len(Q)))
  A_ids <- .nmfkc_dot_sanitize_id(paste0("A_", seq_len(K)))
  scr <- .nmfkc_dot_header(graph_name = "NMF", rankdir = rankdir)
  scr <- paste0(scr, "\n", .nmfkc_dot_cluster_nodes("Y", Y.title, Y_ids, rownames(X), shape = "box", fill = fill,
                                                    fillcolor = "lightblue", line_width = 1.5))
  scr <- paste0(scr, "\n", .nmfkc_dot_cluster_nodes("X", X.title, X_ids, colnames(X), shape = "ellipse", fill = fill,
                                                    fillcolor = "wheat", line_width = 1.0))
  scr <- paste0(scr, "\n", .nmfkc_dot_cluster_nodes("A", A.title, A_ids, colnames(C), shape = "box", fill = fill,
                                                    fillcolor = "lightcoral", line_width = 1.5))
  scr <- paste0(scr, '\n  edge [fontname="Meiryo", fontsize=8, arrowhead=open];\n')
  digits <- .nmfkc_dot_digits_from_threshold(threshold); fmtc <- function(v) .nmfkc_dot_format_coef(v, digits)
  max_C <- max(abs(C))
  for (q in seq_len(Q)) for (k in seq_len(K)) {
    val <- C[q, k]
    if (abs(val) >= threshold) {
      col <- if (val < 0) neg.color else "black"; sty <- if (val < 0) neg.style else "solid"
      scr <- paste0(scr, sprintf('  %s -> %s [label="%s", penwidth=%.2f, color=%s, fontcolor=%s, style=%s];\n', A_ids[k],
                                 X_ids[q], fmtc(val), .nmfkc_dot_penwidth(abs(val), max_C, weight_scale), col, col, sty))
    }
  }
  scr <- paste0(scr, '\n  edge [color="gray0", fontcolor="gray0", style=solid];\n')
  max_X <- max(X)
  for (i in seq_len(P)) for (j in seq_len(Q)) if (X[i, j] >= threshold)
    scr <- paste0(scr, sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n', X_ids[j], Y_ids[i], fmtc(X[i, j]),
                               .nmfkc_dot_penwidth(X[i, j], max_X, weight_scale)))
  paste0(scr, "}\n")
}
Cd <- C[, colnames(C) != "intercept", drop = FALSE]
dot <- nmfkc.DOT.signed(X, Cd)
svg_txt <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
rsvg::rsvg_pdf(charToRaw(svg_txt), file = "Figure2_r2_20261102.pdf")
rsvg::rsvg_png(charToRaw(svg_txt), file = "Figure2_r2_20261102.png", width = 1100, height = 1100)
cat(sprintf("\nFigure 2: %d A->X edges (|Theta| >= 0.05), %d negative (dashed); %d X->Y edges\n",
            sum(abs(Cd) >= 0.05), sum(Cd <= -0.05), sum(X >= 0.05)))
saveRDS(list(X = X, C = C, coefficients = cf, rows_basis = rows_basis, rows_theta = rows_theta, intervals = iv,
             dot = dot), "make_hs_r2_tables_figure_20261102.rds")
cat(sprintf("WROTE Figure2_r2_20261102.pdf/.png\nDONE %s\n", format(Sys.time())))
