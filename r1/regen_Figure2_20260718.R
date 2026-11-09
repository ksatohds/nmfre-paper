# Regenerated on the 2026-07-18 core (row-centering of U removed).  NEW file.
# Regenerate Figure 2 (HS1939, Q=3) with SIGNED A->X (covariate->factor) edges.
# The stock nmfkc.DOT draws only val >= threshold (positive), so under the
# sign-free Theta the NEGATIVE covariate effects (the key result) are dropped.
# This NEW script defines a signed variant nmfkc.DOT.signed (does NOT modify
# nmfkc.R): A->X edges are drawn for |C| >= threshold, negatives in red/dashed
# with the signed label, magnitude-scaled penwidth (consistent max=max|C|).
# Reuses the package DOT helpers (.nmfkc_dot_*) sourced from nmfkc.R.
suppressPackageStartupMessages({
  library(lavaan); library(nmfkc); library(DiagrammeR); library(DiagrammeRsvg); library(rsvg)
})
# The DOT helpers (.nmfkc_dot_*) are internal (non-exported) functions of the
# nmfkc package; import them from the installed package namespace (install with
# remotes::install_github("ksatohds/nmfkc@develop")) so the signed variant below
# can reuse them without sourcing the package source directly.
local({
  ns <- asNamespace("nmfkc")
  for (fn in c(".nmfkc_dot_sanitize_id", ".nmfkc_dot_header",
               ".nmfkc_dot_cluster_nodes", ".nmfkc_dot_penwidth",
               ".nmfkc_dot_format_coef", ".nmfkc_dot_digits_from_threshold"))
    assign(fn, get(fn, envir = ns), envir = globalenv())
})
source("nmfre_core20260718.R")

## ---- signed YXA DOT (helpers come from nmfkc.R) -------------------------
nmfkc.DOT.signed <- function(x, threshold = 0.05, rankdir = "TB", fill = FALSE,
                             weight_scale = 5,
                             Y.title = "Observation (Y)", X.title = "Factor (X)",
                             A.title = "Covariate (A)",
                             neg.color = "black", neg.style = "dashed") {
  X <- x$X; B <- x$B; C <- as.matrix(x$C)
  P <- nrow(X); Q <- ncol(X); A_cols <- ncol(C)
  Y_labels <- rownames(X); X_labels <- colnames(X); A_labels <- colnames(C)
  if (is.null(Y_labels)) Y_labels <- paste0("Y", seq_len(P))
  if (is.null(X_labels)) X_labels <- paste0("Factor", seq_len(Q))
  if (is.null(A_labels)) A_labels <- paste0("A", seq_len(A_cols))
  Y_ids <- .nmfkc_dot_sanitize_id(paste0("Y_", seq_len(P)))
  X_ids <- .nmfkc_dot_sanitize_id(paste0("X_", seq_len(Q)))
  A_ids <- .nmfkc_dot_sanitize_id(paste0("A_", seq_len(A_cols)))

  scr <- .nmfkc_dot_header(graph_name = "NMF", rankdir = rankdir)
  scr <- paste0(scr, "\n  // Output variables (Y)\n",
    .nmfkc_dot_cluster_nodes("Y", Y.title, Y_ids, Y_labels, shape = "box",
                             fill = fill, fillcolor = "lightblue", line_width = 1.5))
  scr <- paste0(scr, "\n  // Latent factors (X)\n",
    .nmfkc_dot_cluster_nodes("X", X.title, X_ids, X_labels, shape = "ellipse",
                             fill = fill, fillcolor = "wheat", line_width = 1.0))
  scr <- paste0(scr, "\n  // Covariates (A)\n",
    .nmfkc_dot_cluster_nodes("A", A.title, A_ids, A_labels, shape = "box",
                             fill = fill, fillcolor = "lightcoral", line_width = 1.5))
  scr <- paste0(scr, '\n  edge [fontname="Meiryo", fontsize=8, arrowhead=open];\n')

  pw <- .nmfkc_dot_penwidth
  digits <- .nmfkc_dot_digits_from_threshold(threshold)
  fmtc <- function(v) .nmfkc_dot_format_coef(v, digits)

  ## A -> X (C), SIGNED: |C| >= threshold; negatives red/dashed, scaled by max|C|
  scr <- paste0(scr, "\n  // A -> X edges (C, signed)\n")
  max_C <- suppressWarnings(max(abs(C), na.rm = TRUE))
  if (is.finite(max_C) && max_C > 0) {
    for (q in seq_len(Q)) for (k in seq_len(A_cols)) {
      val <- C[q, k]
      if (is.finite(val) && abs(val) >= threshold) {
        pen <- pw(abs(val), max_C, weight_scale)
        col <- if (val < 0) neg.color else "black"
        sty <- if (val < 0) neg.style else "solid"
        scr <- paste0(scr, sprintf('  %s -> %s [label="%s", penwidth=%.2f, color=%s, fontcolor=%s, style=%s];\n',
                                   A_ids[k], X_ids[q], fmtc(val), pen, col, col, sty))
      }
    }
  }
  ## X -> Y (X >= 0 always)
  scr <- paste0(scr, '\n  // X -> Y edges (X)\n  edge [color="gray0", fontcolor="gray0", style=solid];\n')
  max_X <- suppressWarnings(max(X, na.rm = TRUE))
  if (is.finite(max_X) && max_X > 0) {
    for (i in seq_len(P)) for (j in seq_len(Q)) {
      val <- X[i, j]
      if (is.finite(val) && val >= threshold) {
        scr <- paste0(scr, sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                                   X_ids[j], Y_ids[i], fmtc(val), pw(val, max_X, weight_scale)))
      }
    }
  }
  paste0(scr, "}\n")
}

## ---- HS data + fit (identical to regen_Figure2_newcore_fullEM.R) --------
data(HolzingerSwineford1939); d <- HolzingerSwineford1939; d <- d[complete.cases(d), ]
d$age.rev <- -(d$ageyr + d$agemo/12); d$sex.2 <- ifelse(d$sex==2,1,0)
d$school.GW <- ifelse(d$school=="Grant-White",1,0)
d[, c("id","sex","ageyr","agemo","school","grade")] <- NULL
d <- nmfkc::nmfkc.normalize(d)
exo <- c("age.rev","sex.2","school.GW"); endo <- setdiff(colnames(d), exo)
Y <- t(d[, endo]); A <- t(d[, exo]); Q0 <- 3
ini <- nmfkc::nmfkc(Y, A, Q = Q0, epsilon = 1e-8, method = "EU", seed = 1)  # match refit_empirical_fullEM.R
fit <- nmf.re(Y, A, Q = Q0, X.ini = ini$X, Theta.ini = ini$C, prefix = "Factor",
              theta_sign = "free", x_update = "seminmf", dfU_control = "off",
              sigma2 = 1, tau2 = 1, sigma2_update = TRUE, tau2_update = TRUE,
              inference = FALSE, verbose = FALSE, post_message = FALSE)
rownames(fit$X) <- endo
if (is.null(colnames(fit$X))) colnames(fit$X) <- paste0("Factor", seq_len(Q0))
rownames(fit$C) <- colnames(fit$X); colnames(fit$C) <- exo
cat("=== Theta (=C): A -> X coefficients (signs matter) ===\n"); print(round(fit$C, 3))
cat("negative |C|>=0.05 paths:\n")
for (q in seq_len(Q0)) for (k in seq_len(ncol(fit$C))) {
  v <- fit$C[q,k]; if (v < -0.05) cat(sprintf("  %s -> %s : %.3f\n", colnames(fit$C)[k], rownames(fit$C)[q], v))
}

## ---- build + export ------------------------------------------------------
dot <- nmfkc.DOT.signed(fit, threshold = 0.05, rankdir = "TB", fill = FALSE,
                        Y.title = "Observation(Y)", X.title = "Factor(X)", A.title = "Covariate(A)")
svg_txt <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
rsvg::rsvg_pdf(charToRaw(svg_txt), file = "Figure2_signed_fullEM.pdf")
rsvg::rsvg_png(charToRaw(svg_txt), file = "Figure2_signed_preview.png", width = 1100, height = 1100)
cat("\nWROTE Figure2_signed_fullEM.pdf and Figure2_signed_preview.png\n")
cat("n negative A->X edges in DOT:", lengths(regmatches(dot, gregexpr("color=red", dot))), "\n")
