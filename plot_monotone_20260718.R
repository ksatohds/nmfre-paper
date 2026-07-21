# Regenerated on the 2026-07-18 core (row-centering of U removed).  NEW file.
# Fixed-lambda monotone-descent trace for the FOUR empirical examples of the paper.
# Produces the 2x2 panel figure  Figure_objective_monotone.pdf  (Appendix A, fig:obj_monotone).
#
# NEW file (2026-07-16); does NOT modify plot_monotone_20260626_fullEM.R.
# WHY IT EXISTS -- the previous version disagreed with the paper on two points:
#   (1) OD-HUB was fitted with FOUR hubs at Q=4 (hub.index = 13,23,27,40), and its panel was
#       labelled "(c) OD-HUB (Q=4, non-negative)".  The paper uses THREE hubs at Q=3
#       (Tokyo 13 / Osaka 27 / Fukuoka 40, P=44) and Appendix E explicitly REJECTS the
#       four-hub Q=4 specification as producing a redundant basis.  A figure labelled Q=4
#       therefore contradicted the text it sits next to.
#   (2) the fixed penalties were lambda = 1.00 / 0.0104 / 6.75 / 3.30, which are not the
#       lambda values the paper reports.  Table tab:example_summary (main text) gives
#       lambda = 0.0401 (Orthodont) / 0.0203 (HS1939) / 0.0164 (OD-HUB) / 0.0019 (Topic)
#       [2026-07-18 refit].
# Both are corrected here.
#
# What the figure shows: monotone descent is a FIXED-penalty property of the block-coordinate
# updates, so lambda = sigma2/tau2 is held fixed at the reported value for each example
# (sigma2_update = tau2_update = FALSE) and the working objective
#   ||Y - X(Theta A + U)||^2 + lambda ||U||^2
# is traced.  This is the diagnostic Reviewer 2 asks for ("empirical evidence that the stated
# objective function decreases monotonically ... after each block update") and is the same
# fixed-variance setting reported in Section 5 of the main text.
suppressMessages({
  library(nlme); library(lavaan); library(nmfkc); library(readxl)
  library(quanteda); library(quanteda.textmodels)
})
source("nmfre_core20260718.R")

fit_trace <- function(Y, A, Q, theta_sign, x_update, lambda, tag) {
  # Self-initialize (NNDSVD-type) so the descent from a non-optimal start is visible.
  fit <- nmf.re(Y, A, Q = Q,
                theta_sign = theta_sign, x_update = x_update,
                dfU_control = "off",
                sigma2 = lambda, tau2 = 1,
                sigma2_update = FALSE, tau2_update = FALSE,
                inference = FALSE, maxit = 2000, epsilon = 1e-10,
                verbose = FALSE, post_message = FALSE, save_obj_trace = TRUE)
  tr <- fit$obj_trace
  tr <- tr[is.finite(tr)]
  inc <- if (length(tr) > 1) max(diff(tr)) else NA_real_          # >0 would mean an increase
  rel <- if (length(tr) > 1) inc / abs(tr[1]) else NA_real_
  cat(sprintf("%-10s Q=%d lambda=%-7g iters=%4d  obj0=%.6g  objN=%.6g  max step-increase=%.3e (rel %.3e)\n",
              tag, Q, lambda, length(tr), tr[1], tr[length(tr)], inc, rel))
  list(tr = tr, tag = tag, maxinc = inc, relinc = rel)
}

## ---------- Orthodont (Q=1, sign-free), lambda = 0.0401 ----------
d <- Orthodont
subjects <- unique(d$Subject); ages <- sort(unique(d$age))
N <- length(subjects); P <- length(ages)
Y <- matrix(d$distance, P, N, dimnames = list(paste0("age", ages), subjects))
male <- as.vector(matrix(ifelse(d$Sex == "Male", 1, 0), P, N)[1, ])
A <- rbind(intercept = rep(1, N), male = male)
ortho <- fit_trace(Y, A, 1, "free", "seminmf", 0.0401, "Orthodont")

## ---------- Holzinger-Swineford 1939 (Q=3, sign-free), lambda = 0.0203 ----------
data(HolzingerSwineford1939)
d <- HolzingerSwineford1939; d <- d[complete.cases(d), ]
d$age.rev <- -(d$ageyr + d$agemo / 12); d$sex.2 <- ifelse(d$sex == 2, 1, 0)
d$school.GW <- ifelse(d$school == "Grant-White", 1, 0)
d[, c("id", "sex", "ageyr", "agemo", "school", "grade")] <- NULL
d <- nmfkc::nmfkc.normalize(d)
exo <- c("age.rev", "sex.2", "school.GW"); endo <- setdiff(colnames(d), exo)
Y <- t(d[, endo]); A <- t(d[, exo])
hs <- fit_trace(Y, A, 3, "free", "seminmf", 0.0203, "HS1939")

## ---------- OD-HUB (Q=3, THREE hubs, non-negative variant), lambda = 0.0164 ----------
dd <- as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx", sheet = 1, skip = 2))
dd <- dd[dd[, 1] == "1" & dd[, 3] == "02" & dd[, 5] != "100" & dd[, 7] != "100", ]
pref <- unique(dd[, 6]); M <- matrix(NA, 47, 47); colnames(M) <- rownames(M) <- pref
dd[, 5] <- as.numeric(dd[, 5]); dd[, 7] <- as.numeric(dd[, 7])
for (i in 1:47) for (j in 1:47) M[i, j] <- dd[which(dd[, 5] == i & dd[, 7] == j), 9]
M <- log(1 + M)
hub.index <- c(13, 27, 40)          # Tokyo / Osaka / Fukuoka -- the paper's 3-hub spec (P=44)
stopifnot(nrow(M[-hub.index, ]) == 44, nrow(M[hub.index, ]) == 3)
od <- fit_trace(M[-hub.index, ], M[hub.index, ], 3, "nonneg", "mu", 0.0164, "OD-HUB")

## ---------- Topic model (Q=3, non-negative variant), lambda = 0.0019 ----------
corp <- corpus(data_corpus_inaugural); tok <- tokens(corp)
tok <- tokens_remove(tok, pattern = stopwords("en", source = "snowball"))
dfm0 <- dfm(tok); dfm0 <- dfm_select(dfm0, min_nchar = 3); dfm0 <- dfm_trim(dfm0, min_termfreq = 100)
dm <- as.matrix(dfm0); dm <- dm[, order(colSums(dm), decreasing = TRUE)]
Ytp <- t(dm)
Atp <- rbind(pre1865   = ifelse(corp$Year <= 1865, 1, 0),
             y1865_45  = ifelse(corp$Year > 1865 & corp$Year <= 1945, 1, 0),
             after1945 = ifelse(corp$Year > 1945, 1, 0))
colnames(Atp) <- colnames(Ytp)
topic <- fit_trace(Ytp, Atp, 3, "nonneg", "mu", 0.0019, "Topic")

## ---------- 2x2 panel figure ----------
panels <- list(ortho, hs, od, topic)
labs <- c("(a) Orthodont (Q=1, sign-free)",
          "(b) Holzinger--Swineford (Q=3, sign-free)",
          "(c) OD-HUB (Q=3, non-negative)",
          "(d) Topic model (Q=3, non-negative)")

# Annotation: report the increase RELATIVE to the initial objective.  The absolute
# figure is misleading across panels whose objectives differ by three orders of
# magnitude (e.g. HS1939's +0.05 is 2e-4 of its own scale, Topic's +0.06 is 2e-6 of its).
ann <- function(p) {
  if (!is.finite(p$maxinc) || p$maxinc <= 0) "monotone: no increase"
  else sprintf("max relative increase = %.1e", p$relinc)
}

draw <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 4.6, 2.6, 1.2), mgp = c(2.6, 0.8, 0), cex.lab = 1.05)
  for (k in seq_along(panels)) {
    tr <- panels[[k]]$tr
    plot(seq_along(tr), tr, type = "l", lwd = 2, col = "black",   # monochrome, as Figures 1-4
         xlab = "Update (iteration)", ylab = "Objective",
         main = labs[k])
    mtext(ann(panels[[k]]), side = 3, line = -1.2, cex = 0.7, col = "gray30")
  }
}

pdf("Figure_objective_monotone_20260718.pdf", width = 9, height = 7); draw(); dev.off()
png("Figure_objective_monotone_20260718.png", width = 1500, height = 1150, res = 150); draw(); dev.off()
cat("\nWrote Figure_objective_monotone_20260718.pdf / .png\n")
cat("Deploy with:  cp Figure_objective_monotone_20260718.pdf ../Figure_objective_monotone.pdf\n")
