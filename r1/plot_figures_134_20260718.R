# Regenerated on the 2026-07-18 core (row-centering of U removed).  NEW file.
# =====================================================================
# plot_figures_134_20260626.R
# Regenerate Figure 1 (Orthodont), Figure 3 (OD-HUB map), Figure 4 (Topic)
# under the CURRENT outer-inner core (nmfre_core20260626.R), PURE ECM
# (variances estimated, no cap), matching the revised manuscript.
#   Fig1: Orthodont Q=1, sign-free Theta.
#   Fig3: OD-HUB Q=3, hubs Tokyo(13)/Osaka(27)/Fukuoka(40), non-negative variant.
#   Fig4: Inaugural-address topics Q=3, non-negative variant, content-relabelled.
# (Figure 2 has its own generator: regen_Figure2_newcore.R.)
# Run from the repository root (sources the core and reads the OD xlsx here; run download_od.R first).
# Writes *_regen.pdf first for verification; copy to ../FigureN.pdf when confirmed.
# =====================================================================
suppressMessages({library(nlme); library(readxl); library(NipponMap)
  library(quanteda); library(nmfkc)})
source("nmfre_core20260718.R")

ECMfit <- function(Y, A, Q, ts, xu, pfx) {
  r0 <- nmfkc(Y = Y, A = A, Q = Q, epsilon = 1e-8, method = "EU")
  nmf.re(Y, A, Q = Q, X.ini = r0$X, Theta.ini = r0$C, prefix = pfx,
         theta_sign = ts, x_update = xu, dfU_control = "off",
         sigma2 = 1, tau2 = 1, sigma2_update = TRUE, tau2_update = TRUE,
         inference = FALSE, verbose = FALSE, post_message = FALSE)
}
Bprob <- function(res, A) { Bf <- pmax(res$Theta %*% A, 0); cs <- colSums(Bf); cs[cs == 0] <- 1; sweep(Bf, 2, cs, "/") }

# ---------------------------------------------------------------- Fig 1
cat("\n===== Figure 1: Orthodont =====\n")
d <- Orthodont; su <- unique(d$Subject); N <- length(su)
Y <- matrix(d$distance, 4, N)
ml <- as.vector(matrix(ifelse(d$Sex == "Male", 1, 0), 4, N)[1, ])
A <- rbind(intercept = rep(1, N), male = ml)
f1 <- ECMfit(Y, A, 1, "free", "seminmf", "Trend")
cat(sprintf("dfU/(NQ)=%.4f  male=%.3f\n", f1$dfU_frac, f1$Theta[1, "male"]))
XB <- f1$XB; XBb <- f1$XB.blup; age <- c(8, 10, 12, 14)
pdf("Figure1_regen_fullEM.pdf", height = 7, width = 10)
par(mar = c(4,4,1,1) + 0.1, cex = 1.2)
plot(age, XB[, 1], type = "n", ylim = range(Y), xlab = "Age", ylab = "Distance")
for (j in 1:N) {
  points(age, Y[, j], pch = ifelse(ml[j] == 1, 4, 1), col = "black", cex = 0.8)
  lines(age, XBb[, j], col = "black", lty = 3, lwd = 1.5)
  lines(age, XB[, j],  col = "black", lwd = 3)
}
legend("topleft",
       legend = c("Fixed effect", "Fixed effect + Random effect (BLUP)", "Male", "Female"),
       lwd = c(3, 1.5, NA, NA), lty = c(1, 3, NA, NA), pch = c(NA, NA, 4, 1),
       col = "black", pt.cex = 1, bg = "white", cex = 0.9)
dev.off(); cat("WROTE Figure1_regen_fullEM.pdf\n")

# ---------------------------------------------------------------- Fig 3
cat("\n===== Figure 3: OD-HUB (Q=3) =====\n")
dd <- as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx", sheet = 1, skip = 2))
dd <- dd[dd[,1]=="1" & dd[,3]=="02" & dd[,5]!="100" & dd[,7]!="100", ]
pref <- unique(dd[,6])
M <- matrix(NA, 47, 47); rownames(M) <- colnames(M) <- pref
dd[,5] <- as.numeric(dd[,5]); dd[,7] <- as.numeric(dd[,7])
for (i in 1:47) for (j in 1:47) M[i,j] <- dd[which(dd[,5]==i & dd[,7]==j), 9]
M <- log(1 + M)
hb <- c(13, 27, 40)          # Tokyo, Osaka, Fukuoka
f3 <- ECMfit(M[-hb, ], M[hb, ], 3, "nonneg", "mu", "Region")
cat(sprintf("dfU/(NQ)=%.4f\n", f3$dfU_frac))
cat("Theta (hub x region):\n"); print(round(f3$Theta, 3))
# dominant fixed-effects component per destination (argmax of Theta a_n; raw scores
# avoid the 0/0 that normalization can produce for low-flow destinations)
Bf3 <- f3$Theta %*% M[hb, ]
B.cluster <- apply(Bf3, 2, which.max)            # length 47, prefecture-code order
pref_en <- c("Hokkaido","Aomori","Iwate","Miyagi","Akita","Yamagata","Fukushima",
 "Ibaraki","Tochigi","Gunma","Saitama","Chiba","Tokyo","Kanagawa","Niigata","Toyama",
 "Ishikawa","Fukui","Yamanashi","Nagano","Gifu","Shizuoka","Aichi","Mie","Shiga","Kyoto",
 "Osaka","Hyogo","Nara","Wakayama","Tottori","Shimane","Okayama","Hiroshima","Yamaguchi",
 "Tokushima","Kagawa","Ehime","Kochi","Fukuoka","Saga","Nagasaki","Kumamoto","Oita",
 "Miyazaki","Kagoshima","Okinawa")
mycol <- c("white", "gray80", "gray45")          # Region1 / Region2 / Region3 (recolored)
pdf("Figure3_regen_fullEM.pdf", height = 7, width = 10)
par(mar = c(1,1,1,1) + 0.1, cex = 1.2)
jp <- JapanPrefMap(col = mycol[B.cluster], main = NULL, border = "gray60")
rownames(jp) <- pref_en
points(jp[hb, ], pch = 19, cex = 0.8)                      # small filled dots at hubs
text(jp[hb, 1], jp[hb, 2], c("Tokyo", "Osaka", "Fukuoka"),
     font = 2, cex = 1.1, pos = c(4, 2, 2))                # bold labels, Tokyo right / others left
legend("bottomright", title = "Dominant component",
       legend = c("Region1 (Tokyo)", "Region2 (Osaka)", "Region3 (Fukuoka)"),
       fill = mycol, bg = "white")
dev.off(); cat("WROTE Figure3_regen_fullEM.pdf\n")

# ---------------------------------------------------------------- Fig 4
cat("\n===== Figure 4: Topic (Q=3) =====\n")
corp <- corpus(data_corpus_inaugural)
tk <- tokens_remove(tokens(corp), pattern = stopwords("en", source = "snowball"))
dfm0 <- dfm_trim(dfm_select(dfm(tk), min_nchar = 3), min_termfreq = 100)
dm <- as.matrix(dfm0); dm <- dm[, order(colSums(dm), decreasing = TRUE)]
At <- rbind(pre1865  = ifelse(corp$Year <= 1865, 1, 0),
            y1865_45 = ifelse(corp$Year > 1865 & corp$Year <= 1945, 1, 0),
            after1945 = ifelse(corp$Year > 1945, 1, 0))
Yt <- t(dm); colnames(At) <- colnames(Yt)
f4 <- ECMfit(Yt, At, 3, "nonneg", "mu", "Topic")
cat(sprintf("dfU/(NQ)=%.4f\n", f4$dfU_frac))
# featured words per LEARNED topic (to determine content relabelling)
Xp <- f4$X / rowSums(f4$X)
for (q in 1:3) {
  cat(sprintf("-- learned topic %d featured words --\n", q))
  print(paste0(rownames(Xp), "(", round(100*Xp[,q],1), "%)")[order(Xp[,q], decreasing=TRUE)][1:8])
}
# content relabelling (from featured words above): learned 1=institutional,
# 2=civic, 3=union/foreign. Paper display order T1 institutional, T2 union/foreign,
# T3 civic  =>  display = learned c(1, 3, 2).
perm <- c(1, 3, 2)
Bp4 <- Bprob(f4, At)[perm, ]
lab4 <- c("Topic1 (institutional)", "Topic2 (union/foreign)", "Topic3 (civic)")
mycol4 <- c("white", "gray80", "gray90")          # T1 / T2 / T3 (stacked bottom->top)
pdf("Figure4_regen_fullEM.pdf", height = 7, width = 10)
par(mar = c(8,4,1,1) + 0.1, cex = 1.0)
barplot(Bp4, col = mycol4, legend = FALSE, las = 3, border = "gray70",
        cex.names = 0.8, ylab = "Topic proportion")
legend("topright", legend = lab4, bg = "white", fill = mycol4)
dev.off(); cat("WROTE Figure4_regen_fullEM.pdf\n")
cat("\nDONE.\n")
