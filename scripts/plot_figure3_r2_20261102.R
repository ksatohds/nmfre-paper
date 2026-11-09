# Figure 3 (OD-HUB map) from the re-fit with the nmfkc default of a343c39 (refit_examples_r2_20261102.rds, C81).
# (NEW file, 2026-11-02; work memo C82.)  Drawing as plot_figure3_r2_20261031.R; the dominant fixed-effects component
# of each destination (argmax_q of Theta a_n) needs the hub rows of A, so this runs on the server where the OD data
# are; it saves only the 47 cluster labels.  Components are labelled Region q by the hub with the largest effect.
# Usage: Rscript plot_figure3_r2_20261102.R
# Output: Figure3_r2_20261102.pdf, odhub_cluster_r2_20261102.rds (cluster labels only), plot_figure3_r2_20261102.log
suppressMessages({library(readxl); library(NipponMap)})
dd <- suppressMessages(as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx", sheet = 1, skip = 2)))
dd <- dd[dd[, 1] == "1" & dd[, 3] == "02" & dd[, 5] != "100" & dd[, 7] != "100", ]; pref <- unique(dd[, 6])
M <- matrix(NA, 47, 47); rownames(M) <- colnames(M) <- pref
dd[, 5] <- as.numeric(dd[, 5]); dd[, 7] <- as.numeric(dd[, 7])
for (i in 1:47) for (j in 1:47) M[i, j] <- dd[which(dd[, 5] == i & dd[, 7] == j), 9]
M <- log(1 + M); hb <- c(13, 27, 40); A <- M[hb, , drop = FALSE]
C <- readRDS("refit_examples_r2_20261102.rds")$odhub$fit$C
p <- order(apply(C, 1, which.max)); stopifnot(identical(unname(sort(apply(C, 1, which.max))), 1:3))
C <- C[p, , drop = FALSE]                                   # row q = component of hub q (Tokyo, Osaka, Fukuoka)
print(round(C, 3))
B.cluster <- as.vector(apply(C %*% A, 2, which.max)); names(B.cluster) <- pref
saveRDS(B.cluster, "odhub_cluster_r2_20261102.rds")
pref_en <- c("Hokkaido","Aomori","Iwate","Miyagi","Akita","Yamagata","Fukushima",
 "Ibaraki","Tochigi","Gunma","Saitama","Chiba","Tokyo","Kanagawa","Niigata","Toyama",
 "Ishikawa","Fukui","Yamanashi","Nagano","Gifu","Shizuoka","Aichi","Mie","Shiga","Kyoto",
 "Osaka","Hyogo","Nara","Wakayama","Tottori","Shimane","Okayama","Hiroshima","Yamaguchi",
 "Tokushima","Kagawa","Ehime","Kochi","Fukuoka","Saga","Nagasaki","Kumamoto","Oita",
 "Miyazaki","Kagoshima","Okinawa")
mycol <- c("white", "gray80", "gray45")
pdf("Figure3_r2_20261102.pdf", height = 7, width = 10)
par(mar = c(1,1,1,1) + 0.1, cex = 1.2)
jp <- JapanPrefMap(col = mycol[B.cluster], main = NULL, border = "gray60")
rownames(jp) <- pref_en
points(jp[hb, ], pch = 19, cex = 0.8)
text(jp[hb, 1], jp[hb, 2], c("Tokyo", "Osaka", "Fukuoka"), font = 2, cex = 1.1, pos = c(4, 2, 2))
legend("bottomright", title = "Dominant component",
       legend = c("Region1 (Tokyo)", "Region2 (Osaka)", "Region3 (Fukuoka)"), fill = mycol, bg = "white")
dev.off()
cat("regions:", paste(sprintf("%s=%d", pref_en, B.cluster), collapse = " "), "\nWROTE Figure3_r2_20261102.pdf\n")
