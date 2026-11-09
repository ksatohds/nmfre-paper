# Projection pursuit regression benchmark for OD-HUB with the three hubs of the paper (Tokyo, Osaka, Fukuoka), as
# Appendix G describes it (M = 3 ridge terms, up to 8 candidate terms).  NEW file, 2026-11-02 (work memo A95, A97).
# The R1 report PPR射影追跡回帰20260111.rmd fitted four hubs (Tokyo, Aichi, Osaka, Fukuoka) with nterms = 4,
# max.terms = 10 (R^2 0.500079568) although the appendix text described three hubs; that Rmd is not modified.
# Part 1 reproduces the four-hub run to check the data pipeline; Part 2 is the three-hub fit used in the paper.
# Data: the e-Stat OD table read as in refit_examples_r2_20261102.R (log(1 + F)); the RDS holds no OD data.
# Usage: Rscript ppr_3hub_r2_20261102.R      Output: ppr_3hub_r2_20261102.{log,rds}
suppressMessages(library(readxl))
dd <- suppressMessages(as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx", sheet = 1, skip = 2)))
dd <- dd[dd[, 1] == "1" & dd[, 3] == "02" & dd[, 5] != "100" & dd[, 7] != "100", ]
M <- matrix(NA, 47, 47)
dd[, 5] <- as.numeric(dd[, 5]); dd[, 7] <- as.numeric(dd[, 7])
for (i in 1:47) for (j in 1:47) M[i, j] <- dd[which(dd[, 5] == i & dd[, 7] == j), 9]
M <- log(1 + M)
pref_en <- c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki", "Tochigi", "Gunma",
  "Saitama", "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu",
  "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama",
  "Hiroshima", "Yamaguchi", "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita",
  "Miyazaki", "Kagoshima", "Okinawa")
rownames(M) <- colnames(M) <- pref_en
fit_ppr <- function(hb, nterms, max.terms) {
  X <- t(M[hb, , drop = FALSE]); Yresp <- t(M[-hb, , drop = FALSE])   # destinations x hubs, destinations x origins
  set.seed(1)
  f <- ppr(x = X, y = Yresp, nterms = nterms, max.terms = max.terms, sm.method = "supsmu")
  Yhat <- predict(f, newdata = X)
  tss <- sum((Yresp - matrix(colMeans(Yresp), nrow(Yresp), ncol(Yresp), byrow = TRUE))^2)
  list(fit = f, R2 = 1 - sum((Yresp - Yhat)^2) / tss, alpha = f$alpha, dims = c(N = nrow(X), P = ncol(Yresp), K = ncol(X)))
}
sink("ppr_3hub_r2_20261102.log", split = TRUE)
cat("ppr_3hub_r2_20261102 |", R.version.string, "|", format(Sys.time()), "\n\n")
cat("== Part 1: four hubs (Tokyo, Aichi, Osaka, Fukuoka), nterms 4, max.terms 10 (the R1 report) ==\n")
r4 <- fit_ppr(c(13, 23, 27, 40), 4, 10)
cat(sprintf("N %d P %d K %d | overall R^2 %.9f (R1 report: 0.500079568)\n", r4$dims[1], r4$dims[2], r4$dims[3], r4$R2))
cat("\n== Part 2: three hubs (Tokyo, Osaka, Fukuoka), nterms 3, max.terms 8 (the paper's OD-HUB specification) ==\n")
r3 <- fit_ppr(c(13, 27, 40), 3, 8)
cat(sprintf("N %d P %d K %d | overall R^2 %.9f | R^2 rounded %.3f\n", r3$dims[1], r3$dims[2], r3$dims[3], r3$R2, r3$R2))
cat("projection directions alpha (rows = hubs, columns = terms):\n")
a <- r3$alpha; if (is.null(dim(a))) a <- matrix(a, ncol = 1); rownames(a) <- c("Tokyo", "Osaka", "Fukuoka")
print(round(a, 3))
cat("largest |alpha| per term:", apply(abs(a), 2, function(z) rownames(a)[which.max(z)]), "\n")
cat("\nsummary of the three-hub fit:\n"); print(summary(r3$fit))
cat("DONE\n")
sink()
saveRDS(list(R2_4hub = r4$R2, R2_3hub = r3$R2, alpha_3hub = a, gofn_3hub = r3$fit$gofn,
             settings = list(three = list(hubs = c("Tokyo", "Osaka", "Fukuoka"), nterms = 3, max.terms = 8),
                             four = list(hubs = c("Tokyo", "Aichi", "Osaka", "Fukuoka"), nterms = 4, max.terms = 10))),
        "ppr_3hub_r2_20261102.rds")
