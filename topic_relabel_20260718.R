# Topic example on the 2026-07-18 core: print the top words of each LEARNED component so the
# content relabelling (institutional/legal = Topic 1, union/constitutional/foreign = Topic 2,
# modern civic rhetoric = Topic 3) can be re-derived rather than assumed, then emit the
# coefficient table in that content order.
suppressMessages({library(nmfkc); library(quanteda); library(quanteda.textmodels)})
source("nmfre_core20260718.R")

corp <- corpus(data_corpus_inaugural)
tk <- tokens_remove(tokens(corp), pattern = stopwords("en", source = "snowball"))
df <- dfm_trim(dfm_select(dfm(tk), min_nchar = 3), min_termfreq = 100)
dm <- as.matrix(df); dm <- dm[, order(colSums(dm), decreasing = TRUE)]
Y <- t(dm)
A <- rbind(pre1865 = ifelse(corp$Year <= 1865, 1, 0),
           y1865_45 = ifelse(corp$Year > 1865 & corp$Year <= 1945, 1, 0),
           after1945 = ifelse(corp$Year > 1945, 1, 0))
colnames(A) <- colnames(Y); terms <- rownames(Y)

r0 <- nmfkc(Y = Y, A = A, Q = 3, epsilon = 1e-8, method = "EU", seed = 1)
f <- nmf.re(Y, A, Q = 3, X.ini = r0$X, Theta.ini = r0$C, prefix = "Topic",
            theta_sign = "nonneg", x_update = "mu", dfU_control = "off",
            sigma2 = 1, tau2 = 1, sigma2_update = TRUE, tau2_update = TRUE,
            inference = TRUE, wild_B = 1000, wild_seed = 123,
            verbose = FALSE, post_message = FALSE)

X <- f$X; rownames(X) <- terms
cat(sprintf("Topic fit: lambda=%.5g  dfU=%.2f  r=%.4f\n\n", f$sigma2/f$tau2, f$dfU, f$dfU_frac))
cat("Top 12 words of each LEARNED component:\n")
for (q in 1:3) {
  o <- order(X[, q], decreasing = TRUE)[1:12]
  cat(sprintf("  learned %d: %s\n", q, paste(terms[o], collapse = ", ")))
}
cat("\nDominant learned component per era (argmax of Theta):\n")
print(apply(f$Theta, 2, function(v) which.max(v)))
cat("\nTheta (rows = learned components, cols = eras):\n")
print(round(f$Theta, 3))

## content order used in the paper: Topic1 institutional/legal, Topic2 union/constitutional/
## foreign, Topic3 modern civic rhetoric.  Set `perm` after reading the word lists above.
perm <- c(1, 3, 2)   # learned -> content (as in the published table); verify against the words
cat(sprintf("\nApplying content permutation learned(%s) -> Topic 1,2,3\n",
            paste(perm, collapse = ",")))
cf <- f$coefficients
lab <- setNames(paste0("Topic ", seq_along(perm)), paste0("Topic", perm))
era <- c(pre1865 = "Early (1775--1865)", y1865_45 = "Industrial (1865--1945)",
         after1945 = "Postwar (1945--)")
cat("\nLaTeX rows in content order:\n")
for (e in names(era)) for (k in seq_along(perm)) {
  q <- perm[k]
  r <- cf[cf$Variable == e & cf$Basis == paste0("Topic", q), ]
  if (nrow(r) == 1)
    cat(sprintf("%s & Topic %d & %7.3f & %6.3f & %6.3f & %5.2f & %s \\\\\n",
                era[[e]], k, r$Estimate, r$SE, r$BSE, r$z_value,
                ifelse(r$p_value < 0.001, "$<0.001$", sprintf("%.4f", r$p_value))))
}
