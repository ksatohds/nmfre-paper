# Topic-model data for the R2 refits (NEW file, 2026-10-31).  Same preprocessing as R1
# (plot_figures_134_20260718.R / check_topic_purity_20261031.R): quanteda data_corpus_inaugural up to
# 2021 (N = 59; quanteda 4.3 also ships 2025), snowball stopwords removed, words of >= 3 characters
# with total frequency >= 100 (P = 67), columns ordered by frequency; A = three era dummies (no
# intercept).  Run where quanteda is installed (.58); the matrices are public data.
# Usage: Rscript make_topic_data_20261031.R      Output: topic_data_20261031.rds
suppressMessages(library(quanteda))
corp <- corpus_subset(data_corpus_inaugural, Year <= 2021)
tk <- tokens_remove(tokens(corp), pattern = stopwords("en", source = "snowball"))
df <- dfm_trim(dfm_select(dfm(tk), min_nchar = 3), min_termfreq = 100)
dm <- as.matrix(df); dm <- dm[, order(colSums(dm), decreasing = TRUE)]
A <- rbind(pre1865 = ifelse(corp$Year <= 1865, 1, 0), y1865_45 = ifelse(corp$Year > 1865 & corp$Year <= 1945, 1, 0),
           after1945 = ifelse(corp$Year > 1945, 1, 0))
Y <- t(dm); colnames(A) <- colnames(Y)
saveRDS(list(Y = Y, A = A, year = corp$Year, quanteda = as.character(packageVersion("quanteda"))), "topic_data_20261031.rds")
cat(sprintf("quanteda %s | P %d, N %d (%d-%d)\n", packageVersion("quanteda"), nrow(Y), ncol(Y), min(corp$Year), max(corp$Year)))
