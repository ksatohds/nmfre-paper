# ============================================================
# save_inaugural.R
# Run this ONCE with quanteda 4.1.0 to save the inaugural
# corpus data locally for reproducibility.
#
# Usage:
#   source("save_inaugural.R")
#
# This saves data/inaugural_corpus.RData containing:
#   inaugural_texts  - character vector of speech texts
#   inaugural_docvars - data.frame with Year, President, etc.
# ============================================================

cat("quanteda version:", as.character(packageVersion("quanteda")), "\n")

library(quanteda)
inaugural_texts   <- as.character(data_corpus_inaugural)
inaugural_docvars <- docvars(data_corpus_inaugural)
names(inaugural_texts) <- names(data_corpus_inaugural)

save(inaugural_texts, inaugural_docvars,
     file = "data/inaugural_corpus.RData")

cat("Saved:", length(inaugural_texts), "documents to data/inaugural_corpus.RData\n")
cat("Columns:", paste(colnames(inaugural_docvars), collapse = ", "), "\n")
