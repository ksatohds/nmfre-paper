# OD-HUB with the non-negative variant (C >= 0) under the final nmfkc default (a343c39): is Q = 3 better than Q = 2 and
# Q = 4?  Author: "今はQ=3ですが，解釈性から「非負の変種」を採用すべきです．Q=2およびQ=４よりも良いことを挙げて欲しい．
# 識別性が最も良い理由になるけど，Q=4では計算できないですか？"  (NEW file, 2026-11-02; work memo C91.)
# Specifications: 3 hubs (Tokyo, Osaka, Fukuoka; P = 44) with Q = 2, 3, 4, and 4 hubs (+ Aichi; P = 43) with Q = 4
# (the R1 exploration).  For each: mix2 fits with seeds 1-3 (the best kept), NLL, degeneracy, variances, X (exact zeros,
# anchor rows = rows with share > 0.999 on one component), Theta with nmfre.inference (wild B = 1000, one-sided p),
# the dominant component of each destination (argmax of Theta a_n), and identification:
#   (a) the package region width (Q <= 3);
#   (b) for every Q, a numerical search of the observationally equivalent transformations O = 11'/Q + V R V' (R
#       orthogonal (Q-1)x(Q-1), V an orthonormal basis of the complement of 1) that keep X O >= 0 and O' Theta >= 0:
#       local (Cayley transforms of 20,000 random skew directions at step sizes 1e-4, 1e-3, 1e-2) and global (100,000
#       Haar draws of R); any feasible O other than a permutation means the region is not a point, and the range of
#       O' Theta over the feasible O found is a lower bound for the coefficient's identification interval.
# Element-wise CV (as check_college_rank_r2_20261102.R: 5 folds, 4 rounds, default start) for the 3 hubs, Q = 1..4.
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript check_odhub_rank_r2_20261102.R <cores> <lib>
# Output: check_odhub_rank_r2_20261102.{log,rds} (estimates, identification and CV summaries; no OD data)
args <- commandArgs(trailingOnly = TRUE); cores <- as.integer(args[1]); lib <- args[2]
suppressMessages({library(nmfkc, lib.loc = lib); library(readxl); library(parallel)})
nmfre <- nmfkc::nmfre; nmfre.inference <- nmfkc::nmfre.inference; nll_of <- nmfkc:::.nmfre.marginal.nll
cat(sprintf("nmfkc %s from %s | R %s.%s | cores %d | %s\n", as.character(packageVersion("nmfkc", lib.loc = lib)), lib,
            R.version$major, R.version$minor, cores, format(Sys.time())))
dd <- suppressMessages(as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx", sheet = 1, skip = 2)))
dd <- dd[dd[, 1] == "1" & dd[, 3] == "02" & dd[, 5] != "100" & dd[, 7] != "100", ]; pref <- unique(dd[, 6])
M <- matrix(NA, 47, 47); rownames(M) <- colnames(M) <- pref
dd[, 5] <- as.numeric(dd[, 5]); dd[, 7] <- as.numeric(dd[, 7])
for (i in 1:47) for (j in 1:47) M[i, j] <- dd[which(dd[, 5] == i & dd[, 7] == j), 9]
M <- log(1 + M)
pref_en <- c("Hokkaido","Aomori","Iwate","Miyagi","Akita","Yamagata","Fukushima","Ibaraki","Tochigi","Gunma","Saitama",
 "Chiba","Tokyo","Kanagawa","Niigata","Toyama","Ishikawa","Fukui","Yamanashi","Nagano","Gifu","Shizuoka","Aichi","Mie",
 "Shiga","Kyoto","Osaka","Hyogo","Nara","Wakayama","Tottori","Shimane","Okayama","Hiroshima","Yamaguchi","Tokushima",
 "Kagawa","Ehime","Kochi","Fukuoka","Saga","Nagasaki","Kumamoto","Oita","Miyazaki","Kagoshima","Okinawa")
specs <- list(h3_Q2 = list(hb = c(13, 27, 40), Q = 2), h3_Q3 = list(hb = c(13, 27, 40), Q = 3),
              h3_Q4 = list(hb = c(13, 27, 40), Q = 4), h4_Q4 = list(hb = c(13, 23, 27, 40), Q = 4))
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
options(width = 220)

## ---- identification search ----
perms <- function(n) if (n == 1) matrix(1L) else do.call(rbind, lapply(seq_len(n), function(i)
  cbind(i, perms(n - 1) + (perms(n - 1) >= i))))
ident_search <- function(X, C, seed = 1) {
  Q <- ncol(X); if (Q < 3) return(list(dim = 0, local = NA, global = NA, n_feasible = 0, range = NULL))
  set.seed(seed)
  V <- qr.Q(qr(cbind(1, diag(Q))))[, 2:Q, drop = FALSE]            # orthonormal basis of the complement of 1
  J <- matrix(1 / Q, Q, Q)
  Omat <- function(R) J + V %*% R %*% t(V)
  tolX <- 1e-10 * max(X); tolC <- 1e-10 * max(abs(C))
  feas <- function(O) all(X %*% O >= -tolX) && all(crossprod(O, C) >= -tolC)
  P <- perms(Q); isperm <- function(O) any(apply(P, 1, function(p) max(abs(O - diag(Q)[, p])) < 1e-6))
  m <- Q - 1; dimF <- m * (m - 1) / 2
  rskew <- function() { A <- matrix(0, m, m); A[upper.tri(A)] <- rnorm(dimF); A - t(A) }
  cayley <- function(A) solve(diag(m) + A, diag(m) - A)
  loc <- sapply(c(1e-4, 1e-3, 1e-2), function(h) sum(replicate(20000, {
    A <- rskew(); A <- A / sqrt(sum(A^2)) * h; feas(Omat(cayley(A))) })))
  names(loc) <- c("1e-4", "1e-3", "1e-2")
  ## global: Haar draws of R in O(m), both determinants
  feasO <- list()
  for (b in 1:100000) {
    Z <- matrix(rnorm(m * m), m, m); qz <- qr(Z); R <- qr.Q(qz) %*% diag(sign(diag(qr.R(qz))), m)
    O <- Omat(R); if (feas(O) && !isperm(O)) feasO[[length(feasO) + 1]] <- O
  }
  rng <- NULL
  if (length(feasO)) {
    Cs <- sapply(feasO, function(O) as.vector(crossprod(O, C)))
    rng <- cbind(lower = pmin(apply(Cs, 1, min), as.vector(C)), upper = pmax(apply(Cs, 1, max), as.vector(C)))
  }
  list(dim = dimF, local = loc, n_feasible = length(feasO), range = rng)
}

## ---- fits ----
res <- list()
for (nm in names(specs)) {
  sp <- specs[[nm]]; hb <- sp$hb; Q <- sp$Q
  Y <- M[-hb, ]; A <- M[hb, , drop = FALSE]
  fs <- mclapply(1:3, function(s) tryCatch(quiet(nmfre(Y, A, rank = Q, C.signed = FALSE, epsilon = 1e-7, maxit = 5e5,
                                                       seed = s, print.trace = FALSE)), error = function(e) conditionMessage(e)),
                 mc.cores = 3)
  ok <- !vapply(fs, is.character, TRUE)
  nl <- sapply(fs[ok], function(f) nll_of(Y, f$X, f$C, A, f$sigma2, f$tau2))
  b <- fs[ok][[which.min(nl)]]
  fi <- quiet(nmfre.inference(b, Y, A, wild.bootstrap = TRUE, wild.B = 1000, wild.seed = 123))
  X <- b$X; rownames(X) <- pref_en[-hb]; C <- b$C; colnames(C) <- pref_en[hb]
  sh <- X / rowSums(X); anchors <- lapply(1:Q, function(q) rownames(X)[sh[, q] > 0.999])
  dom <- apply(C %*% A, 2, which.max); names(dom) <- pref_en
  id <- ident_search(X, C)
  cat(sprintf("\n==== %s: %d hubs (%s), Q = %d, P = %d ====\n", nm, length(hb), paste(pref_en[hb], collapse = ", "), Q, nrow(Y)))
  cat(sprintf("NLL by seed %s | degenerate %s | sigma2 %.4f tau2 %.2f lambda %.5f dfU/(NQ) %.3f\n",
              paste(sprintf("%.3f", nl), collapse = " "), paste(sapply(fs[ok], function(f) isTRUE(f$degeneracy$flag)), collapse = " "),
              b$sigma2, b$tau2, b$lambda, b$dfU.frac))
  cat(sprintf("X: exact zeros %d of %d | anchor rows (share > 0.999): %s\n", sum(X < 1e-12), length(X),
              paste(sprintf("comp%d: %s", 1:Q, sapply(anchors, function(a) if (length(a)) paste(a, collapse = "/") else "none")), collapse = "; ")))
  cat("Theta (rows = components, columns = hubs):\n"); print(round(C, 3))
  print(fi$coefficients[, c("Basis", "Covariate", "Estimate", "SE", "p_value")], row.names = FALSE, digits = 4)
  cat("dominant component of each destination:\n"); print(table(dom))
  cat(paste(sprintf("%s=%d", pref_en, dom), collapse = " "), "\n")
  cat(sprintf("identification: package width %s | family dimension %s | local feasible (of 20000 at 1e-4/1e-3/1e-2): %s | global feasible non-permutation O: %d of 100000\n",
              if (is.null(b$C.ident.lower)) "NA (Q >= 4 or not computed)" else format(signif(max(b$C.ident.upper - b$C.ident.lower), 3)),
              id$dim, if (all(is.na(id$local))) "NA" else paste(id$local, collapse = "/"), id$n_feasible))
  if (!is.null(id$range)) { cat("coefficient ranges over the feasible O found (lower bound of the identification interval):\n")
    print(round(cbind(est = as.vector(C), id$range), 3)) }
  res[[nm]] <- list(hubs = pref_en[hb], Q = Q, nll = nl, X = X, C = C, sigma2 = b$sigma2, tau2 = b$tau2, dfU.frac = b$dfU.frac,
                    coefficients = fi$coefficients, anchors = anchors, dominant = dom, ident = id,
                    pkg_width = if (is.null(b$C.ident.lower)) NA else max(b$C.ident.upper - b$C.ident.lower))
}
cat(sprintf("\n2 x NLL decrease (3 hubs): Q2->Q3 %.1f, Q3->Q4 %.1f\n", 2 * (min(res$h3_Q2$nll) - min(res$h3_Q3$nll)),
            2 * (min(res$h3_Q3$nll) - min(res$h3_Q4$nll))))

## ---- element-wise CV, 3 hubs, Q = 1..4 ----
hb <- c(13, 27, 40); Y <- M[-hb, ]; A <- M[hb, , drop = FALSE]; P <- nrow(Y); N <- ncol(Y)
set.seed(1); fold <- matrix(sample(rep(1:5, length.out = P * N)), P, N); rowm <- rowMeans(Y); ridx <- row(Y)
cvt <- expand.grid(Q = 1:4, k = 1:5)
cv <- mclapply(seq_len(nrow(cvt)), function(i) {
  Q <- cvt$Q[i]; k <- cvt$k[i]; mask <- fold == k; fill <- Y; fill[mask] <- rowm[ridx[mask]]; f <- NULL; pred <- NULL
  for (it in 1:4) {
    f <- tryCatch(quiet(if (is.null(f))
           nmfre(fill, A, rank = Q, C.signed = FALSE, epsilon = 1e-5, epsilon.outer = 1e-3, inner.maxit = 1500,
                 outer.maxit = 80, maxit = 40000, seed = 1, print.trace = FALSE)
         else nmfre(fill, A, rank = Q, C.signed = FALSE, X.init = f$X, C.init = f$C, U.init = f$U, sigma2 = f$sigma2,
                    tau2 = f$tau2, epsilon = 1e-5, epsilon.outer = 1e-3, inner.maxit = 1500, outer.maxit = 80,
                    maxit = 40000, seed = 1, print.trace = FALSE)), error = function(e) NULL)
    if (is.null(f)) break
    pred <- f$XB.blup; fill[mask] <- pred[mask]
  }
  c(Q = Q, k = k, sse = if (is.null(pred)) NA else sum((Y[mask] - pred[mask])^2), n = sum(mask))
}, mc.cores = cores, mc.preschedule = FALSE)
cv <- as.data.frame(do.call(rbind, cv))
cvs <- do.call(rbind, lapply(split(cv, cv$Q), function(d) data.frame(Q = d$Q[1], rmse = sqrt(sum(d$sse) / sum(d$n)),
        fold_sd = sd(sqrt(d$sse / d$n)))))
cat("\n==== element-wise CV (3 hubs, non-negative variant, default start) ====\n"); print(cvs, row.names = FALSE, digits = 5)
cat(sprintf("row-mean baseline RMSE %.5f\n", sqrt(mean((Y - rowm)^2))))
saveRDS(list(fits = res, cv = cv, cv_summary = cvs), "check_odhub_rank_r2_20261102.rds")
cat(sprintf("\nDONE %s\n", format(Sys.time())))
