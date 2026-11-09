# HS1939: how much do the NMF-RE restrictions cost against a free ESEM-within-MIMIC?  (NEW file, 2026-10-31;
# work memo C47 follow-up.)  Author: "その案で進めて" -- on the college data (Q = 2) NMF-RE equals the
# ESEM with equal unique variances and outcome intercepts in the span of the loadings, and both
# restrictions were strongly rejected (C47).  HS is also min-max scaled with a common sigma2, so measure
# the same gap here before deciding anything for the college data.  Nothing in the HS analysis changes.
# For Q = 3 the NMF-RE covariance also needs a rotation of the loadings that is non-negative (equal
# column sums are always reachable by a rotation), so
#   free ESEM  >=  equal unique variances + intercepts in the span  >=  NMF-RE,
# and the last gap is the cost of non-negativity.
# Data as refit_hs_r2_20261031.R (min-max, N = 300); NMF-RE = its main fit (sign-free, intercept).
# ESEM: lavaan, efa block of 3 factors, factors regressed on age.rev, sex.2, school.GW, fixed.x,
# meanstructure; variants free / eqvar / span / eqvar_span; Xi = Lambda Gamma (delta SEs, unrotated).
# Usage: Rscript check_hs_esem_lr_20261031.R      (server; reads refit_hs_r2_20261031.rds)
# Output: check_hs_esem_lr_20261031.{log,rds}
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
suppressMessages(library(lavaan))
h <- readRDS("refit_hs_r2_20261031.rds"); Y <- h$Y; A0 <- h$A0; b <- h$stage1$free_int$best
nll <- min(h$stage1$free_int$table$nll)
ys <- rownames(Y); cv <- rownames(A0); P <- length(ys); Q <- 3; K0 <- length(cv)
cat(sprintf("lavaan %s | P %d, N %d, covariates %d | NMF-RE NLL %.4f | %s\n", as.character(packageVersion("lavaan")),
            P, ncol(Y), K0, nll, format(Sys.time())))
kc <- match(cv, colnames(b$C)); xi_n <- b$X %*% b$C[, kc]; dimnames(xi_n) <- list(ys, cv)
lab <- outer(seq_len(P), 1:Q, function(i, j) paste0("l", i, j))
gam <- outer(1:Q, seq_len(K0), function(j, k) paste0("g", j, k))
dat <- data.frame(t(Y), t(A0)); names(dat) <- c(ys, cv)
build <- function(eqvar, span) {
  m <- paste0(
    paste0('efa("efa1")*f', 1:Q, " =~ ", sapply(1:Q, function(j) paste0(lab[, j], "*", ys, collapse = " + ")), collapse = "\n"), "\n",
    paste0("f", 1:Q, " ~ ", sapply(1:Q, function(j) paste0(gam[j, ], "*", cv, collapse = " + ")), collapse = "\n"))
  if (eqvar) m <- paste0(m, "\n", paste0(ys, " ~~ th*", ys, collapse = "\n"))
  if (span) m <- paste0(m, "\n", paste0(ys, " ~ 0*1", collapse = "\n"), "\n", paste0("f", 1:Q, " ~ 1", collapse = "\n"))
  m
}
fit_variant <- function(eqvar, span) {
  m0 <- build(eqvar, span)
  pt <- parTable(sem(m0, data = dat, rotation = "none", fixed.x = TRUE, meanstructure = TRUE, do.fit = FALSE))
  fixed0 <- pt$label[pt$op == "=~" & pt$free == 0L & nzchar(pt$label)]
  defs <- unlist(lapply(seq_len(P), function(i) lapply(seq_len(K0), function(k) {
    terms <- paste0(lab[i, ], "*", gam[, k]); keep <- !(lab[i, ] %in% fixed0)
    sprintf("xi_%d_%d := %s", i, k, paste(terms[keep], collapse = " + ")) })))
  fe <- sem(paste(c(m0, defs), collapse = "\n"), data = dat, rotation = "none", fixed.x = TRUE, meanstructure = TRUE)
  pe <- parameterEstimates(fe); pe <- pe[pe$op == ":=", ]
  xi <- matrix(pe$est, P, K0, byrow = TRUE); dimnames(xi) <- list(ys, cv)
  th <- diag(lavInspect(fe, "est")$theta)[ys]
  list(xi = xi, ll = as.numeric(logLik(fe)), converged = lavInspect(fe, "converged"), heywood = any(th < 0), theta = th)
}
vars <- list(esem = c(FALSE, FALSE), eqvar = c(TRUE, FALSE), span = c(FALSE, TRUE), eqvar_span = c(TRUE, TRUE))
dfv <- c(esem = 0, eqvar = P - 1, span = P - Q, eqvar_span = 2 * P - 1 - Q)
res <- lapply(vars, function(v) tryCatch(fit_variant(v[1], v[2]), error = function(e) conditionMessage(e)))
options(width = 200)
cat("\nvariant      conv  heywood  logLik   2*dLL(vs free)  df   cor(Xi, NMF-RE)  max|diff|\n")
for (nm in names(res)) {
  r <- res[[nm]]; if (is.character(r)) { cat(nm, "ERROR", r, "\n"); next }
  cat(sprintf("%-12s %-5s %-7s %8.2f %14.2f %4d %16.3f %10.3f\n", nm, r$converged, r$heywood, r$ll, 2 * (res$esem$ll - r$ll),
              dfv[[nm]], cor(as.vector(xi_n), as.vector(r$xi)), max(abs(xi_n - r$xi))))
}
cat(sprintf("NMF-RE       logLik %.2f; 2*dLL vs free ESEM %.2f; vs eqvar_span (cost of non-negativity) %.2f\n",
            -nll, 2 * (res$esem$ll + nll), 2 * (res$eqvar_span$ll + nll)))
cat("\nunique variances (free ESEM):\n"); print(round(res$esem$theta, 4))
cat(sprintf("NMF-RE sigma2 %.4f\n", b$sigma2))
saveRDS(list(variants = res, xi_nmfre = xi_n, nll_nmfre = nll), "check_hs_esem_lr_20261031.rds")
cat(sprintf("\nDONE %s\n", format(Sys.time())))
