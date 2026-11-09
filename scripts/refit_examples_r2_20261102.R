# Re-fit of the five analysis examples with the nmfkc default of a343c39 (the same code as 5255584): X.init = "mix2"
# (k-means and SPA starts, k-means restarts while degenerate), the variances started from the residual covariance
# (probabilistic PCA), the exact X-step and, for C >= 0, the exact non-negative C-step.  Earlier tables were computed
# with sigma2 = tau2 = 1 starts (C40 college, C46 HS, C53 OD-HUB, C55 Orthodont/topic, C69 topic/OD-HUB exact).
# (NEW file, 2026-11-02; work memo C81.)  Author: "サーバーで計算し直して、作業メモに記録して".
# Specifications as in check_variance_start_ppca_examples_20261101.R (R2): hs, college, orthodont sign-free with an
# intercept; topic (era dummies summing to one) and odhub (no intercept) with C >= 0.  epsilon 1e-7, maxit 5e5, seed 1
# (as in the HS simulation v5); geomin orientation for the sign-free fits.  Inference: nmfre.inference with the
# one-step wild bootstrap (B = 1000, seed 123) and, for the sign-free fits, profile-likelihood intervals of the
# covariate effects (the intercept column is not profiled; A24).
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript refit_examples_r2_20261102.R <cores> <lib>
# Output: refit_examples_r2_20261102.{log,rds} (estimates, inference and diagnostics; no Y, no U, no fitted values)
args <- commandArgs(trailingOnly = TRUE); cores <- as.integer(args[1]); lib <- args[2]
suppressMessages({library(nmfkc, lib.loc = lib); library(readxl); library(parallel); library(nlme)})
nll_of <- nmfkc:::.nmfre.marginal.nll
cat(sprintf("nmfkc %s from %s | R %s.%s | cores %d | OPENBLAS_NUM_THREADS %s | %s\n",
            as.character(packageVersion("nmfkc", lib.loc = lib)), lib, R.version$major, R.version$minor, cores,
            Sys.getenv("OPENBLAS_NUM_THREADS"), format(Sys.time())))
## ---- data (copied from check_variance_start_ppca_examples_20261101.R) ----
data(HolzingerSwineford1939, package = "lavaan")
d <- HolzingerSwineford1939; d <- d[complete.cases(d), ]
d$age.rev <- -(d$ageyr + d$agemo / 12); d$sex.2 <- ifelse(d$sex == 2, 1, 0)
d$school.GW <- ifelse(d$school == "Grant-White", 1, 0)
d[, c("id", "sex", "ageyr", "agemo", "school", "grade")] <- NULL; d <- nmfkc.normalize(d)
exh <- c("age.rev", "sex.2", "school.GW"); enh <- setdiff(colnames(d), exh)
z0 <- readRDS("nmfre_college_20261031.rds")
od <- Orthodont; su <- unique(od$Subject); No <- length(su)
Yo <- matrix(od$distance, 4, No); ml <- as.vector(matrix(ifelse(od$Sex == "Male", 1, 0), 4, No)[1, ])
tp <- readRDS("topic_data_20261031.rds")
dd <- suppressMessages(as.data.frame(read_xlsx("OD_estat_R3_000040170612.xlsx", sheet = 1, skip = 2)))
dd <- dd[dd[, 1] == "1" & dd[, 3] == "02" & dd[, 5] != "100" & dd[, 7] != "100", ]; pref <- unique(dd[, 6])
M <- matrix(NA, 47, 47); rownames(M) <- colnames(M) <- pref
dd[, 5] <- as.numeric(dd[, 5]); dd[, 7] <- as.numeric(dd[, 7])
for (i in 1:47) for (j in 1:47) M[i, j] <- dd[which(dd[, 5] == i & dd[, 7] == j), 9]
M <- log(1 + M); hb <- c(13, 27, 40)
ex <- list(hs = list(Y = t(as.matrix(d[, enh])), A = rbind(intercept = 1, t(as.matrix(d[, exh]))), Q = 3, signed = TRUE),
           college = list(Y = z0$Y1, A = rbind(intercept = 1, z0$Y2), Q = 2, signed = TRUE),
           orthodont = list(Y = Yo, A = rbind(intercept = rep(1, No), male = ml), Q = 1, signed = TRUE),
           topic = list(Y = tp$Y, A = tp$A, Q = 3, signed = FALSE),
           odhub = list(Y = M[-hb, ], A = M[hb, , drop = FALSE], Q = 3, signed = FALSE))
## best NLL so far (C71, nine starts, exact X-step) for the check
best <- c(hs = -1199.3415, college = -3614.3527, orthodont = 215.5869, topic = 9882.4675, odhub = 4266.3191)
keep <- c("X", "C", "sigma2", "tau2", "lambda", "dfU", "dfU.frac", "converged", "iter", "X.update", "C.nonneg.update",
          "X.init.search", "variance.start", "degeneracy", "identification", "C.ident.lower", "C.ident.upper")
keep_inf <- c("coefficients", "C.se", "C.se.boot", "C.ci.lower", "C.ci.upper", "C.vec.cov", "C.p.side",
              "C.profile.lower", "C.profile.upper", "C.profile.reached")

run_one <- function(nm) {
  e <- ex[[nm]]; msgs <- character(0)
  hold <- function(expr) withCallingHandlers(expr,
    message = function(m) { msgs <<- c(msgs, paste("message:", conditionMessage(m))); invokeRestart("muffleMessage") },
    warning = function(w) { msgs <<- c(msgs, paste("warning:", conditionMessage(w))); invokeRestart("muffleWarning") })
  t0 <- proc.time()[3]
  f <- tryCatch(hold(nmfre(e$Y, e$A, rank = e$Q, C.signed = e$signed, epsilon = 1e-7, maxit = 5e5, seed = 1,
                           orient = if (e$signed) "geomin" else "none", print.trace = FALSE)),
                error = function(err) conditionMessage(err))
  t_fit <- proc.time()[3] - t0
  if (is.character(f)) return(list(name = nm, ok = FALSE, msg = f, msgs = msgs, t_fit = t_fit))
  nll <- nll_of(e$Y, f$X, f$C, e$A, f$sigma2, f$tau2)
  K <- nrow(e$A); Q <- e$Q
  pw <- if (e$signed) (Q + 1):(Q * K) else NULL          # covariate effects: columns 2..K of C (column 1 = intercept)
  t1 <- proc.time()[3]
  fi <- tryCatch(hold(nmfre.inference(f, e$Y, e$A, wild.bootstrap = TRUE, wild.B = 1000, wild.seed = 123,
                                      profile = e$signed, profile.which = pw)),
                 error = function(err) conditionMessage(err))
  t_inf <- proc.time()[3] - t1
  list(name = nm, ok = TRUE, signed = e$signed, Q = Q, N = ncol(e$Y), P = nrow(e$Y), K = K,
       nll = nll, d2nll_vs_best = 2 * (nll - best[[nm]]), fit = f[intersect(names(f), keep)],
       inf = if (is.character(fi)) fi else fi[intersect(names(fi), keep_inf)],
       msgs = msgs, t_fit = t_fit, t_inf = t_inf)
}
## HS first (the longest: fit about 25 s, profiles of 9 effects)
res <- mclapply(c("hs", "college", "orthodont", "topic", "odhub"), run_one, mc.cores = cores, mc.preschedule = FALSE)
names(res) <- vapply(res, function(z) z$name, "")

options(width = 220)
for (z in res) {
  cat(sprintf("\n==== %s ====\n", z$name))
  if (!isTRUE(z$ok)) { cat("FAILED:", z$msg, "\n"); next }
  f <- z$fit; s <- f$X.init.search
  cat(sprintf("N %d P %d Q %d K %d | %s | NLL %.4f (2dNLL vs best of C71 %+.2e) | converged %s iter %s | fit %.1f s, inference %.1f s\n",
              z$N, z$P, z$Q, z$K, if (z$signed) "sign-free" else "C >= 0", z$nll, z$d2nll_vs_best,
              format(f$converged), format(f$iter), z$t_fit, z$t_inf))
  cat(sprintf("sigma2 %.6g tau2 %.6g lambda %.6g dfU %.4f dfU/(NQ) %.4f | variance start: %s sigma2 %.4g tau2 %.4g | degenerate %s\n",
              f$sigma2, f$tau2, f$lambda, f$dfU, f$dfU.frac, f$variance.start$rule, f$variance.start$sigma2,
              f$variance.start$tau2, format(isTRUE(f$degeneracy$flag))))
  if (!is.null(s)) { cat("starts (mix2):\n"); print(s, row.names = FALSE, digits = 8) }
  if (!is.null(f$C.ident.lower)) cat(sprintf("identification region: max width %.3g\n", max(f$C.ident.upper - f$C.ident.lower)))
  cat("X (rows = variables):\n"); print(round(f$X, 3))
  if (is.character(z$inf)) { cat("inference FAILED:", z$inf, "\n") } else {
    cf <- z$inf$coefficients
    if (!is.null(z$inf$C.profile.lower)) {
      cf$prof_lo <- as.vector(z$inf$C.profile.lower); cf$prof_hi <- as.vector(z$inf$C.profile.upper)
      cf$prof_reached <- as.vector(z$inf$C.profile.reached)
    }
    print(cf, row.names = FALSE, digits = 5)
  }
  if (length(z$msgs)) { cat("messages/warnings:\n"); cat(paste0("  ", unique(z$msgs)), sep = "\n") }
}
saveRDS(res, "refit_examples_r2_20261102.rds")
cat(sprintf("\nDONE %s\n", format(Sys.time())))
