# HS1939 example refitted for R2 (NEW file, 2026-10-31; work memo A10, A23, 改訂文案 §4-5).
# Author: "その案で進めて" -- refit the main-text Holzinger-Swineford example with an intercept row in A,
# by the develop nmfre() (nmfkc 1.0.2.9000: no row-centring of U, C8), with the R2 multistart.
# Data as make_hs_truth_20261031.R: lavaan's HolzingerSwineford1939, complete cases, nmfkc.normalize
# (column-wise min-max), Y = x1-x9 (P = 9), covariates age.rev (= -age in years), sex.2, school.GW.
# Stage 1 (3 configs = sign-free with intercept [main] / sign-free without [R1 setting] /
#   non-negative with intercept): starts default, SPA and NNDSVD with seeds 1-5 (15 per config), the
#   smallest marginal NLL kept; degeneracy; identification range on the equal-likelihood family
#   (family width; Q = 3, one-dimensional); geomin orientation (moves nothing when the family is a
#   point); C1 (sandwich SE, one-step wild bootstrap B = 1000); joint Wald by covariate; row shares of X;
#   sigma2, tau2, lambda, dfU and the saturation ratio dfU / (NQ).
# Stage 2 (main config): C2 with marginal-residual multipliers and a parametric bootstrap from the
#   fitted model, B = 480 each (12 chunks of 40), refits warm-started, SPA restart when degenerate,
#   geomin orientation, columns aligned to the estimate; also Xi = X Theta.
# Usage: Rscript refit_hs_r2_20261031.R <cores>      (server)
# Output: refit_hs_r2_20261031.{log,rds}  (the rds is written after stage 1 and again after stage 2)
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE); cores <- as.integer(args[1])
suppressMessages({library(nmfkc, lib.loc = "~/Rlib_dev"); library(parallel)})
nmfre <- nmfkc::nmfre; nmfre.inference <- nmfkc::nmfre.inference
source("nmfre_r2_helpers_20261031.R")
data(HolzingerSwineford1939, package = "lavaan")
d <- HolzingerSwineford1939; d <- d[complete.cases(d), ]
d$age.rev <- -(d$ageyr + d$agemo / 12); d$sex.2 <- ifelse(d$sex == 2, 1, 0)
d$school.GW <- ifelse(d$school == "Grant-White", 1, 0)
d[, c("id", "sex", "ageyr", "agemo", "school", "grade")] <- NULL; d <- nmfkc::nmfkc.normalize(d)
ex <- c("age.rev", "sex.2", "school.GW"); en <- setdiff(colnames(d), ex)
Y <- t(as.matrix(d[, en])); A0 <- t(as.matrix(d[, ex])); Q <- 3
cat(sprintf("nmfkc %s | R %s.%s | P = %d tests, N = %d pupils, covariates %d | %s\n",
            as.character(packageVersion("nmfkc", lib.loc = "~/Rlib_dev")), R.version$major, R.version$minor,
            nrow(Y), ncol(Y), nrow(A0), format(Sys.time())))
truth <- readRDS("sim_r2_hs_truth_20261031.rds")
cat(sprintf("same Y as the simulation truth: %s\n", isTRUE(all.equal(unname(Y), unname(truth$Y_obs)))))
cfgs <- list(free_int = list(signed = TRUE, int = TRUE), free_noint = list(signed = TRUE, int = FALSE),
             nonneg_int = list(signed = FALSE, int = TRUE))
Aof <- function(cf) if (cf$int) rbind(intercept = 1, A0) else A0

## ---- stage 1: 3 configs x 15 starts, one fit per task ----
tasks1 <- expand.grid(cfg = names(cfgs), start = c("default", "spa", "nndsvd"), seed = 1:5, stringsAsFactors = FALSE)
fits1 <- mclapply(seq_len(nrow(tasks1)), function(i) {
  tk <- tasks1[i, ]; cf <- cfgs[[tk$cfg]]; A <- Aof(cf); t0 <- proc.time()[3]
  f <- tryCatch(suppressMessages(suppressWarnings(nmfre(Y, A, rank = Q, C.signed = cf$signed, epsilon = 1e-7,
         maxit = 5e5, seed = tk$seed, X.init = if (tk$start == "default") NULL else tk$start,
         print.trace = FALSE))), error = function(e) conditionMessage(e))
  list(fit = f, time = proc.time()[3] - t0)
}, mc.cores = cores, mc.preschedule = FALSE)
options(width = 200)
st1 <- list()
for (nm in names(cfgs)) {
  cf <- cfgs[[nm]]; A <- Aof(cf); ii <- which(tasks1$cfg == nm); tab <- NULL; fits <- list()
  for (i in ii) {
    f <- fits1[[i]]$fit
    if (is.character(f)) { cat(sprintf("%s %s seed %d: ERROR %s\n", nm, tasks1$start[i], tasks1$seed[i], f)); next }
    fits[[length(fits) + 1]] <- f
    tab <- rbind(tab, data.frame(start = tasks1$start[i], seed = tasks1$seed[i], nll = r2_nll(Y, f, A),
                                 degenerate = isTRUE(f$degeneracy$flag), cos = f$degeneracy$max.cosine,
                                 canc = f$degeneracy$cancellation, stop = f$stop.reason, outer = f$outer.iter,
                                 time = fits1[[i]]$time))
  }
  best <- fits[[which.min(tab$nll)]]; best$C.signed <- cf$signed
  ir <- if (cf$signed) r2_ident_range(best) else NULL
  bo <- if (cf$signed) r2_orient(best, "geomin")$fit else best
  moved <- max(abs(bo$C - best$C))
  fi <- suppressMessages(suppressWarnings(nmfre.inference(bo, Y, A, wild.bootstrap = TRUE, wild.B = 1000,
                                                          wild.seed = 123)))
  st1[[nm]] <- list(name = nm, cfg = cf, table = tab, best = bo, ident = ir, moved = moved, inference = fi)
  cat(sprintf("\n==== %s ====\n", nm)); print(tab, digits = 7, row.names = FALSE)
  N <- ncol(Y)
  cat(sprintf("best NLL %.4f | starts within 0.01 of it: %d of %d | degenerate %s (cos %.3f, canc %.2f)\n",
              min(tab$nll), sum(tab$nll <= min(tab$nll) + 0.01), nrow(tab), isTRUE(bo$degeneracy$flag),
              bo$degeneracy$max.cosine, bo$degeneracy$cancellation))
  cat(sprintf("sigma2 %.5f tau2 %.5f lambda %.5f dfU %.2f saturation dfU/(NQ) %.3f | converged %s\n",
              bo$sigma2, bo$tau2, bo$lambda, bo$dfU, bo$dfU / (N * Q), isTRUE(bo$converged)))
  if (!is.null(ir)) cat(sprintf("family: arc %.4f deg, family width max %.3g | geomin moved Theta by %.3g\n",
                                ir$arc * 180 / pi, max(ir$range["max", ] - ir$range["min", ]), moved))
  cat("basis-profile X and row shares:\n")
  print(round(cbind(bo$X, share = bo$X / rowSums(bo$X)), 3))
  print(fi$coefficients[, c("Basis", "Covariate", "Estimate", "SE", "BSE", "z_value", "p_value")], digits = 4,
        row.names = FALSE)
  cat("joint Wald by covariate:\n"); print(r2_wald_cols(fi), digits = 4, row.names = FALSE)
}
lr <- function(a, b) 2 * (min(st1[[a]]$table$nll) - min(st1[[b]]$table$nll))
cat(sprintf("\nLR intercept (sign-free): 2*dNLL = %.2f on %d df | sign-free vs non-negative (intercept): 2*dNLL = %.2f\n",
            lr("free_noint", "free_int"), Q, lr("nonneg_int", "free_int")))
cat(sprintf("best fit equals the simulation truth: max|X - X*| %.2g, max|Theta - Theta*| %.2g (after alignment)\n",
            {b <- st1$free_int$best; p <- r2_align(b$X, truth$X); max(abs(r2_permute(b, p)$X - truth$X))},
            {b <- st1$free_int$best; p <- r2_align(b$X, truth$X); max(abs(r2_permute(b, p)$C - truth$C))}))
saveRDS(list(stage1 = lapply(st1, function(z) z[c("name", "cfg", "table", "best", "ident", "moved")]),
             inference = lapply(st1, function(z) z$inference$coefficients),
             C.vec.cov = st1$free_int$inference$C.vec.cov, Y = Y, A0 = A0), "refit_hs_r2_20261031.rds")
cat(sprintf("stage 1 saved %s\n", format(Sys.time())))

## ---- stage 2: C2 (marginal) and parametric bootstrap, main config ----
z <- st1$free_int; A <- Aof(z$cfg); Bc <- 40
tasks2 <- expand.grid(scheme = c("marginal", "parametric"), chunk = 1:12, stringsAsFactors = FALSE)
bt <- mclapply(seq_len(nrow(tasks2)), function(i) {
  tk <- tasks2[i, ]
  r2_full_refit_boot(z$best, Y, A, B = Bc, rule = "geomin", Xref = z$best$X,
                     seed = 5000000 + 100000 * match(tk$scheme, c("marginal", "parametric")) + (tk$chunk - 1) * Bc,
                     epsilon = 1e-7, maxit = 5e5, scheme = tk$scheme, restart.spa = TRUE)
}, mc.cores = cores, mc.preschedule = FALSE)
boot <- list()
for (sc in c("marginal", "parametric")) {
  ii <- which(tasks2$scheme == sc)
  Cb <- do.call(rbind, lapply(bt[ii], `[[`, "C")); Xib <- do.call(rbind, lapply(bt[ii], `[[`, "Xi"))
  dg <- unlist(lapply(bt[ii], `[[`, "degenerate")); rs <- unlist(lapply(bt[ii], `[[`, "restarted"))
  est <- as.vector(z$best$C); q <- apply(Cb, 2, quantile, c(.025, .975), na.rm = TRUE)
  xi <- as.vector(z$best$X %*% z$best$C); qx <- apply(Xib, 2, quantile, c(.025, .975), na.rm = TRUE)
  cat(sprintf("\n==== free_int, %s bootstrap: %d of %d refits ok, flagged degenerate %.3f, restarted %.3f ====\n", sc,
              sum(!is.na(Cb[, 1])), nrow(Cb), mean(dg, na.rm = TRUE), mean(rs)))
  tab <- data.frame(z$inference$coefficients[, c("Basis", "Covariate", "Estimate", "SE")],
                    boot_sd = apply(Cb, 2, sd, na.rm = TRUE), boot_mad = apply(Cb, 2, mad, na.rm = TRUE),
                    low = q[1, ], high = q[2, ], covers_est = q[1, ] <= est & est <= q[2, ])
  tab$mad_over_se <- tab$boot_mad / tab$SE
  print(tab, digits = 4, row.names = FALSE)
  cat(sprintf("covariates only: median MAD/SE %.2f (range %.2f-%.2f) | intercepts: %s\n",
              median(tab$mad_over_se[tab$Covariate != "intercept"]),
              min(tab$mad_over_se[tab$Covariate != "intercept"]), max(tab$mad_over_se[tab$Covariate != "intercept"]),
              paste(sprintf("%.2f", tab$mad_over_se[tab$Covariate == "intercept"]), collapse = " ")))
  cat(sprintf("Xi: bootstrap interval covers the estimate for %d of %d entries\n",
              sum(qx[1, ] <= xi & xi <= qx[2, ]), length(xi)))
  boot[[sc]] <- list(C = Cb, Xi = Xib, X = do.call(rbind, lapply(bt[ii], `[[`, "X")), degenerate = dg,
                     restarted = rs, table = tab)
}
saveRDS(list(stage1 = lapply(st1, function(z) z[c("name", "cfg", "table", "best", "ident", "moved")]),
             inference = lapply(st1, function(z) z$inference$coefficients),
             C.vec.cov = st1$free_int$inference$C.vec.cov, Y = Y, A0 = A0, boot = boot), "refit_hs_r2_20261031.rds")
cat(sprintf("\nDONE %s\n", format(Sys.time())))
