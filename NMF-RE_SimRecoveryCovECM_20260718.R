# ======================================================================
# NMF-RE simulation -- 2026-07-18 re-run on the CORRECTED core.
#   NEW file; does NOT modify NMF-RE_SimRecoveryCovECM_noshift_20260715.R.
#   Changes vs the 07-15 run:
#   (a) sources nmfre_core20260718.R: the row-centering of U is removed (audit item 2 --
#       the fixed point then satisfies U A' = 0, the actual identifying condition, and
#       Thetahat becomes the exact root of the profiled estimating equation that Section 4
#       and Appendix B assume), and the initial Theta is no longer sign-clipped under
#       theta_sign='free' (audit item 17 -- this had been resetting the sign of the effect
#       on every C2 warm start).
#   (b) records the CONVERGENCE of every C2 refit (audit item 24): the 07-15 run accepted a
#       refit whenever Theta existed, without checking fb$converged, while capping it at
#       maxit=300 / outer_maxit=25.  We now count converged / non-converged refits and
#       report the rate, so the C2 column can be read as a full-refit result or not.
# ======================================================================
# ======================================================================
# NMF-RE Revision Simulation -- FULL-ECM RECOVERY + COVERAGE (C1 & C2)  (2026-07-15)
#   NEW file; does NOT modify NMF-RE_SimRecoveryCovECM_20260628.R.
#   CHANGE vs the 20260628 run: the non-negativity guard on the C2 bootstrap
#   sample (Yb <- Yb - min(Yb) + 1e-8) is now applied ONLY to the non-negative
#   arm, which needs it because its multiplicative updates require Y >= 0.
#   The SIGN-FREE arm does not need it: the Ding-Li-Jordan semi-NMF X-step
#   splits Y B' into positive/negative parts by construction and the core makes
#   no Y >= 0 assumption (verified in check_negY_seminmf_20260715.R).  The old
#   guard fired on ~2% of refits (up to 5.31%; shifts up to ~45 on data of scale
#   ~20) and, being triggered by extreme multiplier draws, selectively altered
#   the TAIL of the bootstrap distribution -- see check_C2_shift_rate_20260715.R
#   / C2_shift_rate_20260715.csv.  The sign-free C2 is therefore now a standard
#   wild bootstrap.  Seeds are unchanged, so the non-negative arm reproduces the
#   20260628 numbers bit-for-bit and only the sign-free C2 columns move.
#   n_shift is now recorded per replicate so the residual (non-negative-arm)
#   rate can be reported honestly.
#   Produces, in ONE run, BOTH:
#     - Table D1 (recovery): bias/RMSE of theta_male and composite C_male,
#       NMF-RE(free/nonneg) vs MIMIC vs OLS;
#     - Table D2 (coverage): conditional one-step C1 vs FULL-REFIT C2, now
#       with the variance components (sigma2,tau2) ESTIMATED BY ECM both for
#       the original fit AND re-estimated on every C2 bootstrap refit.
#   This makes D2 the full-ECM coverage table, so the old fixed-variance D2
#   and the separate N=27 full-ECM arm (Table D3) are SUBSUMED -> D3 dropped.
#   NMF-RE starts from sigma2=tau2=1 (truth NOT supplied) -> level footing
#   with MIMIC/OLS. Truth sigma2=1.5, tau2=50 (lambda=0.03, RE-dominant);
#   NMF-RE must RECOVER these. PSOCK-parallel (Windows-safe); C2 is the cost,
#   but ECM at Q=1,P=4 is cheap so the full run finishes in minutes.
# ======================================================================

source("nmfre_core20260718.R")

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(parallel)
  have_lavaan <- requireNamespace("lavaan", quietly = TRUE)
  have_nmfkc  <- requireNamespace("nmfkc",  quietly = TRUE)
})

# ----------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------
SMOKE <- (Sys.getenv("SIM_FULL") != "1")
TINY  <- (Sys.getenv("SIM_TINY") == "1")
if (TINY) {
  N_vec  <- c(27);      R_mc <- 4;   wild_B <- 100
  err_dist_vec_override <- c("gaussian")
} else if (SMOKE) {
  N_vec  <- c(27, 100); R_mc <- 20;  wild_B <- 200
} else {                # publication scale
  N_vec  <- c(27, 108); R_mc <- 500; wild_B <- 1000
}
full_B <- if (TINY) 6L else if (SMOKE) 20L else 40L   # C2 full-refit reps (ECM re-est.; 40 plenty for the SE-based C2 averaged over R_mc)
err_dist_vec <- c("gaussian", "exp_centered")
if (exists("err_dist_vec_override")) err_dist_vec <- err_dist_vec_override
sim_err_env <- Sys.getenv("SIM_ERR")
if (sim_err_env != "") err_dist_vec <- sim_err_env
out_tag <- if (sim_err_env != "") paste0("_", sim_err_env) else ""
scenarios <- list(A_null = 0, B_alt = 9.4285, C_neg = -9.4285)

Q <- 1; P <- 4; alpha <- 0.05
maxit <- 2000; epsilon <- 1e-6
sigma2_true <- 1.5; tau2_true <- 50    # TRUTH (data-generating); NMF-RE must recover it
X_true <- matrix(c(0.2308, 0.2409, 0.2566, 0.2717), nrow = 4, ncol = 1)
stopifnot(abs(sum(X_true) - 1) < 1e-6)
intercept_true <- 90.502
n_cores <- max(1L, parallel::detectCores() - 2L)   # use all cores (PSOCK on Windows)
seed_base_design <- 20260110; seed_base_mc <- 20260623; seed_base_full <- 20260624
maxit_boot <- 300L; outer_maxit_boot <- 25L   # cap the (warm-started) ECM refit: bootstrap theta stabilizes fast

gen_noise <- function(n, dist, sd) {
  if (sd == 0) return(rep(0, n))
  if (dist == "gaussian") rnorm(n, 0, sd) else (rexp(n, 1) - 1) * sd
}
p_male <- { d8 <- subset(nlme::Orthodont, age == 8); mean(as.integer(d8$Sex == "Male")) }
get_A <- function(N) {
  d8 <- subset(nlme::Orthodont, age == 8)
  male27 <- as.integer(d8$Sex == "Male")
  if (N %% 27 == 0) rbind(intercept = rep(1, N), male = rep(male27, N %/% 27))
  else { set.seed(seed_base_design + N); rbind(intercept = rep(1, N), male = rbinom(N, 1, p_male)) }
}

fit_mimic <- function(Y, male) {
  if (!have_lavaan) return(NULL)
  df <- as.data.frame(t(Y)); colnames(df) <- paste0("y", seq_len(nrow(Y))); df$male <- male
  mod <- paste0("f =~ ", paste(colnames(df)[1:nrow(Y)], collapse = " + "), "\n f ~ male")
  fit <- tryCatch(lavaan::sem(mod, data = df, meanstructure = TRUE, warn = FALSE, se = "standard"),
                  error = function(e) NULL)
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) return(NULL)
  est <- lavaan::parameterEstimates(fit)
  lam <- est$est[est$op == "=~"]; gam <- est$est[est$op == "~" & est$rhs == "male"]
  if (length(lam) != nrow(Y) || length(gam) != 1) return(NULL)
  list(C_male = lam * gam)
}
fit_ols_composite <- function(Y, A) t(apply(Y, 1, function(y) coef(lm(y ~ A[2, ]))[2]))

# Full-refit (C2) wild bootstrap for NMF-RE theta_male, with ECM variance
# re-estimation on EACH bootstrap refit (warm-started from the original fit).
# Resamples around the fixed-effect mean XB (= X*Theta*A), matching the C1
# residual base; resampling around XB.blup would undercover at near-saturation.
full_refit_boot <- function(Y, A, fit0, theta_sign, sim_id, N, B) {
  if (B <= 0) return(list(sd = NA_real_, lo = NA_real_, hi = NA_real_, ok = 0L, n_shift = 0L))
  fitted0 <- fit0$XB; R <- Y - fitted0; tb <- rep(NA_real_, B); n_shift <- 0L
  n_conv <- 0L; n_ran <- 0L
  xu <- if (theta_sign == "free") "seminmf" else "mu"
  for (b in 1:B) {
    set.seed(seed_base_full + sim_id + 1000 * N + b)
    w <- rexp(ncol(Y), 1) - 1
    Yb <- fitted0 + sweep(R, 2, w, `*`)
    # Guard needed ONLY by the multiplicative (nonneg/mu) updates, which require
    # Y >= 0.  The sign-free semi-NMF arm handles mixed-sign Y natively, so it now
    # gets a standard, unshifted wild bootstrap sample.
    if (theta_sign != "free") {
      ymin <- min(Yb)
      if (is.finite(ymin) && ymin < 0) { Yb <- Yb - ymin + 1e-8; n_shift <- n_shift + 1L }
    }
    fb <- tryCatch(nmf.re(Yb, A, Q = Q, X.ini = fit0$X, Theta.ini = fit0$Theta,
                          U.ini = fit0$U, theta_sign = theta_sign, x_update = xu,
                          dfU_control = "off", sigma2 = fit0$sigma2, tau2 = fit0$tau2,
                          sigma2_update = TRUE, tau2_update = TRUE,   # ECM re-estimate per refit
                          maxit = maxit_boot, outer_maxit = outer_maxit_boot,
                          epsilon = epsilon, inference = FALSE,
                          verbose = FALSE, post_message = FALSE, save_obj_trace = FALSE,
                          seed = 500000 + sim_id + b), error = function(e) NULL)
    if (!is.null(fb) && !is.null(fb$Theta)) {
      tb[b] <- fb$Theta[1, "male"]
      n_ran <- n_ran + 1L
      if (isTRUE(fb$converged)) n_conv <- n_conv + 1L
    }
  }
  ok <- sum(is.finite(tb))
  if (ok < max(5, 0.6 * B)) return(list(sd = NA_real_, lo = NA_real_, hi = NA_real_, ok = ok, n_shift = n_shift, n_conv = n_conv, n_ran = n_ran))
  list(sd = sd(tb, na.rm = TRUE), lo = quantile(tb, alpha/2, na.rm = TRUE, names = FALSE),
       hi = quantile(tb, 1 - alpha/2, na.rm = TRUE, names = FALSE), ok = ok, n_shift = n_shift, n_conv = n_conv, n_ran = n_ran)
}

# ----------------------------------------------------------------------
# One replicate (ECM variance estimation; no C2)
# ----------------------------------------------------------------------
simulate_one <- function(sim_id, N, scen, theta_male_true, err_dist, u_dist = "gaussian") {
  A <- get_A(N); male <- A[2, ]
  set.seed(seed_base_mc + sim_id + 100000 * N + 7 * match(scen, names(scenarios)))
  Theta_true <- matrix(c(intercept_true, theta_male_true), nrow = 1,
                       dimnames = list("Trend1", c("intercept", "male")))
  U <- matrix(gen_noise(N, u_dist, sqrt(tau2_true)), 1)
  E <- matrix(gen_noise(P * N, err_dist, sqrt(sigma2_true)), P, N)
  Y <- X_true %*% (Theta_true %*% A + U) + E
  ymin <- min(Y); if (is.finite(ymin) && ymin < 0) Y <- Y - ymin + 1e-8
  C_true_male <- as.vector(X_true * theta_male_true)

  ini <- tryCatch(nmfkc::nmfkc(Y = Y, A = A, Q = Q, method = "EU", epsilon = 1e-8,
                               maxit = 2000, seed = 700000 + sim_id,
                               print.trace = FALSE, print.dims = FALSE),
                  error = function(e) NULL)
  if (is.null(ini)) return(NULL)
  Xi <- ini$X; s <- colSums(Xi); s[!is.finite(s) | s <= 0] <- 1
  Xi <- sweep(Xi, 2, s, "/"); Ti <- pmax(sweep(ini$C, 1, s, "*"), 1e-6)

  # ECM: estimate sigma2,tau2 from the data (start at 1,1 -- truth NOT supplied)
  fit_arm <- function(theta_sign) {
    xu <- if (theta_sign == "free") "seminmf" else "mu"
    nmf.re(Y, A, Q = Q, X.ini = Xi, Theta.ini = Ti, prefix = "Trend",
           theta_sign = theta_sign, x_update = xu, verbose = FALSE,
           dfU_control = "off", sigma2 = 1, tau2 = 1,
           sigma2_update = TRUE, tau2_update = TRUE, maxit = maxit,
           epsilon = epsilon, inference = TRUE, wild_B = wild_B,
           wild_seed = 100000 + sim_id, seed = 500000 + sim_id,
           post_message = FALSE, save_obj_trace = TRUE)
  }
  rec <- function(tag, fit, theta_sign) {
    if (is.null(fit) || is.null(fit$coefficients)) return(NULL)
    est <- fit$Theta[1, "male"]; se <- fit$Theta.se[1, "male"]; bse <- fit$Theta.boot.sd[1, "male"]
    ci1L <- fit$Theta.ci.lower[1, "male"]; ci1U <- fit$Theta.ci.upper[1, "male"]
    C_hat_male <- as.vector(fit$X * fit$Theta[1, "male"])
    f2 <- full_refit_boot(Y, A, fit, theta_sign, sim_id, N, full_B)
    data.frame(
      method = tag, sim = sim_id, N = N, scenario = scen, err_dist = err_dist,
      theta_true = theta_male_true, est = est, se = se, bse = bse,
      Cmale_rmse = sqrt(mean((C_hat_male - C_true_male)^2)),
      Cmale_bias = mean(C_hat_male - C_true_male),
      cover_C1_se  = as.integer(theta_male_true >= est - qnorm(1-alpha/2)*se  & theta_male_true <= est + qnorm(1-alpha/2)*se),
      cover_C1_bse = as.integer(theta_male_true >= est - qnorm(1-alpha/2)*bse & theta_male_true <= est + qnorm(1-alpha/2)*bse),
      cover_C1_pct = as.integer(theta_male_true >= ci1L & theta_male_true <= ci1U),
      fbse = f2$sd, full_ok = f2$ok, boot_shift_n = f2$n_shift,
      c2_conv = f2$n_conv, c2_ran = f2$n_ran,
      cover_C2_se  = if (is.finite(f2$sd)) as.integer(theta_male_true >= est - qnorm(1-alpha/2)*f2$sd & theta_male_true <= est + qnorm(1-alpha/2)*f2$sd) else NA_integer_,
      cover_C2_pct = if (is.finite(f2$lo)) as.integer(theta_male_true >= f2$lo & theta_male_true <= f2$hi) else NA_integer_,
      sigma2_hat = fit$sigma2, tau2_hat = fit$tau2, lambda_hat = fit$sigma2 / fit$tau2,
      dfU_frac = fit$dfU_frac, stringsAsFactors = FALSE)
  }
  out <- list()
  ff <- tryCatch(fit_arm("free"),   error = function(e) NULL)
  fn <- tryCatch(fit_arm("nonneg"), error = function(e) NULL)
  out$free   <- rec("NMFRE_free",   ff, "free")
  out$nonneg <- rec("NMFRE_nonneg", fn, "nonneg")

  mm <- fit_mimic(Y, male)
  base_na <- function(tag, C_hat) data.frame(
    method = tag, sim = sim_id, N = N, scenario = scen, err_dist = err_dist,
    theta_true = theta_male_true, est = NA, se = NA, bse = NA,
    Cmale_rmse = sqrt(mean((C_hat - C_true_male)^2)), Cmale_bias = mean(C_hat - C_true_male),
    cover_C1_se = NA, cover_C1_bse = NA, cover_C1_pct = NA,
    fbse = NA, full_ok = NA, boot_shift_n = NA, c2_conv = NA, c2_ran = NA,
    cover_C2_se = NA, cover_C2_pct = NA,
    sigma2_hat = NA, tau2_hat = NA, lambda_hat = NA, dfU_frac = NA, stringsAsFactors = FALSE)
  if (!is.null(mm)) out$mimic <- base_na("MIMIC", mm$C_male)
  out$ols <- base_na("OLS", as.vector(fit_ols_composite(Y, A)))
  do.call(rbind, out[!vapply(out, is.null, logical(1))])
}

# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------
# ---- parallel backend: PSOCK on Windows (uses all cores), fork elsewhere ----
use_psock <- (.Platform$OS.type == "windows") && n_cores > 1L
cl <- NULL
if (use_psock) {
  .export_vars <- ls(globalenv(), all.names = TRUE)   # all.names=TRUE: include .dot core helpers (.seminmf_update_X_pen, .nndsvdar, ...)
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1"))
  # The core functions (nmf.re, .seminmf_*, normalize_X, ...) are already in the
  # master's globalenv (sourced at the top), so clusterExport ships them to the
  # workers -- no per-worker setwd/source needed (avoids the Japanese path failing
  # on PSOCK workers, which was returning n_ok=0).
  parallel::clusterExport(cl, varlist = .export_vars, envir = globalenv())
}
run_reps <- function(N, scen, theta, err_dist) {
  force(N); force(scen); force(theta); force(err_dist)   # avoid lazy-eval promises in the PSOCK closure
  f <- function(i) tryCatch(simulate_one(i, N, scen, theta, err_dist), error = function(e) NULL)
  if (use_psock) parallel::parLapply(cl, seq_len(R_mc), f)
  else mclapply(seq_len(R_mc), f, mc.cores = n_cores)
}
cat(sprintf("[backend] %s n_cores=%d\n", if (use_psock) "PSOCK" else "fork/serial", n_cores))

all_res <- list()
ckpt <- sprintf("sim_recovery_CovECM_20260718_raw_orth%s.csv", out_tag); first_write <- TRUE
for (err_dist in err_dist_vec) for (N in N_vec) for (scen in names(scenarios)) {
  t0 <- Sys.time()
  reps <- run_reps(N, scen, scenarios[[scen]], err_dist)
  reps <- reps[!vapply(reps, is.null, logical(1))]
  if (length(reps)) {
    d <- do.call(rbind, reps); all_res[[length(all_res) + 1]] <- d
    write.table(d, ckpt, sep = ",", row.names = FALSE, col.names = first_write, append = !first_write)
    first_write <- FALSE
  }
  cat(sprintf("[run] err=%s N=%d scen=%s n_ok=%d %.1f min\n",
              err_dist, N, scen, length(reps), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (use_psock) parallel::stopCluster(cl)
res <- do.call(rbind, all_res)

# ----------------------------------------------------------------------
# Summaries (recovery + C1 coverage + recovered variances)
# ----------------------------------------------------------------------
agg <- function(d) {
  th <- d$theta_true[1]
  data.frame(
    method = d$method[1], N = d$N[1], scenario = d$scenario[1], err_dist = d$err_dist[1], n = nrow(d),
    Bias_theta = if (all(is.na(d$est))) NA else mean(d$est) - th,
    RMSE_theta = if (all(is.na(d$est))) NA else sqrt(mean((d$est - th)^2)),
    RMSE_Cmale = mean(d$Cmale_rmse, na.rm = TRUE), Bias_Cmale = mean(d$Cmale_bias, na.rm = TRUE),
    Cover_C1_bse = mean(d$cover_C1_bse, na.rm = TRUE), Cover_C1_pct = mean(d$cover_C1_pct, na.rm = TRUE),
    Cover_C2_se = mean(d$cover_C2_se, na.rm = TRUE), Cover_C2_pct = mean(d$cover_C2_pct, na.rm = TRUE),
    sigma2_hat = mean(d$sigma2_hat, na.rm = TRUE), tau2_hat = mean(d$tau2_hat, na.rm = TRUE),
    lambda_hat = mean(d$lambda_hat, na.rm = TRUE), dfU = mean(d$dfU_frac, na.rm = TRUE),
    boot_shift_pct = 100 * sum(d$boot_shift_n, na.rm = TRUE) / (full_B * nrow(d)),
    c2_conv_pct = 100 * sum(d$c2_conv, na.rm = TRUE) / pmax(sum(d$c2_ran, na.rm = TRUE), 1))
}
key <- with(res, paste(method, N, scenario, err_dist, sep = "|"))
summ <- do.call(rbind, lapply(split(res, key), agg))
summ <- summ[order(summ$err_dist, summ$scenario, summ$N, summ$method), ]; rownames(summ) <- NULL
write.csv(summ, sprintf("sim_recovery_CovECM_20260718_summary_orth%s.csv", out_tag), row.names = FALSE)
cat("\n===== FULL-ECM RECOVERY + COVERAGE (C1 & C2) =====\n"); print(summ, digits = 3)
cat(sprintf("\nSaved: sim_recovery_CovECM_20260718_{raw,summary}_orth%s.csv\n", out_tag))
