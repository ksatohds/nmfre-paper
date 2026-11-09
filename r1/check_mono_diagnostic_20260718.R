# ----------------------------------------------------------------------------
# Re-measure the fixed-variance monotone-decrease diagnostic under the CURRENT core.
#
# WHY: main text (Simulation study) and the response letter quote a largest relative
# per-iteration increase of 3e-7 (sign-free) / 3e-6 (non-negative).  Those were measured
# by NMF-RE_SimRecovery20260627_fullEM.R, which sources nmfre_core20260627_fullEM.R --
# the core that still row-centred U.  The 2026-07-18 core drops that step, so the quoted
# numbers describe an algorithm the paper no longer specifies.  This script re-measures
# them; no existing script is modified.
#
# WHAT IT REPRODUCES: exactly the design, seeds, fixed variances, maxit and epsilon of
# NMF-RE_SimRecovery20260627_fullEM.R, and the same statistic
#     mono_rel = max(diff(obj_trace)) / |obj_trace[1]|
# reported as Mono_relmax = max over replicates within each cell.
# WHAT IT DROPS: the wild bootstrap (inference) and the C2 full refit.  obj_trace is
# produced by the point fit and both run strictly after convergence, so neither can
# change the statistic -- dropping them only makes the run cheap.
#
#   SIM_FULL=1 Rscript check_mono_diagnostic_20260718.R          # both error dists
#   SIM_FULL=1 SIM_ERR=gaussian Rscript check_mono_diagnostic_20260718.R
# ----------------------------------------------------------------------------
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressPackageStartupMessages(library(parallel))
source("nmfre_core20260718.R")

SMOKE <- (Sys.getenv("SIM_FULL") != "1")
if (SMOKE) { N_vec <- c(27, 100); R_mc <- 20 } else { N_vec <- c(27, 108); R_mc <- 500 }
err_dist_vec <- c("gaussian", "exp_centered")
sim_err_env  <- Sys.getenv("SIM_ERR")
if (sim_err_env != "") err_dist_vec <- sim_err_env
out_tag <- if (sim_err_env != "") paste0("_", sim_err_env) else ""
scenarios <- list(A_null = 0, B_alt = 9.4285, C_neg = -9.4285)

# --- identical to NMF-RE_SimRecovery20260627_fullEM.R ------------------------
Q <- 1; P <- 4
maxit <- 2000; epsilon <- 1e-6
sigma2_true <- 1.5; tau2_true <- 50
X_true <- matrix(c(0.2308, 0.2409, 0.2566, 0.2717), nrow = 4, ncol = 1)
stopifnot(abs(sum(X_true) - 1) < 1e-6)
intercept_true <- 90.502
n_cores <- if (.Platform$OS.type == "windows") 1L else max(1, detectCores() - 1)
seed_base_design <- 20260110; seed_base_mc <- 20260623

gen_noise <- function(n, dist, sd) {
  if (sd == 0) return(rep(0, n))
  if (dist == "gaussian") rnorm(n, 0, sd) else (rexp(n, 1) - 1) * sd
}
p_male <- { d8 <- subset(nlme::Orthodont, age == 8); mean(as.integer(d8$Sex == "Male")) }
get_A <- function(N) {
  d8 <- subset(nlme::Orthodont, age == 8)
  male27 <- as.integer(d8$Sex == "Male")
  if (N %% 27 == 0) {
    rbind(intercept = rep(1, N), male = rep(male27, N %/% 27))
  } else {
    set.seed(seed_base_design + N)
    rbind(intercept = rep(1, N), male = rbinom(N, 1, p_male))
  }
}

simulate_one <- function(sim_id, N, scen, theta_male_true, err_dist, u_dist = "gaussian") {
  A <- get_A(N)
  set.seed(seed_base_mc + sim_id + 100000 * N + 7 * match(scen, names(scenarios)))
  Theta_true <- matrix(c(intercept_true, theta_male_true), nrow = 1,
                       dimnames = list("Trend1", c("intercept", "male")))
  U <- matrix(gen_noise(N, u_dist, sqrt(tau2_true)), 1)
  E <- matrix(gen_noise(P * N, err_dist, sqrt(sigma2_true)), P, N)
  Y <- X_true %*% (Theta_true %*% A + U) + E
  ymin <- min(Y); if (is.finite(ymin) && ymin < 0) Y <- Y - ymin + 1e-8

  ini <- tryCatch(nmfkc::nmfkc(Y = Y, A = A, Q = Q, method = "EU", epsilon = 1e-8,
                               maxit = 2000, seed = 700000 + sim_id,
                               print.trace = FALSE, print.dims = FALSE),
                  error = function(e) NULL)
  if (is.null(ini)) return(NULL)
  Xi <- ini$X; s <- colSums(Xi); s[!is.finite(s) | s <= 0] <- 1
  Xi <- sweep(Xi, 2, s, "/"); Ti <- pmax(sweep(ini$C, 1, s, "*"), 1e-6)

  one <- function(theta_sign) {
    xu <- if (theta_sign == "free") "seminmf" else "mu"
    fit <- tryCatch(nmf.re(Y, A, Q = Q, X.ini = Xi, Theta.ini = Ti, prefix = "Trend",
                           theta_sign = theta_sign, x_update = xu, verbose = FALSE,
                           dfU_control = "off", sigma2 = sigma2_true, tau2 = tau2_true,
                           sigma2_update = FALSE, tau2_update = FALSE, maxit = maxit,
                           epsilon = epsilon, inference = FALSE, seed = 500000 + sim_id,
                           post_message = FALSE, save_obj_trace = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ob <- fit$obj_trace
    mono <- if (length(ob) > 1) max(diff(ob)) / abs(ob[1]) else 0
    data.frame(method = paste0("NMFRE_", theta_sign), sim = sim_id, N = N,
               scenario = scen, err_dist = err_dist, mono_rel = mono,
               n_iter = length(ob), stringsAsFactors = FALSE)
  }
  do.call(rbind, list(one("free"), one("nonneg")))
}

grid <- expand.grid(N = N_vec, scen = names(scenarios), err_dist = err_dist_vec,
                    stringsAsFactors = FALSE)
cat(sprintf("cells=%d  R_mc=%d  cores=%d  core=nmfre_core20260718.R\n",
            nrow(grid), R_mc, n_cores))
raw <- list()
for (g in seq_len(nrow(grid))) {
  N <- grid$N[g]; scen <- grid$scen[g]; ed <- grid$err_dist[g]
  t0 <- Sys.time()
  reps <- mclapply(seq_len(R_mc),
                   function(i) tryCatch(simulate_one(i, N, scen, scenarios[[scen]], ed),
                                        error = function(e) NULL),
                   mc.cores = n_cores)
  d <- do.call(rbind, Filter(Negate(is.null), reps))
  raw[[g]] <- d
  for (m in unique(d$method)) {
    v <- d$mono_rel[d$method == m]
    cat(sprintf("  [%2d/%2d] N=%3d %-6s %-12s %-13s  Mono_relmax=%.3e  (n=%d, %.1f min)\n",
                g, nrow(grid), N, scen, ed, m, max(v, na.rm = TRUE), length(v),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
raw <- do.call(rbind, raw)
write.csv(raw, paste0("mono_diagnostic_20260718_raw", out_tag, ".csv"), row.names = FALSE)

agg <- aggregate(mono_rel ~ method + scenario + N + err_dist, data = raw, FUN = max)
names(agg)[names(agg) == "mono_rel"] <- "Mono_relmax"
write.csv(agg, paste0("mono_diagnostic_20260718_summary", out_tag, ".csv"), row.names = FALSE)

cat("\n==== Mono_relmax by cell ====\n"); print(agg[order(agg$method, agg$scenario, agg$N), ])
cat("\n==== HEADLINE (max over all cells) ====\n")
for (m in sort(unique(raw$method)))
  cat(sprintf("  %-14s  %.4e\n", m, max(raw$mono_rel[raw$method == m], na.rm = TRUE)))
cat("\nDONE\n")
