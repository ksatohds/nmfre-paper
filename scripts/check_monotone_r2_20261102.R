# Monotone descent of the marginal negative log-likelihood under the final nmfkc default (a343c39), for the rewrite
# of Algorithm 1 and Appendices A and D (A58).  With the exact X-step (no column renormalization) each inner iteration
# at fixed (sigma2, tau2) is an exact ECM iteration: the U-step is the E-step (posterior mean, and the posterior
# covariance S used by the X-step), and the X- and Theta-steps maximize the expected complete-data log-likelihood
# exactly in their blocks, so the marginal NLL cannot increase within an inner loop; the outer variance updates are
# EM M-steps evaluated at the converged inner point.  Checked here on the five examples of C81.
# (NEW file, 2026-11-02; work memo C87.)  Author: "アルゴリズム1と付録A・Dも直して、作業メモに記録して".
# For each example: (1) the kept start of the C81 fit (X.init.search) is refitted directly with print.trace, so the
#   outer boundaries are known (the "[outer k] inner=n" messages); the NLL must equal C81's;  (2) the per-iteration
#   changes of the marginal NLL (nll.trace) within inner loops and across outer steps, and of the working objective
#   (objfunc.iter) within inner loops;  (3) a run from the same start with (sigma2, tau2) fixed at the estimates, so
#   the whole run is one inner loop.  Figure: the marginal NLL over the whole run (5 panels, monochrome).
# Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript check_monotone_r2_20261102.R <cores> <lib>
# Output: check_monotone_r2_20261102.{log,rds} (traces only; no data), Figure_nll_monotone_r2_20261102.{pdf,png}
args <- commandArgs(trailingOnly = TRUE); cores <- as.integer(args[1]); lib <- args[2]
suppressMessages({library(nmfkc, lib.loc = lib); library(readxl); library(parallel); library(nlme)})
nll_of <- nmfkc:::.nmfre.marginal.nll
cat(sprintf("nmfkc %s from %s | R %s.%s | cores %d | OPENBLAS_NUM_THREADS %s | %s\n",
            as.character(packageVersion("nmfkc", lib.loc = lib)), lib, R.version$major, R.version$minor, cores,
            Sys.getenv("OPENBLAS_NUM_THREADS"), format(Sys.time())))
## ---- data (as refit_examples_r2_20261102.R) ----
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
           topic = list(Y = tp$Y, A = tp$A, Q = 3, signed = FALSE),
           orthodont = list(Y = Yo, A = rbind(intercept = rep(1, No), male = ml), Q = 1, signed = TRUE),
           odhub = list(Y = M[-hb, ], A = M[hb, , drop = FALSE], Q = 3, signed = FALSE))
c81 <- readRDS("refit_examples_r2_20261102.rds")

run_one <- function(nm) {
  e <- ex[[nm]]; s <- c81[[nm]]$fit$X.init.search; k <- which(s$chosen)
  st <- as.character(s$start[k]); sd <- s$seed[k]
  msgs <- character(0)
  t0 <- proc.time()[3]
  f <- withCallingHandlers(suppressWarnings(nmfre(e$Y, e$A, rank = e$Q, C.signed = e$signed, epsilon = 1e-7,
         maxit = 5e5, seed = sd, X.init = st, print.trace = TRUE)),
         message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  t_fit <- proc.time()[3] - t0
  ol <- grep("^\\[outer", msgs, value = TRUE)
  inner_n <- as.integer(sub("^\\[outer [0-9]+\\] inner=([0-9]+) .*$", "\\1", ol))
  ends <- cumsum(inner_n)                                # last iteration of each inner loop
  nll <- f$nll.trace; obj <- f$objfunc.iter; n <- length(nll)
  nll_fit <- nll_of(e$Y, f$X, f$C, e$A, f$sigma2, f$tau2)
  ## step i -> i + 1 crosses an outer step when i is the end of an inner loop
  cross <- seq_len(n - 1) %in% ends
  dn <- diff(nll) / abs(nll[-n]); do <- diff(obj) / abs(obj[-n])
  ## fixed (sigma2, tau2) from the same start: one inner loop to convergence
  t1 <- proc.time()[3]
  g <- suppressMessages(suppressWarnings(nmfre(e$Y, e$A, rank = e$Q, C.signed = e$signed, epsilon = 1e-7,
         maxit = 5e5, seed = sd, X.init = st, sigma2 = f$sigma2, tau2 = f$tau2, sigma2.update = FALSE,
         tau2.update = FALSE, inner.maxit = 5e5, print.trace = FALSE)))
  t_fix <- proc.time()[3] - t1
  gn <- g$nll.trace; go <- g$objfunc.iter; m <- length(gn)
  list(name = nm, start = st, seed = sd, nll_fit = nll_fit, nll_c81 = c81[[nm]]$nll, n = n, ends = ends,
       outer = length(inner_n), stop = f$stop.reason, converged = f$converged, nll = nll, obj = obj,
       max_rel_inc_nll_inner = max(c(dn[!cross], -Inf)), n_inc_nll_inner = sum(dn[!cross] > 0),
       max_rel_inc_nll_outer = max(c(dn[cross], -Inf)), n_inc_nll_outer = sum(dn[cross] > 0),
       max_rel_inc_obj_inner = max(c(do[!cross], -Inf)), n_inc_obj_inner = sum(do[!cross] > 0),
       fixed = list(nll = gn, obj = go, n = m, stop = g$stop.reason, converged = g$converged,
                    nll_end = gn[m], max_rel_inc_nll = max(diff(gn) / abs(gn[-m])), n_inc_nll = sum(diff(gn) > 0),
                    max_rel_inc_obj = max(diff(go) / abs(go[-m])), n_inc_obj = sum(diff(go) > 0)),
       t_fit = t_fit, t_fix = t_fix)
}
res <- mclapply(names(ex), function(nm) tryCatch(run_one(nm), error = function(e) list(name = nm, err = conditionMessage(e))),
                mc.cores = cores, mc.preschedule = FALSE)
names(res) <- names(ex)
for (z in res) {
  if (!is.null(z$err)) { cat(sprintf("\n== %s FAILED: %s\n", z$name, z$err)); next }
  cat(sprintf("\n== %s: start %s seed %d | NLL %.4f (C81 %.4f, diff %.1e) | %d iterations, %d outer steps, stop %s | %.0f s + %.0f s\n",
              z$name, z$start, z$seed, z$nll_fit, z$nll_c81, z$nll_fit - z$nll_c81, z$n, z$outer, z$stop, z$t_fit, z$t_fix))
  cat(sprintf("   marginal NLL, within inner loops: max relative increase %.2e (%d increases of %d steps)\n",
              z$max_rel_inc_nll_inner, z$n_inc_nll_inner, z$n - 1 - (z$outer - 1)))
  cat(sprintf("   marginal NLL, across outer steps: max relative increase %.2e (%d increases of %d)\n",
              z$max_rel_inc_nll_outer, z$n_inc_nll_outer, z$outer - 1))
  cat(sprintf("   working objective at fixed lambda, within inner loops: max relative increase %.2e (%d increases)\n",
              z$max_rel_inc_obj_inner, z$n_inc_obj_inner))
  cat(sprintf("   NLL from %.4f to %.4f\n", z$nll[1], z$nll[z$n]))
  fx <- z$fixed
  cat(sprintf("   fixed (sigma2, tau2): %d iterations, stop %s | NLL end %.4f | NLL max relative increase %.2e (%d) | objective %.2e (%d)\n",
              fx$n, fx$stop, fx$nll_end, fx$max_rel_inc_nll, fx$n_inc_nll, fx$max_rel_inc_obj, fx$n_inc_obj))
}
saveRDS(res, "check_monotone_r2_20261102.rds")
## ---- figure: marginal NLL over the whole run ----
labs <- c(hs = "(a) Holzinger--Swineford (Q=3, sign-free)", college = "(b) College (Q=2, sign-free)",
          topic = "(c) Topic model (Q=3, non-negative)", orthodont = "(d) Orthodont (Q=1, sign-free)",
          odhub = "(e) OD-HUB (Q=3, non-negative)")
draw <- function() {
  par(mfrow = c(3, 2), mar = c(4.2, 5.2, 2.6, 1.0), mgp = c(2.9, 0.8, 0), cex.lab = 1.05)
  for (nm in names(labs)) {
    z <- res[[nm]]; if (!is.null(z$err)) { plot.new(); next }
    it <- seq_len(z$n)
    plot(it, z$nll, type = "l", lwd = 2.2, col = "black", xlab = "Update (iteration)",
         ylab = "Marginal negative log-likelihood", main = labs[[nm]])
    brk <- z$ends[z$ends < z$n]
    if (length(brk)) rug(brk, side = 1, ticksize = 0.035, lwd = 0.8, col = "gray35")
    legend("topright", bty = "n", cex = 0.75, text.col = "gray25",
           legend = c(sprintf("%.1f to %.1f over %d iterations", z$nll[1], z$nll[z$n], z$n),
                      sprintf("%d outer variance steps (ticks)", length(brk))))
  }
}
pdf("Figure_nll_monotone_r2_20261102.pdf", width = 9, height = 10); draw(); dev.off()
png("Figure_nll_monotone_r2_20261102.png", width = 1350, height = 1500, res = 150); draw(); dev.off()
cat(sprintf("\nWROTE Figure_nll_monotone_r2_20261102.pdf/.png\nDONE %s\n", format(Sys.time())))
