# Is the local geomin representative also the global one?  (NEW file, 2026-11-02; work memo A109.)
# The package orients a Q = 3 fit by minimizing geomin over the arc of feasible rotations that contains the fit
# (nmfkc nmfre.R .nmfre.orient.geomin).  All feasible rotations give the same likelihood, so if the feasible set on
# the rotation circle had a second arc that is not a relabeling of the first, a fit could land on it and the local
# geomin member would differ from the population one.  This script scans the whole circle for (a) the HS fit of the
# paper (refit_examples_r2_20261102.rds, C81) and (b) the v5 simulation populations X(e) = X* M(e)
# (sim_r2_hs_truth_20261031.rds, as in sim_r2_hs_v5_20261101.R), counts the feasible arcs up to relabeling (the cyclic
# permutations are the rotations by 2 pi/3; the transpositions are reflections, off the circle), and compares the
# geomin minimum over the local arc with that over the whole feasible set.  Reads saved objects only; no fitting.
# Self-contained copies of the package's rotation and geomin criterion (nmfkc develop 55ba7c9).
# Usage: Rscript check_global_arcs_r2_20261102.R   Output: check_global_arcs_r2_20261102.{log,rds}
V <- qr.Q(qr(cbind(1, diag(3))))[, 2:3, drop = FALSE]
O3 <- function(phi) matrix(1 / 3, 3, 3) + V %*% matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2, 2) %*% t(V)
geomin <- function(X, eps = 0.01) {
  rs <- rowSums(X); keep <- rs > 0; S <- X[keep, , drop = FALSE] / rs[keep]
  sum(apply(S^2 + eps, 1, prod)^(1 / ncol(X)))
}
scan_circle <- function(X, n = 7200L, tol = 1e-10) {
  ph <- -pi + (seq_len(n) - 1) * 2 * pi / n                 # n points, the circle without the duplicate end
  ok <- vapply(ph, function(p) min(X %*% O3(p)) >= -tol, TRUE)
  i0 <- which.min(abs(ph)); ok[i0] <- TRUE
  shift3 <- n / 3                                           # rotation by 2 pi/3 = cyclic relabeling
  sym <- all(ok == ok[(seq_len(n) - 1 + shift3) %% n + 1])
  if (all(ok)) { n_arcs <- 1L; lab <- rep(1L, n) } else {
    st <- which(ok & !ok[c(n, seq_len(n - 1))])              # run starts (circular)
    lab <- integer(n); for (k in seq_along(st)) { j <- st[k]; while (ok[j] && lab[j] == 0L) { lab[j] <- k; j <- j %% n + 1 } }
    n_arcs <- length(st)
  }
  ## arcs up to relabeling: the cyclic shift maps arcs onto arcs
  orb <- if (all(ok)) 1L else length(unique(vapply(seq_len(n_arcs), function(k) {
    idx <- which(lab == k); min(vapply(0:2, function(s) min((idx - 1 + s * shift3) %% n + 1), 0)) }, 0)))
  g <- rep(NA_real_, n); g[ok] <- vapply(ph[ok], function(p) geomin(X %*% O3(p)), 0)
  loc <- lab == lab[i0]
  list(n_arcs = n_arcs, arcs_up_to_relabeling = orb, symmetric = sym,
       local_arc_deg = sum(loc) * 360 / n, feasible_deg = sum(ok) * 360 / n,
       geomin_local = min(g[loc]), geomin_global = min(g, na.rm = TRUE),
       phi_local = ph[loc][which.min(g[loc])], phi_global = ph[which.min(g)])
}
setwd("D:/My Documents/Data/NMF-RE混合線形モデルr2/scripts_results")
hs <- readRDS("refit_examples_r2_20261102.rds")$hs$fit$X
tr <- readRDS("sim_r2_hs_truth_20261031.rds"); Xs <- tr$X; Q <- ncol(Xs)
Mof <- function(e) (1 - e) * diag(Q) + e / (Q - 1) * (matrix(1, Q, Q) - diag(Q))
cases <- c(list(hs_fit = hs), setNames(lapply(c(0, 0.05, 0.1, 0.2), function(e) Xs %*% Mof(e)),
                                       sprintf("population_eps_%.2f", c(0, 0.05, 0.1, 0.2))))
res <- lapply(cases, scan_circle)
tab <- do.call(rbind, lapply(names(res), function(nm) { r <- res[[nm]]
  data.frame(case = nm, n_arcs = r$n_arcs, arcs_up_to_relabeling = r$arcs_up_to_relabeling, symmetric = r$symmetric,
             local_arc_deg = r$local_arc_deg, feasible_deg = r$feasible_deg, geomin_local = r$geomin_local,
             geomin_global = r$geomin_global, same_minimum = abs(r$geomin_local - r$geomin_global) < 1e-10) }))
sink("check_global_arcs_r2_20261102.log", split = TRUE)
cat("check_global_arcs_r2_20261102 |", format(Sys.time()), "\n")
print(tab, row.names = FALSE, digits = 6)
cat(sprintf("feasible set a single arc up to relabeling in %d of %d cases; local geomin = global geomin in %d of %d\n",
            sum(tab$arcs_up_to_relabeling == 1), nrow(tab), sum(tab$same_minimum), nrow(tab)))
cat("DONE\n")
sink()
saveRDS(list(table = tab, detail = res), "check_global_arcs_r2_20261102.rds")
