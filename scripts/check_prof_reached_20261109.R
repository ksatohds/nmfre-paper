# Unrounded count of profile-likelihood intervals whose search reached the cutoff on both sides (NEW, 2026-11-09, A212).
# The summary log summarize_sim_r2_hs_v5_20261102.log prints prof_reached rounded (all 1).  Reads saved summaries only.
# Usage (in scripts_results/): Rscript check_prof_reached_20261109.R
s <- readRDS("summarize_sim_r2_hs_v5_20261102.rds")
sink("check_prof_reached_20261109.log", split = TRUE)
cat("check_prof_reached_20261109 |", format(Sys.time()), "\n")
cat("class:", class(s), "\n"); if (is.list(s)) cat("names:", paste(names(s), collapse = ", "), "\n")
find_col <- function(x, path = "s") {
  if (is.data.frame(x)) {
    hit <- grep("reach", names(x), value = TRUE)
    if (length(hit)) { cat("\n", path, ": columns", paste(hit, collapse = ", "), " nrow", nrow(x), "\n")
      for (h in hit) { v <- x[[h]]; cat("  ", h, ": class", class(v)[1], " min", format(min(v, na.rm = TRUE), digits = 10),
                                     " mean", format(mean(v, na.rm = TRUE), digits = 10), " n<1:", sum(v < 1, na.rm = TRUE),
                                     " NA:", sum(is.na(v)), "\n") }
      grp <- intersect(c("err", "eps"), names(x))
      if (length(grp) == 2) print(aggregate(x[hit], x[grp], function(v) c(n = length(v), not_reached = sum(v < 1, na.rm = TRUE), na = sum(is.na(v)))))
    }
  } else if (is.list(x)) for (nm in names(x)) find_col(x[[nm]], paste0(path, "$", nm))
}
find_col(s)
sink()
