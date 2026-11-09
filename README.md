# NMF-RE: Reproducibility Materials (revision R2)

Reproducibility materials for:

> **Identification and Inference for Covariate Effects in Non-Negative Matrix Factorization with Random Effects**
> Kenichi Satoh (Shiga University), 2026.

This is the second revision (R2, November 2026). The materials of the previous version (R1, 2026-07-22,
titled "Wild Bootstrap Inference for Covariate Effects in Non-Negative Matrix Factorization with Random
Effects") are kept unchanged in [`r1/`](r1/).

## Layout

| Path | Contents |
|---|---|
| `scripts/` | Analysis scripts (R, R Markdown, Python) |
| `logs/` | Their logs; `logs/benchmark_docx_text_20261102/` holds the text of the rendered benchmark reports |
| `results/` | Saved results (`.rds`) from which the tables and figures are drawn, including the per-replicate simulation output `out_hs_v5/`, and the LaTeX table bodies (`tab_*.tex`) that the manuscript reads |
| `download_od.R` | Fetches the e-Stat origin-destination table (not redistributed; see Data) |
| `r1/` | The R1 version |

All files in `scripts/`, `logs/` and `results/` are byte-identical copies from the author's working folder
`scripts_results/`, under the same names. In the manuscript source every computed number carries a comment
`% @src <ID> scripts_results/<file> "<text>"`; the named file is found here under `logs/` or `results/`.
Script headers refer to entries of the author's working notes (for example "work memo C81"); these notes are
not part of the repository. Some headers also name earlier files that are not included (see the last section).

## Method implementation: nmfkc versions

NMF-RE is fitted with `nmfre()` and `nmfre.inference()` of the R package `nmfkc`, development branch
(version 1.0.3). The analyses were run with libraries built from specific commits of the `develop` branch;
the logs print the library path (for example `nmfkc 1.0.3 from /home/<user>/Rlib_dev_a343c39`). The version
number alone does not identify the commit (all commits below report 1.0.3).

| nmfkc commit (develop) | Used by (log) |
|---|---|
| `a343c39abeabebe8397e485a62d276654dc17bd2` (2026-11-01) | `sim_r2_hs_v5_20261101.R` (three server logs), `refit_examples_r2_20261102.R`, `check_monotone_r2_20261102.R`, `check_hs_intercept_r2_20261102.R`, `check_odhub_signfree_r2_20261102.R`, `check_odhub_rank_r2_20261102.R`, `check_odhub_intercept_r2_20261102.R` |
| `55ba7c9f417321c0ec1c21bb5ad679fb6b14f987` (2026-11-02; three commits after a343c39) | `check_topic_signfree_final_r2_20261102.R` |
| `041bf1c` (2026-11-01) | `check_v4_truth_start_20261101.R` (one number in Appendix E) |
| not identifiable (version 1.0.2.9000, 2026-10-31) | `make_hs_truth_20261031.R` (`pkgload::load_all()` of a local working tree; its output, the simulation truth, is included as `results/sim_r2_hs_truth_20261031.rds`); `refit_hs_r2_20261031.R` (server library `~/Rlib_dev`; the log prints the version only) |
| not used | all summary, table and figure scripts, `ppr_3hub_r2_20261102.R`, `check_hs_esem_lr_20261031.R`, `check_hs_positive_manifold_20261031.R`, the benchmark reports |

The header of `sim_r2_hs_v5_20261101.R` still contains the line "Package: nmfkc 1.0.3, develop commit 041bf1c"
from the description of its predecessor (v4); its first lines and its logs give a343c39.

Install the commit used by the main analyses, for example:

```r
remotes::install_github("ksatohds/nmfkc", ref = "a343c39abeabebe8397e485a62d276654dc17bd2")
```

Scripts that take a `<lib>` argument load the package with `library(nmfkc, lib.loc = <lib>)`. To keep several
commits side by side, as was done for the analyses, build each commit and install it into its own library:

```sh
git clone https://github.com/ksatohds/nmfkc.git && cd nmfkc && git checkout a343c39
R CMD build --no-build-vignettes .
mkdir -p ~/Rlib_dev_a343c39 && R CMD INSTALL --library=$HOME/Rlib_dev_a343c39 nmfkc_1.0.3.tar.gz
```

## Data

- **Holzinger-Swineford 1939** (`lavaan`) and **Orthodont** (`nlme`) are read from the R packages.
- **Inaugural addresses** (`quanteda`): `make_topic_data_20261031.R` builds the document-term matrix
  (addresses up to 2021, N = 59, P = 67; its log shows quanteda 4.3.1); the result is included as
  `results/topic_data_20261031.rds`. The STM benchmark report does not filter by year and was rendered with
  quanteda 4.2.0; according to the header of `make_topic_data_20261031.R`, quanteda 4.3 also ships the 2025
  address.
- **Origin-destination (OD) table** (e-Stat, National Road Traffic Census, statInfId=000040170612,
  <https://www.e-stat.go.jp/stat-search/files?stat_infid=000040170612>): third-party official statistics,
  **not redistributed here**. Run `download_od.R` once; it saves `OD_estat_R3_000040170612.xlsx`, the name the
  scripts read. Scripts that read it: `refit_examples_r2_20261102.R`, `check_monotone_r2_20261102.R`,
  `check_odhub_signfree_r2_20261102.R`, `check_odhub_rank_r2_20261102.R`, `check_odhub_intercept_r2_20261102.R`,
  `plot_figure3_r2_20261102.R`, `ppr_3hub_r2_20261102.R`. No file here contains the OD matrix.
- **College** (fifth example, not reported in the manuscript): `refit_examples_r2_20261102.R` and
  `check_monotone_r2_20261102.R` also fit a College example from `nmfre_college_20261031.rds`, a file prepared
  in another project and not included. To run these two scripts without it, remove the lines that read
  `nmfre_college_20261031.rds` and the `college` entries. Each example is fitted by a separate call with its own
  seed, so by the code the other four fits do not depend on it (not checked by a rerun).

## How to run

The scripts were run in one flat folder (`scripts_results/`): they read and write files in the current working
directory. To rerun, copy `scripts/*`, `results/*` and the OD file into one folder and run there. Exceptions:

- `make_hs_truth_20261031.R`, `check_hs_positive_manifold_20261031.R` and `check_global_arcs_r2_20261102.R`
  call `setwd()` with an absolute Windows path of the author's computer, and `make_hs_truth_20261031.R` loads the
  package with `pkgload::load_all("D:/My Documents/GitHub/nmfkc")`; change these lines before running.
- `refit_hs_r2_20261031.R` loads `nmfkc` from the fixed library `~/Rlib_dev`.
- `relabel_tab_topic_theta_20261108.py` reads and writes in the folder of the script itself;
  `extract_benchmark_docx_20261102.py` reads the rendered `.docx` reports from the parent folder of the script's
  folder and writes `benchmark_docx_text_20261102/` next to the script.
- Scripts that call `parallel::mclapply()` with more than one core need Linux or macOS; on Windows pass
  `<cores>` = 1. `check_odhub_rank_r2_20261102.R` sets `mc.cores = 3` in one call regardless of `<cores>`.
- Scripts without `sink()` wrote their log by redirecting standard output (`Rscript x.R > x.log 2>&1`).
- The server runs used single-threaded BLAS: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`.

### Order

1. Data: `download_od.R`; `make_topic_data_20261031.R` (output included).
2. Simulation (Section 5, Appendix E)
   1. `make_hs_truth_20261031.R` → `sim_r2_hs_truth_20261031.rds` (included).
   2. `sim_r2_hs_v5_20261101.R <cores> <lib> out_hs_v5 <specs>` (heavy; see below) → `out_hs_v5/` (included as
      `results/out_hs_v5/`, 2800 files, so the later steps can be rerun without the simulation).
      The full design is `gaussian:<e>:1-500` and `varhet:<e>:1-200` for e = 0, 0.05, 0.1, 0.2 (2800 replicates).
      (The usage line in the header names the v4 file; the arguments are the same.)
   3. `summarize_sim_r2_hs_v5_20261102.R`, `summarize_sim_r2_hs_v5_wald_20261102.R`,
      `summarize_sim_r2_hs_v5_theta_switch_20261102.R` (read `out_hs_v5/`);
      `summarize_sim_r2_hs_v5_extra_20261102.R` (reads the first summary).
   4. `emit_tables_sim_v5_r2_20261102.R` (reads the three summaries and, for the known-anchor arm,
      `out_hs_v5/gaussian_eps_0.00/`) → `emit_tab_sim_q3_nopoint_20261108.R`.
   5. Checks: `check_v5_population_width_20261102.R`, `check_v5_single_start_rds_20261102.R`,
      `check_v5_rep158_starts_20261102.R` (also reads `out_v4_truthstart/`, the output of
      `check_v4_truth_start_20261101.R`), `check_prof_reached_20261109.R`, `check_global_arcs_r2_20261102.R`
      (after step 3.1).
3. Empirical examples (Section 6, Appendix F)
   1. `refit_examples_r2_20261102.R <cores> <lib>` → `refit_examples_r2_20261102.rds` (estimates and inference only).
   2. `compare_examples_r2_20261102.R` (component order of the topic model and OD-HUB as in R1; its inputs, earlier
      fits, are not included, its output is).
   3. `emit_tab_examples_r2_20261102.R` → `emit_tab_examples_boundary_r2_20261102.R` →
      `relabel_tab_topic_theta_20261108.py`; `emit_tab_summary_v4_r2_20261102.R`;
      `emit_tab_orthodont_r2_20261102.R`, `emit_tab_orthodont_v2_r2_20261102.R`.
   4. Figures: `make_hs_r2_tables_figure_20261102.R` (DiagrammeR), `emit_tables_examples_r2_20261102.R`,
      `plot_figure3_r2_20261102.R` (OD file, NipponMap).
   5. Monotone descent: `check_monotone_r2_20261102.R <cores> <lib>` → `plot_nll_monotone_4ex_r2_20261102.R`,
      `plot_nll_monotone_4ex_r2_20261109.R`.
   6. Checks: `check_hs_intercept_r2_20261102.R`, `check_topic_signfree_final_r2_20261102.R`,
      `check_odhub_signfree_r2_20261102.R`, `check_odhub_rank_r2_20261102.R`, `check_odhub_intercept_r2_20261102.R`
      (all `<cores> <lib>`); `check_hs_positive_manifold_20261031.R`; `refit_hs_r2_20261031.R <cores>` →
      `check_hs_esem_lr_20261031.R` (the output of the former, `refit_hs_r2_20261031.rds`, is included).
4. Benchmarks (Appendices F and G): render `LMM線形混合モデル20260112.rmd`, `MIMIC潜在変数回帰20260111.rmd`,
   `STMトピックモデル20260112.rmd` to Word (these reports date from the first submission; the text here is from
   their re-run of 2026-03-01), then `extract_benchmark_docx_20261102.py`; `ppr_3hub_r2_20261102.R` (OD file).

## Manuscript element → script → log / result

Figure files keep their working names: main Figure 1 is `Figure2_r2_20261102.pdf`, main Figure 2 is
`Figure4_r2_20261102.pdf`, Appendix Figure F1 is `Figure3_r2_20261102.pdf`.

| Element | Script(s) | Log / result |
|---|---|---|
| Section 3 (monotone descent, counts) | `check_monotone_r2_20261102.R` | `logs/check_monotone_r2_20261102.log` |
| Table 1 (Q = 3 simulation) | `sim_r2_hs_v5_20261101.R` → `summarize_sim_r2_hs_v5_wald_20261102.R`, `summarize_sim_r2_hs_v5_20261102.R`, `summarize_sim_r2_hs_v5_extra_20261102.R` → `emit_tables_sim_v5_r2_20261102.R` → `emit_tab_sim_q3_nopoint_20261108.R` | `results/tab_sim_q3_nopoint_20261108.tex`; `logs/emit_tab_sim_q3_nopoint_20261108.log` |
| Section 5 text | the summaries above, `summarize_sim_r2_hs_v5_theta_switch_20261102.R`, `check_v5_population_width_20261102.R`, `refit_examples_r2_20261102.R` | `logs/summarize_sim_r2_hs_v5*_20261102.log`, `logs/check_v5_population_width_20261102.log`, `results/tab_app_sim_xi.tex`, `logs/refit_examples_r2_20261102.log` |
| Table 2 (model sizes, variances) | `refit_examples_r2_20261102.R` → `emit_tab_summary_v4_r2_20261102.R` | `results/tab_ex_summary_v4.tex` |
| Table 3 (HS basis) | `refit_examples_r2_20261102.R` → `emit_tab_examples_r2_20261102.R` | `results/tab_hs_basis.tex` |
| Table 4 (HS covariate effects) | same | `results/tab_hs_theta.tex`; intercepts in `logs/emit_tab_examples_r2_20261102.log` |
| Figure 1 (HS path diagram) | `refit_examples_r2_20261102.R` → `make_hs_r2_tables_figure_20261102.R` | `Figure2_r2_20261102.pdf` (not included); `logs/make_hs_r2_tables_figure_20261102.log` |
| Section 6.1 text (HS) | `refit_examples_r2_20261102.R`, `check_hs_intercept_r2_20261102.R`, `emit_tab_examples_r2_20261102.R`, `check_hs_esem_lr_20261031.R`, `check_hs_positive_manifold_20261031.R` | the corresponding logs |
| Table 5 (topic model) | `make_topic_data_20261031.R`, `refit_examples_r2_20261102.R`, `compare_examples_r2_20261102.R` → `emit_tab_examples_r2_20261102.R` → `emit_tab_examples_boundary_r2_20261102.R` → `relabel_tab_topic_theta_20261108.py` | `results/tab_topic_theta_v3_20261108.tex` |
| Figure 2 (topic proportions) | `refit_examples_r2_20261102.R`, `compare_examples_r2_20261102.R` → `emit_tables_examples_r2_20261102.R` | `Figure4_r2_20261102.pdf` (not included); `logs/emit_tables_examples_r2_20261102.log` |
| Section 6.2 text (topic) | `check_topic_signfree_final_r2_20261102.R`, `emit_tables_examples_r2_20261102.R`, `refit_examples_r2_20261102.R` | the corresponding logs |
| Section 7 | `refit_examples_r2_20261102.R` | `logs/refit_examples_r2_20261102.log` |
| Appendix A (numbers) | `summarize_sim_r2_hs_v5_extra_20261102.R`, `check_global_arcs_r2_20261102.R` | the corresponding logs |
| Appendix B, Figure B1 (marginal NLL) | `check_monotone_r2_20261102.R` → `plot_nll_monotone_4ex_r2_20261109.R` (and `..._20261102.R` for one count) | `Figure_nll_monotone_4ex_r2_20261109.pdf` (not included); `logs/plot_nll_monotone_4ex_r2_2026110{2,9}.log`, `logs/check_monotone_r2_20261102.log` |
| Appendix D (profile intervals) | `check_prof_reached_20261109.R` | `logs/check_prof_reached_20261109.log` |
| Appendix E text | `sim_r2_hs_v5_20261101.R` (design constant), `make_hs_truth_20261031.R`, summaries, `check_v5_rep158_starts_20261102.R`, `check_v5_single_start_rds_20261102.R`, `refit_examples_r2_20261102.R` | the corresponding logs |
| Tables E1-E3 (fits, coverage by coefficient, composite effects) | summaries → `emit_tables_sim_v5_r2_20261102.R` | `results/tab_app_sim_fits.tex`, `results/tab_app_sim_coef.tex`, `results/tab_app_sim_xi.tex` |
| Table E4 (recovery of covariate effects) | `summarize_sim_r2_hs_v5_theta_switch_20261102.R` | `results/tab_app_sim_theta.tex` |
| Table F1 (Orthodont) | `refit_examples_r2_20261102.R` → `emit_tab_orthodont_v2_r2_20261102.R` (intercept: `emit_tab_orthodont_r2_20261102.R`) | `results/tab_orthodont_theta_v2.tex`; `logs/emit_tab_orthodont_r2_20261102.log` |
| Appendix F.1 text (LMM benchmark) | `LMM線形混合モデル20260112.rmd` | `logs/benchmark_docx_text_20261102/LMM線形混合モデル20260112.txt` |
| Appendix F.2 text (OD-HUB) | `check_odhub_signfree_r2_20261102.R`, `check_odhub_rank_r2_20261102.R`, `check_odhub_intercept_r2_20261102.R` | the corresponding logs |
| Figure F1 (OD-HUB map) | `refit_examples_r2_20261102.R` + OD file → `plot_figure3_r2_20261102.R` | `Figure3_r2_20261102.pdf` (not included); `logs/plot_figure3_r2_20261102.log` |
| Table F2 (OD-HUB effects) | `refit_examples_r2_20261102.R`, `compare_examples_r2_20261102.R` → `emit_tab_examples_r2_20261102.R` → `emit_tab_examples_boundary_r2_20261102.R` | `results/tab_odhub_theta_v2.tex` |
| Table G1 (MIMIC) | `MIMIC潜在変数回帰20260111.rmd`; NMF-RE values from `refit_examples_r2_20261102.R` | `logs/benchmark_docx_text_20261102/MIMIC潜在変数回帰20260111.txt`, `logs/refit_examples_r2_20261102.log` |
| Appendix G.2 (ESEM) | `refit_hs_r2_20261031.R` → `check_hs_esem_lr_20261031.R` | `logs/check_hs_esem_lr_20261031.log` |
| Appendix G.3 (PPR) | `ppr_3hub_r2_20261102.R` | `logs/ppr_3hub_r2_20261102.log`, `logs/ppr_3hub_r2_20261102_srv60.log` (same script on a Linux server; `_srv57` likewise) |
| Table G2 (STM) | `STMトピックモデル20260112.rmd` | `logs/benchmark_docx_text_20261102/STMトピックモデル20260112.txt` |

## Heavy steps and where they ran

- `sim_r2_hs_v5_20261101.R`: 2800 replicates on three Linux servers with single-threaded BLAS, split by
  replicate range: 90 cores, 1856 replicates, wall 7295 s (R 4.6.1, `logs/sim_r2_hs_v5_20261101_srv60.log`);
  22 cores, 472 replicates, 4885 s (R 4.5.3, `_srv57.log`); 22 cores, 472 replicates, 4092 s (R 4.5.3,
  `_srv58.log`). The median time per replicate was 37 s (`logs/summarize_sim_r2_hs_v5_20261102.log`). The servers
  had different lavaan versions (0.7.2, 0.6.20, 0.6.21; used for the ESEM-within-MIMIC comparison).
- On a Linux server with R 4.6 (from the logs; the usage lines prescribe single-threaded BLAS): `refit_examples_r2_20261102.R`
  (5 cores; the HS fit 12 s and its inference 74 s, the others shorter), `check_monotone_r2_20261102.R` (5 cores),
  `check_odhub_signfree_r2_20261102.R` (10 cores, 450 s), `check_odhub_rank_r2_20261102.R` (20 cores, about 25 s),
  `check_odhub_intercept_r2_20261102.R`, `check_hs_intercept_r2_20261102.R`, `check_topic_signfree_final_r2_20261102.R`,
  `check_v4_truth_start_20261101.R` (80 cores, 92 s), `refit_hs_r2_20261031.R`, `plot_figure3_r2_20261102.R`.
- `ppr_3hub_r2_20261102.R` was run on Windows (R 4.4.1) and on two servers (R 4.6.1, R 4.5.3). The summary, table
  and figure scripts read saved files only and are light.

## Requirements

R (the logs show 4.6.1 and 4.5.3 on the servers and 4.4.1 on Windows); packages `nmfkc` (see above), `lavaan`,
`nlme`, `quanteda`, `stm`, `readxl`, `httr`, `parallel`, `pkgload`, `DiagrammeR`, `DiagrammeRsvg`, `rsvg`,
`NipponMap`, `rmarkdown`; Python 3 for the two `.py` scripts.

## Not included

- `out_v4_truthstart/` (read by `check_v5_rep158_starts_20261102.R`) and the inputs of
  `check_v4_truth_start_20261101.R` (a task file built from the earlier v4 simulation).
- The earlier fits read by `compare_examples_r2_20261102.R`; its output is included.
- `nmfre_college_20261031.rds` (College example, not in the manuscript).
- The OD table and every file derived from it; the `.rds` files of the OD checks (not needed: the numbers are
  in the logs).
- Figures (`.pdf`, `.png`) and the rendered benchmark reports (`.docx`); they are produced by the scripts, and
  `.gitignore` excludes these file types.
- The launcher of the server runs (it contains internal host addresses and a user name); the split is described
  above.
- The manuscript files.
