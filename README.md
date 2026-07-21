# NMF-RE: Reproducibility Materials

Reproducibility materials for:

> **Wild Bootstrap Inference for Covariate Effects in Non-Negative Matrix Factorization with Random Effects**
> Kenichi Satoh (Shiga University), 2026.

## Method implementation

The estimation and inference are implemented in the R package `nmfkc`. The version
used by this paper is the development version, installed with:

```r
remotes::install_github("ksatohds/nmfkc@develop")
```

This repository additionally contains the paper's own analysis scripts. A self-contained
estimation core, `nmfre_core20260718.R`, is included so that the scripts run without the
package's internal build.

## Data

- **Orthodont** (from the `nlme` package), **Holzinger-Swineford 1939** (from the `lavaan`
  package), and the **inaugural-address corpus** (from `quanteda`; a cached copy is provided
  in `data/inaugural_corpus.RData`, regenerable with `save_inaugural.R`) are obtained directly
  from R packages.
- The **origin-destination (OD) table** is THIRD-PARTY official statistics from the Japanese
  e-Stat portal (National Road Traffic Census, statInfId=000040170612,
  <https://www.e-stat.go.jp/stat-search/files?stat_infid=000040170612>). It is **NOT
  redistributed here.** Run `download_od.R` once to fetch it under e-Stat's terms of use
  before running the OD scripts. The file is saved as `OD_estat_R3_000040170612.xlsx`, the
  name the analysis scripts expect.

## Repository structure

| File | Paper element | Description |
|---|---|---|
| `nmfre_core20260718.R` | — | Estimation core (`nmf.re()`): sign-free default, semi-NMF X-step, outer-inner ECM |
| `refit_empirical_20260718.R` | Table 2; Tables 3 & 5; Appendix Tables E1 & E2 | Empirical fits: model-size/variance summary and covariate-effect tables (Orthodont, HS1939, OD-HUB, topic model) |
| `show_HS_X_cluster_20260719.R` | Table 4 | Holzinger-Swineford estimated basis loadings (argmax grouping) |
| `topic_relabel_20260718.R` | Appendix Table E2 | Relabels learned topic columns to content order |
| `NMF-RE_SimRecoveryCovECM_20260718.R` | Table 1; Appendix Tables D1 & D2 | Simulation: recovery (bias/RMSE) and C1-vs-C2 coverage under full ECM |
| `make_tables.R` | Table 1; Appendix Tables D1 & D2 | Aggregate the simulation output into the LaTeX tables (pooled Monte Carlo RMSE) |
| `check_mono_diagnostic_20260718.R` | Section 5 (monotone-decrease diagnostic) | Fixed-variance monotonicity check over the simulation design |
| `plot_figures_134_20260718.R`, `regen_Figure2_20260718.R` | Figures 1-4 | Empirical figures (fits, path diagram) |
| `plot_monotone_20260718.R` | Figure A1 | Fixed-lambda inner-objective descent, four examples |
| `plot_nll_20260718.R` | Figure A2 | Marginal negative log-likelihood over the ECM run |
| `download_od.R` | — | Fetch the e-Stat OD table (run once; not redistributed) |
| `save_inaugural.R` | — | Cache the quanteda inaugural corpus locally |
| `data/inaugural_corpus.RData` | — | Cached inaugural corpus (public; regenerable) |

## Requirements

- R >= 4.3 (tested on 4.4.1).
- Packages: `nmfkc` (from GitHub, `@develop`), `nlme`, `lavaan`, `quanteda`,
  `quanteda.textmodels`, `readxl`, `httr`, `parallel`.
- The simulation uses Linux `mclapply` fork parallelism; on Windows it runs serially with
  identical results.
