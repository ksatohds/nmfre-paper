# Downloads the origin-destination table from the Japanese e-Stat portal
# (National Road Traffic Census, statInfId=000040170612) under e-Stat's terms
# of use, and saves it as OD_estat_R3_000040170612.xlsx for the analysis scripts.
# The table is third-party official statistics and is NOT redistributed in this
# repository; run this once before refit_empirical_20260718.R or regen_Figure2_20260718.R.
library(httr)
url <- "https://www.e-stat.go.jp/stat-search/file-download?statInfId=000040170612&fileKind=0"
GET(url, write_disk("OD_estat_R3_000040170612.xlsx", overwrite = TRUE))
message("Saved OD_estat_R3_000040170612.xlsx")
