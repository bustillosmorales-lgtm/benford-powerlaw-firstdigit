# Data sources and exact queries

Retrieval date of the archives used in the paper: **2026-06-27** (unless noted).
Archive integrity is fixed by `CHECKSUMS.md`.

| Dataset | Provider | Access | Query / selection | Variable |
|---|---|---|---|---|
| Election returns | MIT Election Data and Science Lab, Harvard Dataverse | DOI 10.7910/DVN/VOQCHQ (V13) | County Presidential Election Returns 2000–2020; 2024 county vote counts | vote count |
| Financial statements | SEC EDGAR | Financial Statement Data Sets, 2026Q1 (`num.txt`) | numeric `value` field | reported value |
| Physical constants | NIST CODATA 2022 | `https://physics.nist.gov/cuu/Constants/Table/allascii.txt` | all numeric constant values | value |
| World Bank population | World Bank API | `https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL?date=2022&format=json&per_page=400` | indicator SP.POP.TOTL, date=2022, all economies | population |
| Astronomical catalogue | HYG Database (astronexus) | `hygdata_v41.csv` (Codeberg CURRENT) / `hyg_v37.csv` (GitHub) | `dist` column, `0 < dist < 1e5` (drops the 1e5 placeholder) | distance |
| Seismic catalogue | USGS ComCat | `https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=2024-01-01&endtime=2024-01-31&minmagnitude=4&limit=500` | first 500 events, M>=4, Jan 2024 | magnitude |
| Fibonacci | OEIS A000045 | generated in `benford_powerlaw.py` | first 500 terms | term |
| Pruitt / Laskowski (retracted) | Dryad | DOI 10.5061/dryad.33f0n | continuous measurement columns | measured value |

## Retracted / restricted cases (not fetched by script)

| Dataset | Reference | Note |
|---|---|---|
| Gino moral-virtue | Gino et al. 2015, Psychol. Sci. 26:983, retracted 2023 | Likert 1–7: first-digit test does not apply (see paper, Domain of applicability) |
| Handshake data | Schroeder et al. 2019, JPSP 116:743; correction 2024 | raw records behind restricted access; future work |
| Stapel | — | no data deposited; cited only as a retracted article |

The exact `first_digit()` extraction rule and the handling of zeros, negatives
and missing values are documented in `CLEANING.md`.

## Live vs archived

`recompute_table.py` loads seven rows from local archives (SEC, HYG, Dataverse,
Dryad) or regenerates them from a seed/formula (simulated baselines, Fibonacci).
Three rows — physical constants (NIST), World Bank and seismic (USGS) — are
downloaded **live** from third-party APIs at run time and are not snapshotted
here; they reproduce the exact N and estimators reported, but the endpoints are
versioned and should be cached with a checksum for a permanent archive.
