# Data dictionary

Variables actually used from each dataset (all others are ignored).

| Dataset | Field used | Type | Meaning |
|---|---|---|---|
| Election returns | county vote count | integer | votes cast for a candidate in a county |
| Financial statements | `value` (num.txt) | float | reported numeric value of a financial statement line |
| Physical constants | constant value | float | CODATA 2022 recommended value |
| World Bank population | `value` | integer | total population of an economy, 2022 |
| Astronomical catalogue | `dist` | float | stellar distance (parsecs), `0 < dist < 1e5` |
| Seismic catalogue | magnitude (`mag`) | float | event magnitude, M>=4 |
| Fibonacci | term | integer | n-th Fibonacci number, first 500 terms |
| Pruitt / Laskowski | continuous measurement | float | behavioural/performance measurements |

## Derived quantities (all scripts)

| Symbol | Code | Definition |
|---|---|---|
| `P_k(d)` | `Pk(k)` | power-law first-digit law, eq. (5) |
| `psi` | `psi(obs)` | classical structure index, chi-squared divergence from Benford |
| `k_hat` | `khat(obs)` | minimum-chi-squared exponent estimate (bound 0.05–50; 50 = saturation, reported as inf) |
| `psi_corr` | `psi(obs, Pk(k_hat))` | corrected index, deviation from the best-fitting power law |
| `N*psi_corr` | — | calibrated statistic, ~ chi-squared with 7 df under the fitted null |
| `p_boot` | `bootstrap.py` | parametric-bootstrap p-value (B=3000) |
| `N_min(k)` | `Nmin(k)` | reference sample size ceil(15.51 / psi_k) |
