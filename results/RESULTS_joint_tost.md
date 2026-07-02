# Joint (k, phi) fit + TOST recompute of Table 4 (2026-07-02)

Scripts: `joint_phi_tost.py` (estimator + P_{k,phi}), `recompute_table.py` (loads
the 10 datasets from ../repro2, fits, bootstraps, TOST). Seed 20260702, B=600,
delta=0.10, alpha=0.05, k>=1. Every dataset's observed first-digit vector is
loaded from its archive (checksums in ../repro2), so the table is reproducible.

## Why joint (k, phi) with k >= 1

The old phi=0 "normalization" was a no-op (dividing by an integer power of ten
leaves every mantissa, hence every first digit, unchanged). Estimating phi
jointly with k fixes the scale-phase bias, BUT with phi free over the nine cells
the two-parameter family P_{k,phi} over-fits: a degenerate small k with a phase
can mimic bounded-support data. Empirically, free phi made the seismic
catalogue (bounded, not a power law) drop from psi_corr=7.53 to 0.002. Centering
by the mean of log10 also fails (biases k for a genuine power law).

**Fix: restrict k >= 1** (the physical range already assumed in Section 2.2).
Then the joint fit (a) recovers a genuine power law with any phase
(k=2, phi=0.2 -> k_hat=2.01, phi_hat=0.20, psi_corr=0) and (b) preserves
discrimination (seismic psi_corr=2.71, uniform 0.245, Pruitt 0.130 all flagged).

## Table 4 (joint fit, 6 df, TOST at delta=0.10)

| Dataset | N | psi | k_hat | phi_hat | psi_corr | u_0.95 | verdict |
|---|---|---|---|---|---|---|---|
| Simulated power law k=2 | 2000 | 0.107 | 1.91 | 0.99 | 0.002 | 0.013 | consistent |
| Simulated uniform | 2000 | 0.391 | 1.69 | 0.66 | 0.245 | 0.279 | anomalous |
| Election returns | 20272 | 0.001 | 31.19 | 0.95 | 0.000 | 0.001 | consistent |
| Financial statements | 3236586 | 0.002 | 15.01 | 0.65 | 0.001 | 0.001 | consistent |
| Physical constants | 355 | 0.036 | 4.89 | 0.93 | 0.021 | 0.063 | consistent |
| World Bank population | 264 | 0.015 | 12.77 | 0.55 | 0.013 | 0.060 | consistent |
| Astronomical catalogue | 109400 | 0.018 | 3.10 | 0.17 | 0.002 | 0.003 | consistent |
| Fibonacci scaled | 500 | 0.000 | inf | -- | 0.000 | 0.022 | consistent |
| Pruitt (retracted) | 10269 | 0.373 | 1.38 | 0.75 | 0.130 | 0.147 | anomalous |
| Seismic catalogue | 500 | 7.466 | 1.00 | 0.60 | 2.707 | 2.980 | anomalous |

Verdict = equivalence decision (Prop. TOST): consistent iff u_0.95 < delta=0.10.

## Reading

- **Clean separation.** All 7 authentic datasets are certified *consistent*
  (u_0.95 < 0.10); the 3 genuinely irregular cases (fabricated uniform,
  retracted Pruitt, bounded-support seismic) are *anomalous*.
- **TOST vs point-null.** N*psi_corr and its p-value reject any large authentic
  sample (Financial N*psi_corr=1625, Astronomical 232, all p<0.01) because no
  real dataset is exactly a power law; the equivalence decision instead certifies
  closeness within delta, which is the correct large-N reading.
- **Financial N = 3,236,586** reproduces the original manuscript value exactly
  with the cleaned SEC loader. World Bank N=264 (was 265); HYG v40 (local) gives
  N=109400.
