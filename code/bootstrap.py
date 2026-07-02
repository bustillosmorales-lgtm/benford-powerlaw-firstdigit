#!/usr/bin/env python3
"""
Parametric-bootstrap calibration of the corrected first-digit test, addressing
reviewer points 3.5 (report N*psi_corr, df, p-values), 3.6 (fixed thresholds
are not inferential rules), and 3.7 (a full parametric bootstrap, including the
saturated k_hat = inf rows).

For each dataset the fitted null is P_{k_hat} with the scale phase fixed at
phi = 0 by normalization, so only k is estimated and the asymptotic reference
is chi-squared with 9 - 1 - 1 = 7 degrees of freedom. The bootstrap repeats the
full estimate-and-test cycle, which is the honest calibration when k_hat
saturates (k_hat = inf, fitted null = Benford) or expected cell counts are
small.

Algorithm (per dataset):
  1. fitted null p0 = P_{k_hat}  (Benford if k_hat = inf)
  2. observed statistic  T_obs = N * psi(obs, base = p0) = N * psi_corr
  3. for b = 1..B:  draw counts ~ Multinomial(N, p0); re-estimate k on the
     replicate; recompute T_b = N * psi(freq_b, base = P_{k_hat_b})
  4. bootstrap p-value = (1 + #{T_b >= T_obs}) / (B + 1)
  5. compare with the asymptotic chi-squared_7 p-value

Run:  python bootstrap.py
"""
import numpy as np
from math import ceil
from scipy.stats import chi2
from scipy.optimize import minimize_scalar

DIGITS = np.arange(1, 10)
PB = np.log10(1 + 1 / DIGITS)
DF = 7                      # 9 cells - 1 (sum) - 1 (k estimated), phi fixed = 0
KMAX = 50.0                # numerical saturation bound; k_hat >= KMAX == "inf"


def Pk(k):
    r = 10.0 ** (1.0 / k)
    return r / (r - 1) * (DIGITS ** (-1.0 / k) - (DIGITS + 1) ** (-1.0 / k))


def psi(p, base):
    p = np.asarray(p, float)
    return float(np.sum((p - base) ** 2 / base))


def khat(obs):
    """Minimum-chi-squared exponent; returns (k, saturated?)."""
    f = lambda k: psi(obs, base=Pk(k)) if k > 1e-6 else 1e9
    res = minimize_scalar(f, bounds=(0.05, KMAX), method="bounded")
    return res.x, res.x >= KMAX - 1e-3


def fitted_null(k_hat, saturated):
    return PB if saturated else Pk(k_hat)


def bootstrap_row(N, k_hat, psi_corr, B, rng):
    saturated = (k_hat is None) or (not np.isfinite(k_hat)) or (k_hat >= KMAX - 1e-3)
    p0 = fitted_null(k_hat, saturated)
    T_obs = N * psi_corr
    # cap the simulated N for very large datasets: the null distribution of the
    # statistic depends on N only through cell-count discreteness, negligible
    # above ~2e5, so we simulate at min(N, NSIM) and scale T_obs is kept as given.
    NSIM = int(min(N, 200_000))
    T = np.empty(B)
    for b in range(B):
        counts = rng.multinomial(NSIM, p0)
        f = counts / NSIM
        kb, satb = khat(f)
        pb = fitted_null(kb, satb)
        T[b] = NSIM * psi(f, base=pb)
    p_boot = (1 + np.sum(T >= T_obs)) / (B + 1)
    p_asym = float(chi2.sf(T_obs, DF))
    return dict(N=N, k_hat=k_hat, psi_corr=psi_corr, Npsi=T_obs,
                p_boot=p_boot, p_asym=p_asym,
                null_mean=float(T.mean()), null_q95=float(np.quantile(T, 0.95)),
                saturated=saturated)


# ---- Table 4 rows as currently in the manuscript (tab:realdata) -------------
# (dataset, N, psi_classical, k_hat[inf->np.inf], psi_corr)
ROWS = [
    ("Simulated power law k=2", 2000,    0.107, 2.01,     0.004),
    ("Simulated uniform",       2000,    0.391, np.inf,   0.410),
    ("Election returns",        20272,   0.001, np.inf,   0.001),
    ("Financial statements",    3236586, 0.002, np.inf,   0.002),
    ("Physical constants",      355,     0.036, 7.71,     0.029),
    ("World Bank population",    265,     0.016, 40.0,     0.016),
    ("Astronomical catalogue",  109400,  0.018, 8.57,     0.012),
    ("Fibonacci scaled",        500,     0.000, np.inf,   0.001),
    ("Pruitt (retracted)",      10269,   0.374, np.inf,   0.386),
    ("Seismic catalogue",       500,     7.466, np.inf,   7.526),
]

B = 3000
rng = np.random.default_rng(20260702)

print(f"Parametric bootstrap, B={B}, df={DF}, phi=0 (scale-normalized).\n")
hdr = f"{'Dataset':<26}{'N':>9}{'k_hat':>7}{'psi_corr':>9}{'N*psi_corr':>11}{'p_boot':>9}{'p_asym':>9}{'null95':>8}"
print(hdr)
print("-" * len(hdr))
out = []
for name, N, psi_c, kh, psi_corr in ROWS:
    r = bootstrap_row(N, kh, psi_corr, B, rng)
    khs = "inf" if r["saturated"] else f"{kh:.2f}"
    print(f"{name:<26}{N:>9}{khs:>7}{psi_corr:>9.3f}{r['Npsi']:>11.1f}"
          f"{r['p_boot']:>9.3f}{r['p_asym']:>9.1e}{r['null_q95']:>8.1f}")
    out.append((name, r))

print("\nReading:")
print("- p_boot / p_asym near 0  => data incompatible with ANY power-law baseline")
print("  at that N (the calibrated test rejects); psi_corr is then read as an")
print("  effect size, small = close to the baseline.")
print("- null_q95 is the 95% point of N*psi_corr under the fitted null; compare")
print("  N*psi_corr against it. For huge N even tiny misspecification exceeds it,")
print("  which is exactly why verdicts must be effect-size based, not p-based.")
