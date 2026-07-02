#!/usr/bin/env python3
"""
Reproducibility script for "When the Benford law must fail:
A corrected first-digit test for power-law distributions" (Bustillos & Leiva).

Implements the closed-form power-law first-digit law P_k(d), the structure
index psi, the minimum-chi-squared exponent estimator k-hat, the corrected
index psi_corr, the reference sample size N_min(k), and the contamination
scaling. Reproduces the analytical reference tables and the deterministic /
simulated validation rows. Real-data forensic rows require their own
externally downloaded datasets (see data/ and the README); this script
covers everything that is reproducible without a third-party download.

Run:  python benford_powerlaw.py
"""
import numpy as np
from math import log10, ceil
from scipy.stats import chi2
from scipy.optimize import minimize_scalar

DIGITS = np.arange(1, 10)
PB = np.log10(1 + 1 / DIGITS)            # Benford first-digit law
CHI2_05_8 = chi2.ppf(0.95, 8)            # 15.507...
CHI2_01_8 = chi2.ppf(0.99, 8)


def Pk(k):
    """Power-law first-digit distribution, Theorem (power-law first-digit)."""
    r = 10 ** (1.0 / k)
    return r / (r - 1) * (DIGITS ** (-1.0 / k) - (DIGITS + 1) ** (-1.0 / k))


def psi(p, base=PB):
    """Structure index: chi-squared divergence of distribution p from base."""
    p = np.asarray(p, float)
    return float(np.sum((p - base) ** 2 / base))


def psi_k(k):
    return psi(Pk(k))


def Nmin(k):
    return ceil(CHI2_05_8 / psi_k(k))


def first_digits(values):
    """Leading significant digit of |x| for x != 0."""
    v = np.abs(np.asarray(values, float))
    v = v[v > 0]
    d = (v / 10.0 ** np.floor(np.log10(v))).astype(int)
    return d


def observed_freq(values):
    d = first_digits(values)
    counts = np.array([(d == j).sum() for j in DIGITS], float)
    return counts / counts.sum(), len(d)


def khat(obs):
    """Minimum-chi-squared estimator of the exponent."""
    f = lambda k: psi(obs, base=Pk(k)) if k > 1e-6 else 1e9
    res = minimize_scalar(f, bounds=(0.05, 50), method="bounded")
    return res.x, res.fun


def report(title):
    print("\n" + "=" * 66 + "\n" + title + "\n" + "=" * 66)


# ---------------------------------------------------------------- analytics
report("Table: analytical reference values (N = 1000)")
print(f"{'k':>4} {'P_k(1)':>9} {'psi_k':>9} {'N*psi_k':>9}")
for k in [1, 2, 3, 4, 5]:
    p = Pk(k)
    print(f"{k:>4} {p[0]:>9.3f} {psi_k(k):>9.3f} {1000*psi_k(k):>9.1f}")
print(f"{'inf':>4} {PB[0]:>9.3f} {0.0:>9.3f} {0.0:>9.1f}")

report("Table: reference sample sizes N_min(k) = ceil(15.51/psi_k)")
for k in [1, 2, 3, 4, 5]:
    print(f"  k={k}: psi_k={psi_k(k):.4f}  N_min={Nmin(k)}")
print(f"  check: chi2_0.05,8 = {CHI2_05_8:.3f}")

# ------------------------------------------------------------- Fibonacci row
report("Validation row: Fibonacci scaled (deterministic, N=500)")
fib = [1, 1]
while len(fib) < 500:
    fib.append(fib[-1] + fib[-2])
obs, n = observed_freq(fib[:500])
print(f"  N={n}  classical psi (vs Benford) = {psi(obs):.4f}")
kh, _ = khat(obs)
print(f"  k_hat={kh:.3f}  psi_corr={psi(obs, base=Pk(kh)):.4f}")
print("  (Fibonacci follows Benford closely -> classical psi small -> clean)")

# ---------------------------------------------------- simulated baselines
report("Validation rows: simulated baselines (seed=20260316)")
rng = np.random.default_rng(20260316)

# power-law baseline k=2, N=2000: draw first digits from P_k(2)
p2 = Pk(2)
draw = rng.choice(DIGITS, size=2000, p=p2)
obs2 = np.array([(draw == j).mean() for j in DIGITS])
kh2, _ = khat(obs2)
print(f"  Simul k=2 (N=2000): classical psi={psi(obs2):.4f}  "
      f"k_hat={kh2:.3f}  psi_corr={psi(obs2, base=Pk(kh2)):.4f}")

# fabricated/uniform baseline N=2000: flat first digits
draw_u = rng.choice(DIGITS, size=2000, p=np.full(9, 1/9))
obs_u = np.array([(draw_u == j).mean() for j in DIGITS])
khu, _ = khat(obs_u)
print(f"  Simul uniform (N=2000): classical psi={psi(obs_u):.4f}  "
      f"k_hat={khu:.3f}  psi_corr={psi(obs_u, base=Pk(khu)):.4f}")

# ------------------------------------------------------------- contamination
report("Contamination scaling check: psi_mix = alpha^2 * psi_k (population)")
k = 2
for a in [0.25, 0.5, 0.75, 1.0]:
    pmix = a * Pk(k) + (1 - a) * PB
    print(f"  alpha={a:>4}: psi_mix={psi(pmix):.5f}  alpha^2*psi_k={a**2*psi_k(k):.5f}")

print("\nDone. Worked example N_min(2) =", Nmin(2))
