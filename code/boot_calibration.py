#!/usr/bin/env python3
"""
Bootstrap equivalence-decision calibration study (referee point 4.2a), at paper
resolution. Fast, grid-vectorised: P_{k,phi} is precomputed on a fine (k,phi) grid
(200 x 120) and the minimum-chi-squared fit is a vectorised min over that grid, so
each bootstrap of B resamples is a single array reduction. The grid step (~0.25 in k,
~0.008 in phi) is fine enough that the grid minimum tracks the continuous optimiser;
coarsening the grid would bias psi_corr up and flatter the level, so it is kept fine.

Checks, for the decision "declare CONSISTENT iff u_0.95 < delta = 0.10":
  (A) LEVEL   : populations with true psi_corr = delta (H0 boundary) -> P(consistent) <= alpha
  (B) interior: psi_corr > delta -> P(consistent) ~ 0
  (C) POWER   : genuine power law (psi_corr ~ 0) -> P(consistent) ~ 1
  (D) boundary k=1 and saturation (Benford) -> decision well behaved
  (E) COVERAGE: empirical coverage of the upper bound u_0.95 for psi_corr
Run:  python boot_calibration.py
"""
import numpy as np

DIGITS = np.arange(1, 10)
PB = np.log10(1 + 1.0 / DIGITS)
L_LO = np.log10(DIGITS); L_HI = np.log10(DIGITS + 1)
rng = np.random.default_rng(20260704)
DELTA, ALPHA, N = 0.10, 0.05, 2000

def Pkphi(k, phi):
    beta = -1.0 / k
    def Ghat(x):
        fx = np.floor(x)
        return fx + (10.0 ** (beta * (x - fx)) - 1.0) / (10.0 ** beta - 1.0)
    return Ghat(L_HI - phi) - Ghat(L_LO - phi)

# ---- precomputed fine (k,phi) grid for a vectorised min-chi-squared fit ----
KG = np.linspace(1.0, 50.0, 200)
PG = np.linspace(0.0, 1.0, 120, endpoint=False)
PMAT = np.array([[Pkphi(k, ph) for ph in PG] for k in KG]).reshape(-1, 9)  # (G,9)
PMAT = np.clip(PMAT, 1e-12, None)
INVP = (1.0 / PMAT).T          # (9,G); chi^2 = sum_d res_d^2/P_d - 1 for prob vectors

def psi_grid(obs):
    """min-chi-squared distance of a probability vector obs to the family (grid)."""
    return float(((obs * obs) @ INVP).min() - 1.0)

def u95_boot(counts, B=600):
    """vectorised nonparametric multinomial bootstrap upper 95% bound on psi_corr.
    Uses chi^2(res, P) = sum_d res_d^2/P_d - 1 (res and P both sum to 1), so the
    fit over the grid is a single (B x 9)(9 x G) matmul -- low memory, fast."""
    obs = counts / counts.sum()
    res = rng.multinomial(N, obs, size=B) / N            # (B,9)
    psi = ((res * res) @ INVP).min(1) - 1.0              # (B,)
    return float(np.quantile(psi, 0.95))

def decide_rate(q, R=400, B=400):
    return np.mean([u95_boot(rng.multinomial(N, q), B) < DELTA for _ in range(R)])

def make_pop(kb, pb, target):
    base = Pkphi(kb, pb)
    v = np.array([1, -1, 1, -1, 1, -1, 1, -1, 0.0]); v -= v.mean(); v /= np.linalg.norm(v)
    lo, hi = 0.0, 0.5
    for _ in range(60):
        mid = (lo + hi) / 2; q = base + mid * v
        if np.any(q <= 1e-6): hi = mid; continue
        q = q / q.sum()
        if psi_grid(q) < target: lo = mid
        else: hi = mid
    q = np.clip(base + lo * v, 1e-6, None); return q / q.sum()

print(f"N={N} delta={DELTA} alpha={ALPHA} grid=200x120 B=400/600 R=400 seed=20260704")
print("=" * 66)
print("(A) LEVEL  P(consistent | true psi_corr = delta)   target <= 0.05")
for kb, pb in [(2.0, 0.3), (1.5, 0.6), (3.0, 0.0)]:
    q = make_pop(kb, pb, DELTA)
    print(f"    base k={kb},phi={pb}: psi_corr={psi_grid(q):.4f}  "
          f"level={decide_rate(q):.3f}", flush=True)
print("(B) interior H0  psi_corr>delta -> P(consistent) ~ 0")
q = make_pop(2.0, 0.3, 0.18)
print(f"    psi_corr={psi_grid(q):.4f}  P(consistent)={decide_rate(q):.3f}", flush=True)
print("(C) POWER  genuine power law psi_corr~0 -> P(consistent) ~ 1")
for kb, pb in [(2.0, 0.3), (4.0, 0.7)]:
    q = Pkphi(kb, pb)
    print(f"    k={kb},phi={pb}: psi_corr={psi_grid(q):.4f}  "
          f"P(consistent)={decide_rate(q):.3f}", flush=True)
print("(D) boundary k=1 / saturation Benford")
for nm, q in [("k=1,phi=0.6", Pkphi(1.0, 0.6)), ("Benford k=inf", PB.copy())]:
    print(f"    {nm}: psi_corr={psi_grid(q):.4f}  "
          f"P(consistent)={decide_rate(q):.3f}", flush=True)
print("(E) COVERAGE of u_0.95 for psi_corr (nominal 0.95)")
for kb, pb in [(2.0, 0.3), (1.5, 0.6)]:
    q = Pkphi(kb, pb); psi_true = psi_grid(q)
    cov = np.mean([u95_boot(rng.multinomial(N, q), 600) >= psi_true for _ in range(400)])
    print(f"    k={kb},phi={pb}: psi_true={psi_true:.4f}  coverage={cov:.3f}", flush=True)
print("=" * 66); print("Done.")
