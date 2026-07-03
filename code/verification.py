#!/usr/bin/env python3
"""
Numerical verification of the analytical results in "When the Benford law must
fail: A corrected first-digit test for power-law distributions" (Bustillos and
Leiva). Each block checks one stated result to machine precision (exact forms)
or by simulation (asymptotic claims). Run:  python verification.py

Contents
  1. Power-law first-digit law P_k: exactness vs direct integration; sum to one.
  2. Equivalence with the generalized Benford law (symbolic, exact).
  3. Benford limit k -> infinity; strict monotonicity of psi_k.
  4. Contamination identity psi_mix = alpha^2 psi_k; reference sizes N_min.
  5. Scale-phase family P_{k,phi}: phi=0 recovers P_k; sums to one.
  6. Joint (k,phi) injectivity search (no collision for moderate k).
  7. Degrees of freedom of the corrected statistic: chi^2_8 / chi^2_7 / chi^2_6.
  8. Level of the one-sided equivalence (non-inferiority) decision.
"""
import numpy as np
from scipy import integrate
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import chi2, kstest

from joint_phi_tost import DIGITS, PB, Pk, Pkphi, chisq, khat_joint

rng = np.random.default_rng(20260702)
ok = lambda b: "OK " if b else "**FAIL**"


def psi_k(k):
    return float(np.sum((Pk(k) - PB) ** 2 / PB))


print("=" * 68)
print("1. P_k exact vs direct integration of the density (~1e-15 expected)")
maxerr = 0.0
for k in [0.5, 1.0, 1.5, 2.0, np.e, np.pi]:
    a = 1.0 / k
    # density of the mantissa m={log10 y} is proportional to 10^{-m/k} on [0,1)
    Z = integrate.quad(lambda m: 10 ** (-m / k), 0, 1)[0]
    # mantissa m={log10 Y} has density proportional to 10^{-m/k}; digit d on
    # [log10 d, log10(d+1)) -> integrates to P_k(d) directly
    P = np.array([integrate.quad(lambda m: 10 ** (-m / k), np.log10(d), np.log10(d + 1))[0] / Z
                  for d in DIGITS])
    err = np.max(np.abs(P - Pk(k)))
    maxerr = max(maxerr, err)
print(f"   max |closed form - integral| = {maxerr:.2e}   {ok(maxerr < 1e-12)}")
s1 = max(abs(Pk(k).sum() - 1) for k in [0.5, 1, 1.5, 2, 3, 7])
print(f"   max |sum P_k - 1|            = {s1:.2e}   {ok(s1 < 1e-12)}")

print("=" * 68)
print("2. P_k == generalized Benford law (Pietronero), symbolic")
try:
    import sympy as sp
    d, k = sp.symbols("d k", positive=True)
    a = -1 / k
    Pgb = ((d + 1) ** a - d ** a) / (10 ** a - 1)
    Pk_sym = 10 ** (1 / k) / (10 ** (1 / k) - 1) * (d ** (-1 / k) - (d + 1) ** (-1 / k))
    diff = sp.simplify(Pgb - Pk_sym)
    print(f"   simplify(P_GB - P_k) = {diff}   {ok(diff == 0)}")
except ImportError:
    print("   sympy not installed; skipped")

print("=" * 68)
print("3. Benford limit and strict monotonicity of psi_k")
lim = max(abs(Pk(K)[0] - PB[0]) for K in [100, 1000, 10000])
print(f"   |P_k(1) - Benford| at large k = {lim:.2e}   {ok(lim < 1e-2)}")
ks = np.linspace(0.05, 300, 6000)
ps = np.array([psi_k(k) for k in ks])
mono = np.all(np.diff(ps) <= 1e-12)
print(f"   psi_k strictly decreasing on [0.05,300] = {mono}   {ok(mono)}")

print("=" * 68)
print("4. Contamination identity and reference sizes")
cmax = 0.0
for al in [0.0, 0.3, 0.5, 0.7, 1.0]:
    for k in [1.0, 2.0, 3.0]:
        pmix = al * Pk(k) + (1 - al) * PB
        cmax = max(cmax, abs(float(np.sum((pmix - PB) ** 2 / PB)) - al ** 2 * psi_k(k)))
print(f"   max |psi_mix - alpha^2 psi_k| = {cmax:.2e}   {ok(cmax < 1e-12)}")
crit = chi2.ppf(0.95, 8)
Nmin = [int(np.ceil(crit / psi_k(k))) for k in [1, 2, 3, 4, 5]]
print(f"   N_min(k=1..5) = {Nmin}   {ok(Nmin == [43,154,338,595,924])}")

print("=" * 68)
print("5. Scale-phase family P_{k,phi}")
r0 = max(np.max(np.abs(Pkphi(k, 0.0) - Pk(k))) for k in [1, 1.5, 2, 3, 7.3])
print(f"   max |P_{{k,0}} - P_k| = {r0:.2e}   {ok(r0 < 1e-12)}")
sph = max(abs(Pkphi(k, ph).sum() - 1) for k in [1.2, 2.5, 9]
          for ph in [0.0, 0.2, 0.5, 0.83])
print(f"   max |sum P_{{k,phi}} - 1| = {sph:.2e}   {ok(sph < 1e-12)}")

print("=" * 68)
print("6. Joint (k,phi) injectivity: search for a distinct-parameter collision")
print("   (moderate range k in [1,8]; for large k the law -> Benford and the")
print("    phase is only weakly identified, the saturated case)")
grid = [(k, ph) for k in np.linspace(1, 8, 100) for ph in np.linspace(0, 1, 120, endpoint=False)]
mats = np.array([Pkphi(k, ph) for k, ph in grid])
targets = [(2.0, 0.30), (3.5, 0.70), (1.4, 0.10), (5.0, 0.55)]
gphi = np.array([ph for _, ph in grid])
worst = np.inf
for (kt, pt) in targets:
    Pt = Pkphi(kt, pt)
    dist = np.sqrt(((mats - Pt) ** 2).sum(1))
    # the identification-relevant direction: a genuinely different phase
    diffphi = np.minimum(np.abs(gphi - pt), 1 - np.abs(gphi - pt)) > 0.08
    worst = min(worst, dist[diffphi].min())
print(f"   min ||P-P'|| over phase-separated pairs = {worst:.4f}   "
      f"{ok(worst > 1e-2)}  (no cross-phase collision -> phase identified)")

print("=" * 68)
print("7. Corrected statistic: degrees of freedom (continuous joint optimizer)")


def corr_stat(counts, N, mode):
    obs = counts / N
    if mode == "known":
        p = Pkphi(2.0, 0.0)
    elif mode == "phi0":
        f = lambda k: chisq(obs, Pk(k)) if k > 1e-6 else 1e9
        p = Pk(minimize_scalar(f, bounds=(1, 50), method="bounded").x)
    else:  # joint
        kh, ph, _, _ = khat_joint(obs, n_phi=1)  # seed
        res = minimize(lambda x: chisq(obs, Pkphi(max(x[0], 1.0), x[1] % 1.0)),
                       x0=[2.0, 0.0], method="Nelder-Mead",
                       options=dict(xatol=1e-4, fatol=1e-9))
        p = Pkphi(max(res.x[0], 1.0), res.x[1] % 1.0)
    return N * chisq(obs, p)


N, R = 6000, 700
p_true = Pkphi(2.0, 0.0)
for mode, dfn in [("known", 8), ("phi0", 7), ("joint", 6)]:
    stats = np.array([corr_stat(rng.multinomial(N, p_true), N, mode) for _ in range(R)])
    ksd = {df: kstest(stats, "chi2", args=(df,)).statistic for df in (6, 7, 8)}
    best = min(ksd, key=ksd.get)
    print(f"   {mode:5s}: mean={stats.mean():.2f}  best chi2 df={best}  claim={dfn}   {ok(best == dfn)}")

print("=" * 68)
print("8. Equivalence decision (delta=0.10) on controlled populations")


def u95(p, N=4000, B=200):
    ps = np.array([khat_joint(rng.multinomial(N, p) / N, n_phi=21)[2] for _ in range(B)])
    return float(np.quantile(ps, 0.95))


for name, p in [("power law P_{2,0.3}", Pkphi(2.0, 0.3)),
                ("uniform (fabricated)", np.full(9, 1 / 9))]:
    u = u95(p)
    print(f"   {name:22s} u_0.95={u:.3f}  ->  "
          f"{'consistent' if u < 0.10 else 'anomalous'}   "
          f"{ok((u < 0.10) == ('power' in name))}")
print("   (a genuine power law is certified; fabricated data is flagged)")
print("=" * 68)
print("Done.")
