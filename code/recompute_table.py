#!/usr/bin/env python3
"""
Full recomputation of Table 4 with the JOINT (k, phi) minimum-chi-squared
estimator (6 df), the parametric bootstrap p-value, and a TOST equivalence
decision. Loads every dataset's observed first-digit vector from the archives in
../repro2 (or regenerates seeded/API data), so every row is reproducible.

Outputs one row per dataset:
  N, psi (classical), k0/psic0 (old phi=0 fit), khat/phihat/psi_corr (joint, 6 df),
  N*psi_corr, p_boot (6 df), psi_corr upper 95% bound, TOST verdict at delta.
"""
import io, os, gzip, zipfile, json, urllib.request
import numpy as np
from scipy.optimize import minimize_scalar

from joint_phi_tost import (DIGITS, PB, KMAX, Pk, Pkphi, chisq,
                            khat_phi0, khat_joint)

REPRO2 = os.path.join(os.path.dirname(__file__), "..", "repro2")
UA = {"User-Agent": "Mozilla/5.0 (research reproducibility)"}
DELTA = 0.10          # TOST equivalence margin on psi_corr (effect size)
ALPHA = 0.05
B = 600
NPHI_BOOT = 46


def first_digits(vals):
    v = np.abs(np.asarray(vals, float))
    v = v[(v > 0) & np.isfinite(v)]
    d = (v / 10.0 ** np.floor(np.log10(v))).astype(int)
    return d[(d >= 1) & (d <= 9)]


def obs_from(vals):
    d = first_digits(vals)
    c = np.array([(d == j).sum() for j in DIGITS], float)
    return c / c.sum(), int(c.sum()), c


# ----------------------------------------------------------------- loaders
def load_sim_powerlaw():
    rng = np.random.default_rng(20260316)
    draw = rng.choice(DIGITS, size=2000, p=Pk(2))
    c = np.array([(draw == j).sum() for j in DIGITS], float)
    return c / c.sum(), 2000, c


def load_sim_uniform():
    rng = np.random.default_rng(20260316)
    rng.choice(DIGITS, size=2000, p=Pk(2))            # advance to match sequence
    draw = rng.choice(DIGITS, size=2000, p=np.full(9, 1 / 9))
    c = np.array([(draw == j).sum() for j in DIGITS], float)
    return c / c.sum(), 2000, c


def load_fibonacci():
    fib = [1, 1]
    while len(fib) < 500:
        fib.append(fib[-1] + fib[-2])
    return obs_from(fib[:500])


def get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90).read()


def load_constants():
    import re
    raw = get("https://physics.nist.gov/cuu/Constants/Table/allascii.txt").decode("utf-8", "ignore")
    vals = []
    for ln in raw.splitlines():
        m = re.match(r"^.{60}(.{25})", ln)
        if not m:
            continue
        v = m.group(1).replace(" ", "").replace("...", "")
        if re.match(r"^[-+]?\d", v):
            try:
                vals.append(float(v.replace("e", "E")))
            except ValueError:
                pass
    return obs_from(vals)


def load_worldbank():
    raw = get("https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL"
              "?date=2022&format=json&per_page=400").decode("utf-8")
    data = json.loads(raw)[1]
    return obs_from([d["value"] for d in data if d["value"]])


def load_seismic():
    raw = get("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv"
              "&starttime=2024-01-01&endtime=2024-01-31&minmagnitude=4&limit=500").decode("utf-8")
    rows = [r.split(",") for r in raw.strip().splitlines()[1:]]
    return obs_from([float(r[4]) for r in rows if r[4]])


def load_election():
    import pandas as pd
    with zipfile.ZipFile(os.path.join(REPRO2, "dataverse_files.zip")) as z:
        with z.open("countypres_2000-2024.csv") as f:
            df = pd.read_csv(f, usecols=["candidatevotes", "year"], low_memory=False)
    v = df[df["year"] == 2024]["candidatevotes"].to_numpy()
    return obs_from(v)


def load_financial():
    """Stream the 'value' column (index 8) of SEC num.txt from the zip."""
    path = os.path.join(REPRO2, "2026q1.zip")
    counts = np.zeros(10, float)
    with zipfile.ZipFile(path) as z:
        with z.open("num.txt") as raw:
            for i, bline in enumerate(io.BufferedReader(raw, 1 << 20)):
                if i == 0:
                    continue
                parts = bline.split(b"\t")
                if len(parts) < 9:
                    continue
                s = parts[8].strip()
                if not s:
                    continue
                try:
                    x = abs(float(s))
                except ValueError:
                    continue
                if x > 0 and np.isfinite(x):
                    d = int(x / 10.0 ** np.floor(np.log10(x)))
                    if 1 <= d <= 9:
                        counts[d] += 1
    c = counts[1:]
    return c / c.sum(), int(c.sum()), c


def load_astro():
    import pandas as pd
    with zipfile.ZipFile(os.path.join(REPRO2, "HYG-Database-main.zip")) as z:
        name = "HYG-Database-main/hyg/CURRENT/hygdata_v40.csv.gz"
        with z.open(name) as f:
            df = pd.read_csv(gzip.open(io.BytesIO(f.read())), usecols=["dist"], low_memory=False)
    d = df["dist"].to_numpy()
    d = d[(d > 0) & (d < 1e5)]
    return obs_from(d)


def load_pruitt():
    import pandas as pd
    with zipfile.ZipFile(os.path.join(REPRO2, "doi_10_5061_dryad_33f0n__v20200121.zip")) as z:
        data = z.read("Laskowski+et+al_social+niche+disruption_DATA.xlsx")
    xl = pd.ExcelFile(io.BytesIO(data))
    vals = []
    for sh in xl.sheet_names:
        df = xl.parse(sh)
        for col in df.columns:
            s = pd.to_numeric(df[col], errors="coerce").to_numpy()
            s = s[np.isfinite(s)]
            # continuous columns: many distinct non-integer values
            if s.size and (np.mean(np.abs(s - np.round(s)) > 1e-9) > 0.3):
                vals.extend(s.tolist())
    return obs_from(vals)


LOADERS = [
    ("Simulated power law k=2", load_sim_powerlaw),
    ("Simulated uniform",       load_sim_uniform),
    ("Election returns",        load_election),
    ("Financial statements",    load_financial),
    ("Physical constants",      load_constants),
    ("World Bank population",    load_worldbank),
    ("Astronomical catalogue",  load_astro),
    ("Fibonacci scaled",        load_fibonacci),
    ("Pruitt (retracted)",      load_pruitt),
    ("Seismic catalogue",       load_seismic),
]


def bootstrap_and_tost(counts, N, rng):
    """Nonparametric bootstrap of psi_corr (joint fit) for the equivalence CB,
    and a parametric-null p_boot for the 6-df goodness-of-fit reference."""
    obs = counts / counts.sum()
    khat, phihat, psic, sat = khat_joint(obs)
    # nonparametric bootstrap: resample observed counts -> CI of psi_corr
    psis = np.empty(B)
    p0 = obs
    for b in range(B):
        cc = rng.multinomial(N, p0) / N
        _, _, psis[b], _ = khat_joint(cc, n_phi=NPHI_BOOT)
    ucb = float(np.quantile(psis, 1 - ALPHA))         # upper 95% bound on psi_corr
    # parametric null p_boot (fit -> simulate -> refit) with 6-df reference
    fit = Pkphi(khat, phihat) if not sat else PB
    NSIM = int(min(N, 200000))
    Tnull = np.empty(B)
    for b in range(B):
        cc = rng.multinomial(NSIM, fit) / NSIM
        _, _, s, _ = khat_joint(cc, n_phi=NPHI_BOOT)
        Tnull[b] = NSIM * s
    p_boot = (1 + np.sum(Tnull >= N * psic)) / (B + 1)
    tost = "equivalent" if ucb < DELTA else "not-equiv"
    return khat, phihat, psic, sat, ucb, p_boot, tost


if __name__ == "__main__":
    rng = np.random.default_rng(20260702)
    hdr = (f"{'Dataset':<26}{'N':>9}{'psi':>8}{'k0':>7}{'psic0':>8}"
           f"{'khat':>7}{'phih':>6}{'psi_c':>8}{'Npsic':>9}{'p6':>7}{'ucb':>7}{'TOST':>11}")
    print(hdr); print("-" * len(hdr))
    for name, fn in LOADERS:
        try:
            obs, N, counts = fn()
        except Exception as e:
            print(f"{name:<26} LOAD FAILED: {e}")
            continue
        psi = chisq(obs, PB)
        k0, c0, s0 = khat_phi0(obs)
        kh, ph, psic, sat, ucb, p6, tost = bootstrap_and_tost(counts, N, rng)
        k0s = "inf" if s0 else f"{k0:.2f}"
        khs = "inf" if sat else f"{kh:.2f}"
        print(f"{name:<26}{N:>9}{psi:>8.3f}{k0s:>7}{c0:>8.3f}"
              f"{khs:>7}{ph:>6.2f}{psic:>8.3f}{N*psic:>9.1f}{p6:>7.3f}{ucb:>7.3f}{tost:>11}")
