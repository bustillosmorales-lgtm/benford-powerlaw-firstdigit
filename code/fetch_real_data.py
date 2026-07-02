#!/usr/bin/env python3
"""
Fetch the publicly-accessible validation datasets from their stable APIs and
recompute the real first-digit structure index. Documents the exact query for
each so the numbers are reproducible. Forensic-case datasets behind portals
(Pruitt/Stapel/Gino/Schroeder) are not fetched here.
"""
import json, urllib.request, re
import numpy as np
from math import log10
from scipy.optimize import minimize_scalar

DIGITS = np.arange(1, 10)
PB = np.log10(1 + 1 / DIGITS)
UA = {"User-Agent": "Mozilla/5.0 (research reproducibility script)"}


def get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60).read().decode("utf-8", "ignore")


def Pk(k):
    r = 10 ** (1.0 / k)
    return r / (r - 1) * (DIGITS ** (-1.0 / k) - (DIGITS + 1) ** (-1.0 / k))


def psi(p, base=PB):
    p = np.asarray(p, float)
    return float(np.sum((p - base) ** 2 / base))


def first_digit(vals):
    v = np.abs(np.asarray(vals, float))
    v = v[(v > 0) & np.isfinite(v)]
    return (v / 10.0 ** np.floor(np.log10(v))).astype(int)


def analyze(name, vals, query):
    d = first_digit(vals)
    obs = np.array([(d == j).mean() for j in DIGITS])
    cl = psi(obs)
    res = minimize_scalar(lambda k: psi(obs, Pk(k)) if k > 1e-6 else 1e9,
                          bounds=(0.05, 50), method="bounded")
    kh, corr = res.x, res.fun
    print(f"\n{name}: N={len(d)}")
    print(f"  query: {query}")
    print(f"  classical psi = {cl:.4f}   k_hat = {kh:.3f}   psi_corr = {corr:.4f}")
    return cl, kh, corr, len(d)


# ---- seismic (USGS ComCat) : first digit of magnitudes -------------------
try:
    q = "fdsnws M>=4, 2024-01, first 500"
    raw = get("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv"
              "&starttime=2024-01-01&endtime=2024-01-31&minmagnitude=4&limit=500")
    rows = [r.split(",") for r in raw.strip().splitlines()[1:]]
    mags = [float(r[4]) for r in rows if r[4]]
    analyze("Seismic catalogue (magnitudes)", mags, q)
    # also depth, which spans more orders of magnitude
    deps = [float(r[3]) for r in rows if r[3]]
    analyze("Seismic catalogue (depths)", deps, q + ", depth column")
except Exception as e:
    print("seismic FAILED:", e)

# ---- physical constants (NIST CODATA 2022) -------------------------------
try:
    raw = get("https://physics.nist.gov/cuu/Constants/Table/allascii.txt")
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
    analyze("Physical constants (CODATA 2022 values)", vals, "NIST allascii.txt, all numeric values")
except Exception as e:
    print("constants FAILED:", e)

# ---- World Bank population (SP.POP.TOTL, 2022) ---------------------------
try:
    raw = get("https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL"
              "?date=2022&format=json&per_page=400")
    data = json.loads(raw)[1]
    pops = [d["value"] for d in data if d["value"]]
    analyze("World Bank population (2022)", pops, "SP.POP.TOTL date=2022 all economies")
except Exception as e:
    print("worldbank FAILED:", e)

# ---- astronomical catalogue (HYG) : distances ---------------------------
for url in ["https://raw.githubusercontent.com/astronexus/HYG-Database/master/hyg/v3/hyg_v37.csv",
            "https://codeberg.org/astronexus/hyg/raw/branch/main/data/hyg/CURRENT/hygdata_v41.csv"]:
    try:
        raw = get(url)
        lines = raw.splitlines()
        hdr = lines[0].split(",")
        di = hdr.index("dist") if "dist" in hdr else None
        if di is None:
            continue
        ds = []
        for ln in lines[1:]:
            parts = ln.split(",")
            try:
                x = float(parts[di])
                if 0 < x < 1e5:        # drop the 100000 placeholder
                    ds.append(x)
            except (ValueError, IndexError):
                pass
        analyze("Astronomical catalogue (HYG distances)", ds, url.split('/')[-1] + ", dist column")
        break
    except Exception as e:
        print("HYG try FAILED:", e)
