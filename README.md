# When the Benford law must fail — reproducibility repository

Code and documentation reproducing every computed value in

> F. Bustillos and V. Leiva, *When the Benford law must fail: A corrected
> first-digit test for power-law distributions*, Chilean Journal of Statistics
> (ChJS), 17(1), 2026.

The paper proves that the classical chi-squared Benford test rejects authentic
inverse power-law data (the false positive theorem) and supplies a corrected,
calibrated first-digit test built on the power-law first-digit law
`P_k(d) = (10^{1/k}/(10^{1/k}-1)) (d^{-1/k} - (d+1)^{-1/k})`.

## Contents

```
code/
  benford_powerlaw.py   analytical tables, N_min(k), contamination law,
                        seeded simulated baselines, Fibonacci case
  joint_phi_tost.py     scale-phase law P_{k,phi}, the JOINT (k,phi)
                        minimum-chi-squared estimator (k>=1, 6 df), and the
                        TOST equivalence decision
  recompute_table.py    loads every dataset from ../repro2, runs the joint fit,
                        bootstrap and TOST -> the forensic Table (Table 4)
  bootstrap.py          earlier phi=0 parametric-bootstrap calibration (df=7)
  fetch_real_data.py    downloads the public validation datasets by API and
                        recomputes their first-digit indices
docs/
  QUERIES.md            exact source, URL, date and query for each dataset
  CLEANING.md           first-digit extraction rule; handling of zeros,
                        negatives and missing values
  DATA_DICTIONARY.md    variables used from each dataset
  CHECKSUMS.md          SHA-256 of the raw downloaded archives
results/
  RESULTS_bootstrap.md      calibrated Table (matches paper Table 4)
  RESULTS_real_recompute.md first-digit recomputation on the real data
requirements.txt
LICENSE
```

## Reproduce

```bash
python -m pip install -r requirements.txt
python code/benford_powerlaw.py     # analytical + simulated rows (fixed seed)
python code/bootstrap.py            # calibrated Table 4 (fixed seed 20260702)
python code/fetch_real_data.py      # re-downloads public data and recomputes
```

`benford_powerlaw.py` and `bootstrap.py` are fully self-contained and
deterministic (fixed seeds). `fetch_real_data.py` reaches live public APIs, so
its exact counts depend on the archive version on the day it is run; the
versions used in the paper are pinned by the checksums in `docs/CHECKSUMS.md`
and the queries in `docs/QUERIES.md`.

## Raw data

The large public archives (SEC EDGAR, HYG, Dataverse, Dryad, OSF) are **not**
committed here; they are redistributed by their maintainers under their own
terms. `docs/QUERIES.md` gives the exact retrieval for each and
`docs/CHECKSUMS.md` records the SHA-256 of the archives used, so any download
can be verified against what the paper used.

## License

Code released under the MIT License (`LICENSE`). The public datasets remain
under the licenses of their respective providers (see `docs/QUERIES.md`).
