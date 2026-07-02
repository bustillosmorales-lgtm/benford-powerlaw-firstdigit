# Cleaning rules and first-digit extraction

## First significant digit

For a value `x`, the leading significant digit is

```python
d = int(abs(x) / 10 ** floor(log10(abs(x))))     # d in {1,...,9}
```

implemented in `first_digit()` / `first_digits()` in the scripts. This is
scale-free: it depends only on the mantissa of `|x|`, not on its magnitude or
sign.

## Inclusion / exclusion

Applied uniformly to every dataset before extracting digits:

- **Zeros** are dropped (`x > 0` after taking the absolute value): `log10(0)`
  is undefined and a zero has no significant digit.
- **Negatives** are folded to their absolute value (`abs(x)`), then treated
  like any positive value; the sign carries no first-digit information.
- **Missing / non-finite** values (`NaN`, `inf`, empty fields) are dropped
  (`np.isfinite`).
- No rounding, binning or truncation is applied to the retained values.

`N` in the paper's tables is the count of retained values after these rules.

## Scale phase (jointly estimated, k >= 1)

The corrected test does NOT normalize the scale (dividing by an integer power of
ten leaves every first digit unchanged, so it cannot fix the phase). Instead the
exponent `k` and the scale phase `phi = {log10 C}` are estimated jointly by
minimum chi-squared against the scale-phase family `P_{k,phi}`, over the
physical range `k >= 1`. The `k >= 1` restriction is essential: with `phi` free
and `k` unrestricted, the two-parameter family over-fits the nine digit cells (a
degenerate small `k` with a phase mimics bounded-support data). Estimating both
parameters gives 9 - 1 - 2 = 6 degrees of freedom. See `joint_phi_tost.py`.

## Forensic verdict (TOST equivalence)

Because the corrected statistic is a consistent goodness-of-fit statistic, a
point-null p-value rejects any large sample. The verdict is instead an
equivalence decision (TOST): the data are consistent with a power law when the
upper 95% nonparametric-bootstrap bound on `psi_corr` falls below the margin
`delta = 0.10`, and anomalous otherwise.

## Determinism

`benford_powerlaw.py` uses seed `20260316`; `bootstrap.py` uses seed
`20260702`. Both are fixed in the source, so their output is bit-reproducible.
