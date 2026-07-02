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

## Scale normalization (phi = 0)

The corrected test is applied to scale-normalized data. Each dataset is divided
by `10 ** round(mean(log10(|x|)))` before the digits are read, which fixes the
scale phase `phi = 0` and leaves the first-digit deviations unchanged. Only the
exponent `k` is then estimated, so the corrected test carries seven degrees of
freedom.

## Determinism

`benford_powerlaw.py` uses seed `20260316`; `bootstrap.py` uses seed
`20260702`. Both are fixed in the source, so their output is bit-reproducible.
