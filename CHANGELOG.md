# Changelog

All notable changes to **spmixw** are documented in this file. Versions
follow [Semantic Versioning](https://semver.org/).

## [0.1.0] — 2026-05-08

Initial release.

### Models

- `spmixw, model(ols)` — Bayesian fixed-effects OLS panel
- `spmixw, model(sar)` — spatial autoregressive
- `spmixw, model(sem)` — spatial error
- `spmixw, model(sdm)` — spatial Durbin
- `spmixw, model(sdem)` — spatial Durbin error
- `spmixw, model(slx)` — spatial lag of X
- `spmixw, model(sar_conv)` — SAR with convex combination of M weight matrices
- `spmixw, model(sem_conv)` — SEM convex
- `spmixw, model(sdm_conv)` — SDM convex
- `spmixw, model(sdem_conv)` — SDEM convex

### Postestimation

- `estat impact` — LeSage-Pace direct / indirect / total effects for
  SAR, SDM, SDEM, SLX and their `_conv` counterparts. Stochastic
  cross-trace estimator for the SDM-conv decomposition.

### Bayesian Model Averaging

- `spmixw_bma` — averages over the 2^M − 1 non-empty subsets of M
  candidate weight matrices, returning posterior model probabilities,
  per-subset coefficients, and BMA-weighted point estimates.

### Other features

- Heteroscedasticity-robust posterior via Geweke (1993) student-t errors
  (`rval()` option) for non-convex models.
- Informative priors on β via `bprior()` and `bvar()`.
- Region, time, and twoway fixed effects.
- Bayesian column headers (Post. Mean / Post. SD / Cred. Interval) with
  the spatial parameter rendered as a labelled row in the coefficient
  table (`Wy` for SAR/SDM, `We` for SEM/SDEM, `gam_m` for convex
  weights).

### Validation

- 120 validation checks pass on Stata 13: 9 OLS + 25 single-W spatial +
  35 convex + 26 BMA + 7 convex-impact + 18 misc.
- End-to-end demo (`examples/spmixw_demo.do`) generates a 50-region
  synthetic panel, builds three distinct W matrices, and exercises every
  user-facing command.

### Requirements

- Stata 13.0 or later. No Stata 14+ features required.
