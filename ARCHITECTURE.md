# spmixw (Stata) — Architecture

Design notes. These are open decisions until they're locked in Phase 1.
Everything here should be revisited before Phase 2 begins.

## 1. Command surface

**Proposal:** single top-level command `spmixw` with a `model()` option.

```stata
spmixw depvar indepvars [if] [in], model(name) [options]
```

`model()` values: `ols`, `sar`, `sem`, `sdm`, `sdem`, `slx`,
`sar_conv`, `sem_conv`, `sdm_conv`, `sdem_conv`,
`sar_bma`, `sdm_bma`, `sdem_bma`.

Internally, `spmixw.ado` parses common options (panel, MCMC, priors) and
dispatches to a model-specific subroutine `_spmixw_<model>.ado` or directly to
a Mata function. The subroutines share a common Mata MCMC engine.

Alternative: one ado per model (`spmixw_sar`, `spmixw_sdm`, ...). Closer to
`xsmle`'s style, easier help-file separation, but more parsing duplication.

**Decision deferred to Phase 1.**

## 2. Mata layer

The numerical core is **Mata**, compiled into `build/lspmixw.mlib`. Ado is the
parsing/return shell.

Public Mata functions (prefix `_spmixw_`):

```
_spmixw_demean()           // within-transformation
_spmixw_ols()              // pooled / FE OLS
_spmixw_logdet_exact()     // ln|I - rho W| via Cholesky
_spmixw_logdet_mc()        // Monte Carlo log-det
_spmixw_logdet_taylor()    // 4th-order Taylor for convex models
_spmixw_griddy_rho()       // griddy Gibbs sampler for rho
_spmixw_mh_gamma()         // Metropolis-Hastings for convex weights
_spmixw_lmarg()            // log-marginal likelihood
_spmixw_sar_mcmc()         // SAR MCMC main loop (also used by SDM)
_spmixw_sem_mcmc()         // SEM MCMC main loop (also used by SDEM)
_spmixw_sar_conv_mcmc()    // convex SAR (also used by SDM convex)
_spmixw_sem_conv_mcmc()    // convex SEM (also used by SDEM convex)
_spmixw_bma_enumerate()    // enumerate W subsets, compute posterior probs
_spmixw_impacts()          // direct / indirect / total effects
```

## 3. Sparse W storage

Mata's native sparse support is thin. **Plan:**

- Accept W as a Stata matrix name (or names, for convex models).
- Inside the Mata routine, convert to a sparse triplet representation:
  three vectors `i`, `j`, `v`. Provide helper Mata functions:
  - `_spmixw_sparse_mvmult(i, j, v, x)` — sparse × dense multiply
  - `_spmixw_sparse_quadform(i, j, v, x)` — x'Wx
- For Cholesky-based exact log-det, fall back to dense `cholesky()` if W is
  small enough; revisit for large W in a later release.

If Mata performance is unacceptable, escalate to a C plugin. Out of scope for
v0.1.

## 4. RNG and reproducibility

- Stata's `set seed N` should fully control all random draws.
- Inside Mata, use `rseed()` synchronised to the Stata seed at the start of
  each MCMC routine, and let Mata's `rnormal()` / `runiform()` drive draws.
- All MCMC routines must accept an explicit seed argument so tests are
  deterministic regardless of the caller's global state.

## 5. `e()` returns

Standard set, with spatial additions:

| Return            | Type     | Description                                |
| ----------------- | -------- | ------------------------------------------ |
| `e(b)`            | matrix   | posterior means of β                       |
| `e(V)`            | matrix   | posterior covariance of β                  |
| `e(rho)`          | scalar   | posterior mean of ρ                        |
| `e(rho_se)`       | scalar   | posterior std dev of ρ                     |
| `e(gamma)`        | matrix   | posterior means of convex weights γ        |
| `e(draws_b)`      | matrix   | (optional) raw β draws if `saving()`       |
| `e(draws_rho)`    | matrix   | (optional) raw ρ draws                     |
| `e(draws_gamma)`  | matrix   | (optional) raw γ draws                     |
| `e(impacts)`      | matrix   | direct / indirect / total effects          |
| `e(lmarg)`        | scalar   | log-marginal likelihood                    |
| `e(model_probs)`  | matrix   | posterior model probabilities (BMA only)   |
| `e(N)` `e(T)`     | scalar   | panel dimensions                           |
| `e(cmd)`          | string   | `"spmixw"`                                 |
| `e(model)`        | string   | `"sar"`, `"sdm"`, ...                      |
| `e(ndraw)`        | scalar   | total MCMC draws                           |
| `e(nomit)`        | scalar   | burn-in                                    |
| `e(acceptrate)`   | scalar   | MH acceptance rate for γ (convex models)   |

Large draw matrices are saved only when `saving(filename)` is given;
otherwise `e()` carries summaries only, to avoid bloating `estimates save`.

## 6. Post-estimation

- `estat impact` — direct / indirect / total effects, posterior summaries.
  Mirrors `xsmle`'s API so users don't have to learn a new idiom.
- `estat mcmcdiag` — Geweke z-scores, effective sample size, acceptance
  rates. Implemented in Mata; uses `coda` algorithms re-derived (no
  dependency on R).
- `predict` — leave for a later release. Spatial prediction in panels is
  non-trivial and the R port doesn't expose it cleanly either.

## 6.1 Display conventions (locked 2026-05-08)

Output formatting follows a **hybrid** of the LeSage MATLAB toolbox and the R
port, deliberately diverging from both where it best serves Stata users:

- **ρ lives in the coefficient table as a labelled row** (MATLAB convention).
  Label is `Wy` for SAR/SDM (matches `xsmle`); `We` for SEM/SDEM when added
  in Phase 2.B; γ_k for convex models in Phase 3.
- **Column headers use Bayesian semantics** (closer to R port): `Post. Mean`,
  `Post. SD`, `[level% Cred. Interval]`. The `z` and `P>|z|` columns are
  dropped — they imply a frequentist sampling distribution we don't have.
- **Credible intervals are quantile-based** — read directly from the chain
  via `_spmixw_summary()`, not built from a normal approximation.
- **Scalar block** below the table: `sigma^2 (post. mean)` and elapsed time
  for v0.1; `R-squared`, `corr-squared`, `log-likelihood`, `min/max rho` to
  follow when the underlying quantities are computed (additive change).

Reasoning: putting ρ in the table makes downstream tooling (`estout`,
`esttab`, `coefplot`) work without users learning a special-case extraction;
quantile-based credible intervals are honest under skewed posteriors;
Bayesian column labels prevent the t-stat/p-value misreading.

## 7. Plotting

Stata-native graphics (`twoway kdensity`, `tsline` for trace plots). Avoid
external dependencies. Provide thin wrappers:

- `spmixw_plot, trace`     — trace plots for ρ, γ, β
- `spmixw_plot, post(rho)` — posterior density of ρ
- `spmixw_plot, post(gamma)` — posterior densities of γ_k

## 8. Build & distribution

- `build/build_mlib.do` — compiles `mata/*.mata` into `lspmixw.mlib`
- `build/package.do` — assembles `spmixw.pkg` and `stata.toc` for
  `net install`-able distribution
- Distribution: GitHub release first (`net install spmixw, from(...)`), SSC
  submission once Phase 4 is stable

## 9. Minimum Stata version

**Development environment:** Stata 13 (what the author has installed).
**Declared minimum (v0.1):** Stata 13.0 — the author must be able to run
their own package, and `version 14.0` source directives won't execute on 13.
**Public release minimum (Phase 4 / SSC):** to be lifted to Stata 14 or 15
before submission, since the wider community has moved on.
**Modern path (v0.2+, additive):** Stata 17 features behind a capability check.

### Rationale

- Stata 13 (2013) can compile and run everything we need: Mata `cholesky()`,
  `solvelower()` / `solveupper()`, `panelsetup()`, `rseed()` / `rnormal()`,
  `lmbuild` for `.mlib`, `cscript` for tests. No code-blocking limitation.
- Stata 14 (2015) adds Unicode and minor polish that make `.sthlp` authoring
  cleaner. No code rewrites required between 13 and 14, so the cost of moving
  the declared minimum from 13 → 14 is zero.
- Stata 17+ features (`frames`, faster Mata, raised `matsize`) are real
  improvements. We don't need them for v0.1, but we will *want* them once N·T
  gets large or users want to juggle several W matrices in memory at once.

### Compatibility strategy: dual-path code

All ado entry points start with a capability probe:

```stata
program _spmixw_caps, rclass
    return scalar has_frames = 0
    capture which frames
    if (_rc == 0) return scalar has_frames = 1
end
```

…and any code path that *could* exploit a modern feature degrades gracefully:

| Concern               | Stata 13/14 path (v0.1)                         | Stata 17+ path (v0.2+, optional)                |
| --------------------- | ----------------------------------------------- | ----------------------------------------------- |
| Multi-W input         | Space-separated matrix names; copy into Mata    | `frames` — keep each W in its own frame         |
| Large W storage       | Hand-rolled sparse triplets in Mata             | Same triplets; can additionally back with frame |
| `matsize`             | Document the ~11000 cap; warn in ado            | Skipped (Stata 16+ has no practical cap)        |
| MCMC throughput       | Plain Mata loops                                | Same Mata, but newer Mata is intrinsically faster |
| Help-file Unicode     | ASCII only                                      | UTF-8 allowed                                   |

The Mata numerical core never branches on version. Only the ado-layer
plumbing (input parsing, storage of multiple W matrices) has dual paths, and
the modern path is *additive* — it's a fast-path optimisation, not a
re-implementation.

### Version directives

Every `.ado` and `.mata` file begins with:

```stata
version 13.0
```

That makes the source runnable on the author's Stata 13 *and* on every
later release (newer Statas accept lower `version` directives). When the
declared minimum is lifted at SSC time, this single directive moves to
`version 14.0` (or higher) in one sweep.

## 10. Testing

`tests/` holds `cscript` do-files. One per Mata function and one per ado.
Validation against the R port lives in `validation/` (see PLAN.md §
"Validation strategy").
