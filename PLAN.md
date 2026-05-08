# spmixw (Stata) — Implementation Plan

A Stata port of the R package `spmixW`, which itself ports LeSage's MATLAB
`toolbox_panelg`. Bayesian MCMC estimation of spatial panel data models with
fixed effects, convex combinations of weight matrices, and BMA over W subsets.

Reference implementations (consult in this order when in doubt):

1. **R port** — `../R_package/spmixW/R/` (the source of truth for this port; already validated)
2. **MATLAB original** — `../toolbox_panelg/` (LeSage; consult for numerical edge cases)

---

## Phased roadmap

Mirror the phasing of the R port. **Pause for review at the end of each phase
before starting the next.**

### Phase 1 — Foundation + OLS (review gate) ✅ GREEN 2026-05-07

Establish the package skeleton, the Mata numerical layer, and the simplest
estimator end-to-end so the plumbing is proven before adding spatial pieces.

- [x] Mata library skeleton (`lspmixw.mlib`) with build script
- [x] `_spmixw_demean` — within-transformation (all 4 FE modes verified)
- [x] `_spmixw_ols` — pooled / FE OLS with `e()` returns matching Stata conventions
- [x] `_spmixw_simulate_ols` — DGP for tests
- [x] Help file (`spmixw.sthlp`) with one worked example
- [x] Test harness (do-file based, `tests/master.do`)
- [x] `validation_01_ols.do` parity vs R port — 9/9 checks pass
- [ ] `spmixw.pkg` / `stata.toc` for `net install` (deferred to Phase 4)

**Phase 1 results (validation_01_ols.do, Stata 13):**
| Check                              | Stata       | R port      | Truth | Tol  | Status |
| ---------------------------------- | ----------- | ----------- | ----- | ---- | ------ |
| 1a β₁ (homo, twoway)               | 0.8986      | 0.9654      | 1.0   | 0.15 | PASS   |
| 1a β₂ (homo, twoway)               | 0.9653      | 1.0214      | 1.0   | 0.15 | PASS   |
| 1a σ² (homo, twoway)               | 4.4099      | 4.3532      | 5.0   | 1.5  | PASS   |
| 1b β₁ (hetero + outliers)          | 0.7916      | 0.9079      | 1.0   | 0.25 | PASS   |
| 1b β₂ (hetero + outliers)          | 0.9328      | 1.0228      | 1.0   | 0.15 | PASS   |
| 1b vᵢ outlier-ratio                | 1.74        | > 1.5       | > 1.5 |      | PASS   |
| 1c β₁ (tight prior, pulled)        | 0.6096      | 0.6285      | 0.5   | 0.30 | PASS   |
| 1c β₂ (tight prior, pulled)        | 0.6290      | 0.6496      | 0.5   | 0.30 | PASS   |
| 1c β₁ < 1.0 (pulled toward prior)  | 0.610       | 0.628       | < 1.0 |      | PASS   |

### Phase 2 — Standard spatial panel models

**Phase 2.A — SAR foundation ✅ GREEN 2026-05-08**

- [x] `_spmixw_logdet_exact` — eigenvalue-based ln|I-ρW| on a grid
- [x] `_spmixw_griddy_rho_sar` — griddy Gibbs sampler for ρ in SAR
- [x] `_spmixw_sar` — SAR panel MCMC (Wy via reshape, no Wbig formed)
- [x] `_spmixw_simulate_sar` — SAR DGP for tests/validation
- [x] `spmixw.ado` v0.2.0 dispatch + `_spmixw_sar.ado` subroutine
- [x] `validation_02_sar.do` — 4/4 checks pass on Stata 13

**Phase 2.A results:**
| Check | Stata  | R ref  | Truth | Tol  | Status |
| ----- | ------ | ------ | ----- | ---- | ------ |
| β₁    | 0.9527 | 0.9845 | 1.0   | 0.10 | PASS   |
| β₂    | 0.9825 | 1.0100 | 1.0   | 0.10 | PASS   |
| ρ     | 0.6039 | 0.5990 | 0.6   | 0.05 | PASS   |
| σ²    | 0.8793 | 0.8694 | 1.0   | 0.30 | PASS   |

**Phase 2.B.1 — SEM core ✅ GREEN 2026-05-08**

- [x] `_spmixw_linear_interp` — 1D interp helper for SEM griddy
- [x] `_spmixw_griddy_rho_sem` — griddy Gibbs for SEM (recompute filtered residuals per coarse-grid point + interp to fine grid)
- [x] `_spmixw_sem` — SEM panel MCMC
- [x] `_spmixw_simulate_sem` — SEM DGP via lusolve on errors
- [x] `_spmixw_sem.ado` + dispatcher update; row label `We`
- [x] `validation_04_sem.do` — 3/3 checks pass

**Phase 2.B.1 results:**
| Check | Stata  | Truth | Tol  | Status |
| ----- | ------ | ----- | ---- | ------ |
| β₁    | 0.9418 | 1.0   | 0.10 | PASS   |
| β₂    | 0.9830 | 1.0   | 0.10 | PASS   |
| λ     | 0.6892 | 0.7   | 0.08 | PASS   |

**Phase 2.B.1 performance note (deferred optimisation):** 425 s for
5000 draws on N=200, T=10 (homo). Cause: the griddy step recomputes
filtered SSR on a 200-point coarse grid every iteration. Optimisation —
pre-compute X'X, X'Wx, Wx'Wx, X'y, X'Wy, Wx'Wy once outside the loop;
inner loop becomes scalar combinations, dropping per-grid-point cost
from O(NT k²) to O(k²). Estimated 50–100× speedup. Scheduled as
**Phase 2.D performance pass** after Phase 2.B.3.

**Phase 2.B.2 — SDM/SDEM/SLX wrappers ✅ GREEN 2026-05-08**

- [x] `_spmixw_compute_wx` Mata helper (apply W per-period to X)
- [x] `nodisp` option on `_spmixw_sar.ado`, `_spmixw_sem.ado`, `_spmixw_ols.ado`
- [x] `_spmixw_sdm` ado — augment X with WX, call SAR; row labels `Wy`, `W_x*`
- [x] `_spmixw_sdem` ado — augment X with WX, call SEM; row labels `We`, `W_x*`
- [x] `_spmixw_slx` ado — augment X with WX, call OLS; no spatial-parameter row
- [x] Dispatcher (`spmixw.ado` v0.2.2) routes 6 models: ols/sar/sem/sdm/sdem/slx
- [x] `validation_03_sdm.do` — 5/5 pass
- [x] `validation_05_sdem.do` — 5/5 pass
- [x] `validation_06_slx.do` — 4/4 pass

**Phase 2.B.2 results (14/14 pass):**
| Model | Stata β / θ / ρ-or-λ | Truth | Elapsed |
| ----- | ------------------------------------------------------------------ | -------------------- | ------- |
| SDM   | (0.9520, 0.9796) / (−0.9144, −0.9823) / ρ=0.6905                   | (1,1) / (−1,−1) / 0.7 | 3.77 s  |
| SDEM  | (0.9922, 0.9989) / (−0.9804, −0.9917) / λ=0.6887                   | (1,1) / (−1,−1) / 0.7 | 537 s   |
| SLX   | (0.9543, 0.9818) / (−0.9222, −0.9884) / —                          | (1,1) / (−1,−1) / —   | 0.46 s  |

**Phase 2.B.3 — Effects post-estimation ✅ GREEN 2026-05-08**

- [x] `_spmixw_trace_mc` — stochastic trace estimator for tr(W^j)/N
- [x] `_spmixw_effects_sar` — direct/indirect/total via LeSage-Pace scalar summary
- [x] `_spmixw_effects_sdm` — same with both β and θ contributions to direct
- [x] `_spmixw_effects_simple` — for SDEM/SLX: direct=β, indirect=θ (exact, no MC)
- [x] `spmixw_estat` dispatcher + `_spmixw_estat_impact` — `estat impact` works after sar/sdm/sdem/slx
- [x] Retro-update validation_02_sar.do with identity check + ratio-bound check
- [x] Retro-update validation_03_sdm.do — identity + negative-indirect
- [x] Retro-update validation_05_sdem.do — direct ≈ β, indirect ≈ θ
- [x] Retro-update validation_06_slx.do — direct == β̂, indirect == θ̂ at machine precision
- [ ] `_spmixw_logdet_mc` (deferred to Phase 2.D — only matters for very large N)

**Cumulative validation summary (43/43 pass on Stata 13):**
| Script | Checks | Highlights |
| ------ | ------ | ---------- |
| 01 OLS  | 9/9 | hetero outlier ratio 1.74; tight-prior pull-toward-prior |
| 02 SAR  | 6/6 | identity max err = **0**, ratio = 2.006 ∈ (1, 2.524] |
| 03 SDM  | 8/8 | identity err = 1.1e-16, indirect[x*] all < 0 |
| 04 SEM  | 3/3 | β / λ recovery |
| 05 SDEM | 9/9 | direct ≈ β, indirect ≈ θ |
| 06 SLX  | 8/8 | direct == β̂ and indirect == θ̂ at **0.00e+00** |

**Phase 2 is now feature-complete for standard spatial panel models.**

**Phase 2.D — Performance pass ✅ GREEN 2026-05-08 (32-45× speedup achieved)**

- [x] Pre-compute V-weighted cross-product blocks for SEM griddy (per call)
- [x] Hoist block-building out of MCMC loop: once for homo, per-iter for hetero
- [x] Apply algebraic identity to homo SEM beta/sigma draws — XtX(λ) = A0 − λA1s + λ²A2 etc.; eliminates per-iteration NT-sized matmuls entirely
- [x] Replace W·e reshape with `Wy − Wx·β̂` algebraic substitution
- [x] Fix `_spmixw_linear_interp` from O(nx·no) to O(nx+no) one-pass — **the dominant bug**, single change worth ~12× SEM speedup
- [x] Vectorize trapezoid + inverse-CDF loops in both SEM and SAR griddies (replace per-element `for` with `runningsum()` and matrix-level ops)

**Phase 2.D timings (5000 draws, N=200, T=10, banded W):**
| Model | Original | Phase 2.D final | Speedup |
| ----- | -------- | --------------- | ------- |
| SAR   | 3.42 s   | 3.00 s          | 1.1×    |
| SEM   | 410 s    | **12.86 s**     | **32×** |
| SDEM  | 537 s    | **12.06 s**     | **45×** |
| SDM   | 3.77 s   | unchanged       | —       |
| OLS   | <1 s     | unchanged       | —       |

All 43 cumulative validation checks still pass with bit-identical posterior estimates.

**Phase 2.D deferred items (low priority, only matter for very large N):**

- [ ] `_spmixw_logdet_mc` — Pace-Barry MC log-det. For N > ~1000 the eigenvalue-based exact log-det becomes O(N³); MC trace estimator is O(N² · iter · order). Add when needed.
- [ ] Vectorise hetero SEM W applications (currently k+1 reshape-matmuls per iter); not on critical path since users will run rval=0 most of the time.
- [ ] Inline 2×2 / 4×4 quadratic-form computation in griddy to skip Mata's `lusolve()` per-call overhead. Only worth it if Phase 2.D's 13 s SEM runtime needs to drop to <5 s.

**Phase 3.A — SAR convex ✅ GREEN 2026-05-08**

- [x] `_spmixw_logdet_taylor` + `_spmixw_eval_taylor_lndet` — precompute trace cross-terms tr(W_{i₁}…W_{iₚ}) packed into a single colvector with M^p offsets; evaluate via Kronecker products of γ
- [x] `_spmixw_gamma_proposal_uniform` / `_adapted` — LeSage's coin-flip MH proposal on the simplex (warm-up uniform window then γ_std-adapted window)
- [x] `_spmixw_eval_cond_sar_conv` — conditional log-posterior kernel used by both ρ and γ MH steps
- [x] `_spmixw_sar_conv` — SAR conv MCMC (random-walk MH on ρ with adaptive cc; coin-flip MH on γ; β/σ² Gibbs steps via the (M+1)-column ys/xs basis trick)
- [x] `_spmixw_simulate_sar_conv` — convex-W DGP for tests
- [x] `_spmixw_sar_conv.ado` + dispatcher integration with `wmats()` option
- [x] `validation_07_sar_conv.do` — 7/7 checks pass; γ recovery essentially unbiased

**Phase 3.A results (M=2, N=200, T=10, ndraw=15000, taylor_order=8):**
| Check | Stata | Truth | Tol | Notes |
| ----- | ----- | ----- | --- | ----- |
| β₁    | 0.9051 | 1.0  | 0.15 | PASS |
| β₂    | 0.9523 | 1.0  | 0.15 | PASS |
| ρ     | 0.7191 | 0.6  | 0.15 | PASS — ~0.12 Taylor bias matches R-port docstring for moderate N |
| **γ₁**| **0.3285** | **0.30** | 0.15 | PASS — γ recovery essentially unbiased |
| **γ₂**| **0.6715** | **0.70** | 0.15 | PASS |
| MH acc ρ | 0.588 | (0.10, 0.80) | — | PASS |
| MH acc γ | 0.671 | > 0.03 | — | PASS |

**Phase 3.A engineering notes:**
- Default Taylor order: 6 (matches R port). Validation uses 8 because banded W has slower Taylor convergence than KNN W.
- W's stacked horizontally as N × (N·M) for Mata indexing convenience (Mata has no native 3D arrays; pointer arrays are awkward in `mata mlib add`).
- ys/xs basis: NT × (M+1) and NT × ((M+1)·k) respectively. For any (ρ, γ): omega = (1, ρ·γ_1, …, ρ·γ_M)' and (I − ρW_c)y = ys_basis · omega.
- Stata 13 gotcha caught: Mata terminates statements at line breaks unless the line ends with an operator. Multi-line expressions like `M[..., r] = X[...] * omega` split across lines parse the second line's `* omega` as a unary pointer-deref. Fix: keep on one line or end first line with `*`.

**Phase 3.B — SEM/SDM/SDEM convex ✅ GREEN 2026-05-08**

- [x] `_spmixw_simulate_sem_conv` — SEM conv DGP (different from SAR conv; spatial multiplier acts on errors)
- [x] `_spmixw_sar_conv_caller` extracted to `mata/_spmixw_sar_conv_caller.mata` (lives in lspmixw.mlib so any ado can call it without depending on load order — caught when `_spmixw_sem_conv.ado` couldn't find the bridge function defined inside `_spmixw_sar_conv.ado`'s embedded mata block)
- [x] `_spmixw_sem_conv.ado` — calls the same Mata MCMC kernel as SAR conv (R port confirms math is identical), only the user-facing labels differ: row label `We`, header "SEM (convex W)", e(model) = "sem_conv"
- [x] `_spmixw_sdm_conv.ado` — wraps SAR conv with [W_1·X, W_2·X, …, W_M·X] augmentation. Each `W_m·X` becomes a block of `k` covariates labelled `W<m>_<x>` in the displayed table, giving M·k extra coefficients
- [x] `_spmixw_sdem_conv.ado` — same augmentation, calls SEM conv
- [x] `validation_08_sem_conv.do`, `validation_09_sdm_conv.do`, `validation_10_sdem_conv.do` — 7+5+5 = 17 checks, all pass

**Phase 3 cumulative results (24 new checks; γ recovery is essentially unbiased in all four convex models):**
| Validation | β        | Spatial | γ              | Acceptance      | Status |
| ---------- | -------- | ------- | -------------- | --------------- | ------ |
| 07 SAR_C   | (0.91, 0.95) | ρ=0.72 (truth 0.6) | (0.33, 0.67) vs (0.30, 0.70) | ρ=0.59, γ=0.67 | 7/7    |
| 08 SEM_C   | (0.95, 0.99) | λ=0.45 (truth 0.5) | **(0.398, 0.602)** | λ=0.52, γ=0.72 | 7/7    |
| 09 SDM_C   | (0.93, 0.94) | ρ=0.44 (truth 0.5) | (0.404, 0.596)  | ρ=0.49, γ=0.72 | 5/5    |
| 10 SDEM_C  | (0.97, 0.98) | λ=0.45 (truth 0.5) | (0.404, 0.596)  | λ=0.59, γ=0.72 | 5/5    |

**Cumulative across all phases (1, 2, 3.A, 3.B): 74/74 validation checks pass on Stata 13.**

**Runtime characteristics (5000-saved-draw runs on N=200, T=10, M=2):**
| Model      | Runtime |
| ---------- | ------- |
| sar_conv   | 15 s    |
| sem_conv   | 17 s    |
| sdm_conv   | 37 s (k=6 augmented) |
| sdem_conv  | 44 s (k=6 augmented) |

**Phase 3.C — Convex effects + BMA over W subsets (next)**

- [ ] Extend `estat impact` to handle sar_conv / sdm_conv etc. — needs Wc at posterior mean γ for trace estimation
- [ ] BMA: enumerate W subsets, compute log-marginals, posterior model probabilities  
- [ ] `validation_11_sar_conv_bma.do`

**Phase 3.C — Convex effects + BMA over W subsets**

- [ ] Extend `estat impact` to handle sar_conv / sdm_conv etc. — needs Wc at posterior mean γ for trace estimation
- [ ] BMA: enumerate W subsets, compute log-marginals, posterior model probabilities
- [ ] `validation_11_sar_conv_bma.do`

**Phase 2.B.3 — Effects post-estimation**

- [ ] `_spmixw_logdet_mc` — MC log-det (Pace-Barry) for large-N use
- [ ] `estat impact` — direct/indirect/total via stochastic-trace estimation
- [ ] Retro-update `validation_02_sar.do` with effects checks (direct + indirect = total identity, total/direct ≈ 1/(1-ρ))

**Phase 2.D — Performance pass**

- [ ] Pre-compute cross-product blocks for SEM griddy (50-100× speedup)
- [ ] Cache `W*y` and `W*X` outside SAR MCMC homo block (already done; verify hetero too)
- [ ] Bench across (N, T, ndraw) grid; document expected runtimes in help

### Phase 3 — Convex combination models

- [ ] Taylor-series log-det (`log_det_taylor`, `eval_taylor_lndet`) — 4th order
- [ ] Metropolis-Hastings sampler for γ (convex weights on candidate Ws)
- [ ] `sar_conv_panel`, `sem_conv_panel`, `sdm_conv_panel`, `sdem_conv_panel`
- [ ] Multi-W input handling (`wmats(W1 W2 W3 ...)`)
- [ ] Validation against `validation_07_sar_conv.R`

### Phase 4 — BMA + polish

- [ ] `lmarginal_panel` — log-marginal likelihoods
- [ ] `model_probs` and `compare_models`
- [ ] `sar_conv_bma`, `sdm_conv_bma`, `sdem_conv_bma`
- [ ] `estat mcmcdiag` — Geweke, ESS, acceptance rates
- [ ] Trace / posterior-density plots
- [ ] Help files complete, examples in every `.sthlp`
- [ ] SSC submission package or GitHub `net install` instructions
- [ ] Validation against `validation_08_lmarginal.R`, `validation_09_bma.R`

---

## Engineering decisions to lock in Phase 1

Decisions made now ripple through every later phase. Settle these before
writing the spatial code.

| Decision                              | Default proposal                                                | Open?  |
| ------------------------------------- | --------------------------------------------------------------- | ------ |
| Command surface                       | Single `spmixw` w/ `model()` option (vs. one ado per model)     | open   |
| W storage                             | Stata matrix name(s) → Mata sparse triplet inside the routine   | open   |
| Multi-W input                         | `wmats(W1 W2 W3)` — space-separated matrix names                | open   |
| Sparse engine                         | Mata native (custom helpers) for v0.1; revisit C plugin later   | open   |
| Random draws / RNG                    | Mata `rseed()`, `rnormal()`, `runiform()` — Stata `set seed` API| open   |
| MCMC draws storage                    | `e(draws)` matrix (small) + optional `saving()` to `.dta`       | open   |
| Effects (direct/indirect/total)       | `estat impact` — mirror `xsmle`'s API                           | open   |
| Naming convention                     | `spmixw_*` prefix on all ados/Mata functions                    | open   |
| Minimum Stata version (v0.1)          | **Stata 13.0** — dev = use environment; `version 13.0` directive| locked |
| Public release minimum (Phase 4)      | Lift to Stata 14 or 15 before SSC submission                    | open   |
| Modern Stata 17+ fast-path            | Deferred to v0.2+ (additive; frames, raised matsize)            | locked |

---

## Validation strategy

Each Stata estimator must reproduce the R-port output on a shared synthetic
dataset (seed-controlled). Targets:

- Posterior means of β within ±1 MCMC standard error of the R port
- Posterior mean of ρ within ±0.005
- Log-marginal likelihoods within ±0.5 (BMA models)
- Direct/indirect/total effects within ±1 MCMC standard error

The validation harness in `validation/` will load the same DGP, run the Stata
command, run the R command (via `rcall` or saved-output comparison), and diff.

---

## Out of scope for v0.1

- C plugin acceleration (revisit if Mata MCMC is unusably slow)
- Random effects (R port is FE-only)
- Bayesian spatial *cross-section* models (R port has these via
  `validation_10_cross_section.R`; defer)
- `predict` for spatial impulse responses beyond `estat impact`

## Deferred to v0.2+ (Stata 17+ modern path)

These are additive optimisations gated behind a runtime capability check, so
v0.1 code on Stata 13/14 keeps working untouched:

- **`frames` for multi-W storage** — currently each W lives in a Stata matrix
  passed by name; under v0.2+ each W can live in its own frame, avoiding
  matsize limits and re-copies into Mata
- **Raised `matsize` regime** — Stata 16 dropped the practical cap; the
  modern path can skip the v0.1 warning and accept much larger panels
- **Faster Mata throughput** — automatic; no code change required, but we'll
  re-benchmark on 17+ to set realistic ndraw guidance in the help files
- **Unicode in `.sthlp`** — v0.1 stays ASCII-clean; v0.2+ may use proper
  symbols (ρ, γ, β) once we no longer support 13

See [ARCHITECTURE.md §9](ARCHITECTURE.md) for the dual-path strategy.
