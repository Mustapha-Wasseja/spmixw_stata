# spmixw (Stata)

**Bayesian spatial panel data models with convex combinations of weight
matrices.**

Stata port of the R package `spmixW`, which itself ports James P. LeSage's
MATLAB `Panel Data Toolbox` (`toolbox_panelg`). MCMC estimation of SAR, SEM,
SDM, SDEM, and SLX panel models with fixed effects, plus convex combinations
of multiple W matrices (Debarsy & LeSage 2021) and Bayesian Model Averaging
over W subsets.

> **Status:** v0.1.0 — feature-complete, 120/120 validation checks passing.
> Pre-SSC; install via local copy or `net install` from a release archive.

## Models supported

| `model()`     | Description                                              |
| ------------- | -------------------------------------------------------- |
| `ols`         | Bayesian fixed-effects OLS (no spatial component)        |
| `sar`         | Spatial autoregressive: y = ρ W y + Xβ + ε               |
| `sem`         | Spatial error: y = Xβ + u, u = λ W u + ε                 |
| `sdm`         | Spatial Durbin: y = ρ W y + Xβ + W X θ + ε               |
| `sdem`        | Spatial Durbin error: y = Xβ + W X θ + u, u = λ W u + ε  |
| `slx`         | Spatial lag of X: y = Xβ + W X θ + ε                     |
| `sar_conv`    | SAR with W_c(γ) = Σ_m γ_m W_m  (convex combination)      |
| `sem_conv`    | SEM convex combination                                   |
| `sdm_conv`    | SDM convex combination                                   |
| `sdem_conv`   | SDEM convex combination                                  |

Plus `spmixw_bma` for Bayesian Model Averaging over the 2^M − 1 non-empty
subsets of M weight matrices, and `estat impact` for LeSage-Pace direct /
indirect / total effects.

## Installation

### From a local copy (development)

```stata
* Add the package's ado/ directory to Stata's adopath
adopath ++ "/path/to/spmixw_stata/ado"
```

### From a release archive (planned for v0.1.0)

```stata
net install spmixw, from("https://example.org/spmixw")
```

### From SSC (planned)

```stata
ssc install spmixw
```

The Mata library `lspmixw.mlib` ships pre-built in `ado/`. To rebuild from
source:

```stata
do "spmixw_stata/build/build_mlib.do"
```

## Quick start

```stata
* Stata 13.0 or later
xtset region year                       // balanced panel required

* Single-W SAR
spmixw y x1 x2, model(sar) w(W) effects(twoway) ndraw(5000) nomit(1000)
estat impact                            // direct/indirect/total effects

* Convex combination of three W matrices
spmixw y x1 x2, model(sar_conv) wmats(W_contig W_knn W_distance) ///
    effects(twoway) ndraw(15000) nomit(5000) seed(42)

* BMA over W subsets
spmixw_bma y x1 x2, model(sdm_conv) wmats(W1 W2 W3) effects(twoway) ///
    ndraw(8000) nomit(2000) seed(42)
```

A complete worked example with synthetic regions and three W matrices ships
in `examples/`. Run `do examples/spmixw_demo.do` after adding the package to
your adopath.

## Folder layout

```
spmixw_stata/
├── README.md             # this file
├── LICENSE               # MIT
├── PLAN.md               # phased development roadmap
├── ARCHITECTURE.md       # design notes (command surface, Mata layer, RNG)
├── ado/                  # user-facing .ado commands + lspmixw.mlib
├── mata/                 # Mata source (.mata files, compiled into lspmixw.mlib)
├── help/                 # .sthlp help files
├── tests/                # cscript test suite
├── validation/           # validation against R port + DGP-based checks
├── build/                # build_mlib.do, archival lspmixw.mlib copy
├── examples/             # worked example: synthetic data + demo .do
└── docs/                 # design memos, derivations
```

## Requirements

- **Stata 13.0** or later (no Stata 14+ features required for v0.1)
- Balanced panel data (warning issued if unbalanced)

## Output conventions

The coefficient table uses Bayesian column headers:

```
            y  |  Post. Mean   Post. SD   [95% Cred. Interval]
```

Frequentist columns (`z`, `P>|z|`) are omitted because they imply a sampling
distribution we don't have. Credible intervals are quantile-based (no
normality assumption). The spatial parameter ρ (SAR / SDM) or λ (SEM / SDEM)
appears as a labelled row (`Wy` or `We`) in the same table, matching `xsmle`
and supporting downstream tools (`estout`, `coefplot`, …).

## References

- Debarsy, N. and LeSage, J. P. (2021). "Bayesian model averaging for
  spatial autoregressive models based on convex combinations of different
  types of connectivity matrices." *Journal of Business & Economic
  Statistics* 40(2), 547-558.
- LeSage, J. P. and Pace, R. K. (2009). *Introduction to Spatial
  Econometrics.* Boca Raton, FL: Chapman & Hall/CRC.
- Geweke, J. (1993). "Bayesian treatment of the independent Student-t
  linear model." *Journal of Applied Econometrics* 8(S1), S19-S40.

## Reference implementations

- **R port (source of truth):** `../R_package/spmixW/`
- **MATLAB original:** `../toolbox_panelg/`

When numerical behaviour is ambiguous, the R port wins. The MATLAB original
is consulted only for edge cases the R port itself flags.

## Author

Mustapha Wasseja Mohammed — <muswaseja@gmail.com>

## License

MIT. See [LICENSE](LICENSE).
