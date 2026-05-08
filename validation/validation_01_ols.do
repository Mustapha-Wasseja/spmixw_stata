*! Validation 1: OLS Panel — Stata port of validation_01_ols.R
*!
*! Mirrors the three runs in the R script:
*!   (a) Homoscedastic, two-way FE
*!   (b) Heteroscedastic + outliers in t=5,6, region FE
*!   (c) Tight informative prior pulling beta toward 0.5, region FE
*!
*! NOTE on cross-port parity: Stata and R use different RNG implementations,
*! so draw-for-draw equality is impossible. We compare *posterior summaries*
*! against the same true-DGP values used in the R script, with tolerances
*! generous enough to absorb both Monte Carlo error and RNG differences.
*! Reference column ("R ref") shows the R-port estimates for visual context.
version 13.0

if "${SPMIXW_PKG}" == "" global SPMIXW_PKG ".."
quietly do "${SPMIXW_PKG}/tests/_setup.do"

clear all
quietly do "${SPMIXW_PKG}/tests/_setup.do"
quietly do "${SPMIXW_PKG}/validation/validation_helpers.do"

set seed 10203040

// -- DGP (matches the R script structurally) ------------------------------
local N    = 200
local T    = 10
local nobs = `N' * `T'

mata:
    N      = strtoreal(st_local("N"))
    TT     = strtoreal(st_local("T"))
    nobs   = N * TT
    bt     = (1 \ 1)
    sigt   = 5

    rseed(10203040)
    X = rnormal(nobs, 2, 0, 1)

    sfe = (1::N)  :/ N
    tfe = (1::TT) :/ TT
    fe  = (J(TT,1,1) # sfe) + (tfe # J(N,1,1))

    eps = rnormal(nobs, 1, 0, sqrt(sigt))
    y   = X * bt + fe + eps

    region = J(TT, 1, 1) # (1::N)
    year   = (1::TT) # J(N, 1, 1)
    M = (region, year, y, X)

    st_addvar("int", "region")
    st_addvar("int", "year")
    st_addvar("double", "y")
    st_addvar("double", "x1")
    st_addvar("double", "x2")
    st_addobs(rows(M))
    st_store(., ., M)
end

xtset region year

local NPASS = 0
local NFAIL = 0

// =============================================================
// (a) Homoscedastic, two-way FE
// =============================================================
print_header Validation 1a: OLS Panel - Homoscedastic, Two-way FE

spmixw y x1 x2, model(ols) effects(twoway) rval(0) ///
    ndraw(2500) nomit(500) seed(10203040)

check_param, name("beta_1") estimate(`=_b[x1]') truth(1.0) tol(0.15) reference("0.9654")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2") estimate(`=_b[x2]') truth(1.0) tol(0.15) reference("1.0214")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("sigma^2") estimate(`=e(sige)') truth(5.0) tol(1.5) reference("4.3532")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// =============================================================
// (b) Heteroscedastic + outliers in t = 5, 6 (region FE)
// =============================================================
print_header Validation 1b: OLS Panel - Heteroscedastic, Outliers in T=5,6

quietly replace y = y + 10 if inlist(year, 5, 6)

spmixw y x1 x2, model(ols) effects(region) rval(5) ///
    ndraw(2500) nomit(500) seed(10203040)

check_param, name("beta_1") estimate(`=_b[x1]') truth(1.0) tol(0.25) reference("0.9079")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2") estimate(`=_b[x2]') truth(1.0) tol(0.15) reference("1.0228")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// vmean on outlier rows must be larger than on non-outlier rows
matrix vmean = e(vmean)
mata:
    vm = st_matrix("vmean")
    NN = strtoreal(st_local("N"))
    TT = strtoreal(st_local("T"))
    // Data is sorted by (year, region) so rows for year y are
    //   ((y-1)*N + 1) .. (y*N)
    idx5 = ((5-1)*NN + 1)::(5*NN)
    idx6 = ((6-1)*NN + 1)::(6*NN)

    // Build an explicit 0/1 mask. Avoids reliance on Mata broadcasting,
    // which Stata 13's Mata refuses on (NT x 1) :== (1 x 2N).
    out_mask = J(rows(vm), 1, 0)
    out_mask[(idx5 \ idx6), 1] = J(rows(idx5) + rows(idx6), 1, 1)

    vi_out   = sum(select(vm, out_mask))         / sum(out_mask)
    vi_other = sum(select(vm, 1 :- out_mask))    / sum(1 :- out_mask)

    st_local("vi_ratio", strofreal(vi_out / vi_other))
end

di as txt "  vi(t=5,6) / vi(other) = " as res `vi_ratio' ///
   as txt " (should be > 1.5)"
local pass_v = real("`vi_ratio'") > 1.5
local NPASS = `NPASS' + `pass_v'
local NFAIL = `NFAIL' + (1 - `pass_v')

// Restore y for run (c)
quietly replace y = y - 10 if inlist(year, 5, 6)

// =============================================================
// (c) Informative prior pulling beta toward 0.5
// =============================================================
print_header Validation 1c: OLS Panel - Tight Prior beta=(0.5,0.5)

matrix b0 = (0.5 \ 0.5)
matrix B0 = I(2) * 0.001

spmixw y x1 x2, model(ols) effects(region) rval(0) ///
    ndraw(2500) nomit(500) seed(10203040) ///
    bprior(b0) bvar(B0)

check_param, name("beta_1 (biased)") estimate(`=_b[x1]') truth(0.5) tol(0.3) reference("0.6285")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2 (biased)") estimate(`=_b[x2]') truth(0.5) tol(0.3) reference("0.6496")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// Beta_1 should be pulled below 1.0 by the tight prior
local pulled = (_b[x1] < 1.0)
di as txt "  beta_1 < 1.0 (pulled toward prior): " ///
   as `=cond(`pulled', "txt", "err")' cond(`pulled', "PASS", "FAIL")
local NPASS = `NPASS' + `pulled'
local NFAIL = `NFAIL' + (1 - `pulled')

// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Validation 1 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
