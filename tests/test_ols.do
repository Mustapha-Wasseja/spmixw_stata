*! Test: spmixw, model(ols) recovers known betas on synthetic data
version 13.0

if "${SPMIXW_PKG}" == "" global SPMIXW_PKG ".."

clear all
quietly do "${SPMIXW_PKG}/tests/_setup.do"

di as txt _newline "----- test_ols: spmixw model(ols) recovers known betas -----"

set seed 20260507

// -- Synthetic panel via the in-package DGP -------------------------------
local N    = 100
local T    = 8
local k    = 2

mata:
    N  = strtoreal(st_local("N"))
    TT = strtoreal(st_local("T"))
    beta_true = (1.5 \ -0.8)
    sigma2    = 0.25
    y = .; X = .
    _spmixw_simulate_ols(N, TT, beta_true, sigma2, "twoway", y, X)

    // Build a Stata dataset, sorted by time then region
    region = J(TT, 1, 1) # (1::N)
    year   = (1::TT) # J(N, 1, 1)
    M      = (region, year, y, X)
    st_addvar("int", "region")
    st_addvar("int", "year")
    st_addvar("double", "y")
    for (j = 1; j <= cols(X); j++) {
        st_addvar("double", "x" + strofreal(j))
    }
    st_addobs(rows(M))
    st_store(., ., M)
end

xtset region year

// -- Run the estimator ----------------------------------------------------
spmixw y x1 x2, model(ols) effects(twoway) rval(0) ///
    ndraw(2000) nomit(500) seed(20260507)

// -- Recovery checks ------------------------------------------------------
local b1 = _b[x1]
local b2 = _b[x2]
local s2 = e(sige)

di as txt "  beta_1 (true 1.5) = " as res %8.4f `b1'
di as txt "  beta_2 (true -0.8) = " as res %8.4f `b2'
di as txt "  sigma^2 (true 0.25) = " as res %8.4f `s2'

// Tolerances: with N=100, T=8 the standard error on beta is small; 0.1 is
// generous slack to absorb MCMC variability across machines.
assert abs(`b1' - 1.5) < 0.10
assert abs(`b2' - (-0.8)) < 0.10
assert abs(`s2' - 0.25) < 0.15

di as txt "ols recovery test passed"
