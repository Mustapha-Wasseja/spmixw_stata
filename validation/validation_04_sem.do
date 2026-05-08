*! Validation 4: SEM Panel — Stata port of validation_04_sem.R
*!
*! DGP: u = (I - λ*Wbig)^{-1} eps;  y = X*β + FE + u
*! λ = 0.7,  β = (1, 1),  σ² = 1
*! W: banded contiguity (each region has neighbours {i-1, i+1}), row-normalised.
version 13.0

if "${SPMIXW_PKG}" == "" global SPMIXW_PKG ".."
quietly do "${SPMIXW_PKG}/tests/_setup.do"

clear all
quietly do "${SPMIXW_PKG}/tests/_setup.do"
quietly do "${SPMIXW_PKG}/validation/validation_helpers.do"

set seed 10203040

local N    = 200
local T    = 10

// -- Build a deterministic banded W (neighbours {i-1, i+1}, row-normalised) -
mata:
    N = strtoreal(st_local("N"))
    W = J(N, N, 0)
    for (i = 1; i <= N; i++) {
        if (i > 1) W[i, i-1] = 1
        if (i < N) W[i, i+1] = 1
    }
    rs = rowsum(W)
    for (i = 1; i <= N; i++) {
        W[i, .] = W[i, .] :/ rs[i]
    }
    st_matrix("W", W)
end

// -- Generate SEM DGP via the package simulator ---------------------------
mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    lambda_true = 0.7
    beta_true   = (1 \ 1)
    sigma2      = 1
    W = st_matrix("W")
    y = .; X = .
    _spmixw_simulate_sem(N, TT, W, lambda_true, beta_true, sigma2, "twoway", y, X)

    region = J(TT, 1, 1) # (1::N)
    year   = (1::TT) # J(N, 1, 1)
    M = (region, year, y, X)

    st_addvar("int",    "region")
    st_addvar("int",    "year")
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
print_header Validation 4: SEM Panel - Homoscedastic, Two-way FE

spmixw y x1 x2, model(sem) w(W) effects(twoway) rval(0) ///
    ndraw(5000) nomit(1500) seed(10203040)

check_param, name("beta_1") estimate(`=_b[x1]') truth(1.0) tol(0.10)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2") estimate(`=_b[x2]') truth(1.0) tol(0.10)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("lambda") estimate(`=e(rho)') truth(0.7) tol(0.08)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Validation 4 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
