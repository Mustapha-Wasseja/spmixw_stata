*! Validation 6: SLX Panel — Stata port of validation_06_slx.R
*!
*! DGP: y = X*beta + WX*theta + FE + eps  (no spatial inverse)
*! beta = (1, 1),  theta = (-1, -1),  sige = 1
*! For SLX: direct = beta and indirect = theta exactly (no spatial multiplier),
*! so the parameter checks below double as effects checks.
version 13.0

if "${SPMIXW_PKG}" == "" global SPMIXW_PKG ".."
quietly do "${SPMIXW_PKG}/tests/_setup.do"

clear all
quietly do "${SPMIXW_PKG}/tests/_setup.do"
quietly do "${SPMIXW_PKG}/validation/validation_helpers.do"

set seed 10203040

local N    = 200
local T    = 10

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

mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    W    = st_matrix("W")
    nobs = N * TT
    beta_true  = (1 \ 1)
    theta_true = (-1 \ -1)
    sigma2     = 1

    rseed(10203040)
    X  = rnormal(nobs, 2, 0, 1)
    WX = _spmixw_compute_wx(X, W, N, TT)

    sfe = (1::N)  :/ N
    tfe = (1::TT) :/ TT
    fe  = (J(TT,1,1) # sfe) + (tfe # J(N,1,1))

    eps = rnormal(nobs, 1, 0, sqrt(sigma2))
    y   = X * beta_true + WX * theta_true + fe + eps

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
print_header Validation 6: SLX Panel - Homoscedastic, Two-way FE

spmixw y x1 x2, model(slx) w(W) effects(twoway) rval(0) ///
    ndraw(5000) nomit(1500) seed(10203040)

check_param, name("beta_1")  estimate(`=_b[x1]')   truth(1.0)  tol(0.10)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2")  estimate(`=_b[x2]')   truth(1.0)  tol(0.10)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("theta_1") estimate(`=_b[W_x1]') truth(-1.0) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("theta_2") estimate(`=_b[W_x2]') truth(-1.0) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// =============================================================
// Effects checks (Phase 2.B.3): for SLX, direct = beta and
// indirect = theta EXACTLY (no spatial multiplier on X), so the
// effects table reproduces the coefficient table to machine precision.
// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Effects checks via estat impact (SLX)"
di as txt "{hline 70}"

estat impact

matrix dir_t = r(direct)
matrix ind_t = r(indirect)

local b1 = _b[x1]
local b2 = _b[x2]
local t1 = _b[W_x1]
local t2 = _b[W_x2]

local d1 = dir_t[1, 1]
local d2 = dir_t[2, 1]
local i1 = ind_t[1, 1]
local i2 = ind_t[2, 1]

// Exact equality (within tiny float tolerance) since SLX has no MC step
// for effects.
local pass_d1 = (abs(`d1' - `b1') < 1e-10)
local pass_d2 = (abs(`d2' - `b2') < 1e-10)
local pass_i1 = (abs(`i1' - `t1') < 1e-10)
local pass_i2 = (abs(`i2' - `t2') < 1e-10)

di as txt "  direct[x1]   == β̂_1   : err = `: di %9.2e abs(`d1' - `b1')'  " ///
   as `=cond(`pass_d1', "txt", "err")' cond(`pass_d1', "PASS", "FAIL")
di as txt "  indirect[x1] == θ̂_1   : err = `: di %9.2e abs(`i1' - `t1')'  " ///
   as `=cond(`pass_i1', "txt", "err")' cond(`pass_i1', "PASS", "FAIL")
di as txt "  direct[x2]   == β̂_2   : err = `: di %9.2e abs(`d2' - `b2')'  " ///
   as `=cond(`pass_d2', "txt", "err")' cond(`pass_d2', "PASS", "FAIL")
di as txt "  indirect[x2] == θ̂_2   : err = `: di %9.2e abs(`i2' - `t2')'  " ///
   as `=cond(`pass_i2', "txt", "err")' cond(`pass_i2', "PASS", "FAIL")

local NPASS = `NPASS' + `pass_d1' + `pass_d2' + `pass_i1' + `pass_i2'
local NFAIL = `NFAIL' + (4 - `pass_d1' - `pass_d2' - `pass_i1' - `pass_i2')

// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Validation 6 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
