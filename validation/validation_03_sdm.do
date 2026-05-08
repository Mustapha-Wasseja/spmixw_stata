*! Validation 3: SDM Panel — Stata port of validation_03_sdm.R
*!
*! DGP: y = (I - rho*Wbig)^{-1} ( X*beta + WX*theta + FE + eps )
*! rho = 0.7,  beta = (1, 1),  theta = (-1, -1),  sige = 1
*! Effects identity / negative-indirect checks deferred to Phase 2.B.3.
version 13.0

if "${SPMIXW_PKG}" == "" global SPMIXW_PKG ".."
quietly do "${SPMIXW_PKG}/tests/_setup.do"

clear all
quietly do "${SPMIXW_PKG}/tests/_setup.do"
quietly do "${SPMIXW_PKG}/validation/validation_helpers.do"

set seed 10203040

local N    = 200
local T    = 10

// -- Build banded W -------------------------------------------------------
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

// -- Generate SDM DGP ------------------------------------------------------
mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    W    = st_matrix("W")
    nobs = N * TT
    rho_true   = 0.7
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
    mu  = X * beta_true + WX * theta_true + fe + eps

    // y = (I_N - rho*W)^{-1} mu, applied per period
    A     = I(N) - rho_true :* W
    M_mat = rowshape(mu, TT)'
    Y_mat = lusolve(A, M_mat)
    y     = vec(Y_mat)

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
print_header Validation 3: SDM Panel - Homoscedastic, Two-way FE

spmixw y x1 x2, model(sdm) w(W) effects(twoway) rval(0) ///
    ndraw(5000) nomit(1500) seed(10203040)

check_param, name("beta_1")  estimate(`=_b[x1]')   truth(1.0)  tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2")  estimate(`=_b[x2]')   truth(1.0)  tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("theta_1") estimate(`=_b[W_x1]') truth(-1.0) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("theta_2") estimate(`=_b[W_x2]') truth(-1.0) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("rho")     estimate(`=e(rho)')   truth(0.7)  tol(0.08)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// =============================================================
// Effects checks (Phase 2.B.3): identity holds, indirect should be
// negative because theta = -1 dominates the spatial-multiplier effect.
// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Effects checks via estat impact (SDM)"
di as txt "{hline 70}"

estat impact, seed(20260508)

matrix dir_t = r(direct)
matrix ind_t = r(indirect)
matrix tot_t = r(total)

mata:
    dir_pm = st_matrix("dir_t")[., 1]
    ind_pm = st_matrix("ind_t")[., 1]
    tot_pm = st_matrix("tot_t")[., 1]
    max_id = max(abs(dir_pm + ind_pm - tot_pm))
    st_local("identity_err",   strofreal(max_id))
    st_local("indirect_x1",    strofreal(ind_pm[1]))
    st_local("indirect_x2",    strofreal(ind_pm[2]))
end

local pass_id = (real("`identity_err'") < 1e-8)
di as txt "  direct + indirect = total identity: max err = `identity_err'  " ///
   as `=cond(`pass_id', "txt", "err")' cond(`pass_id', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_id'
local NFAIL = `NFAIL' + (1 - `pass_id')

local pass_neg1 = (real("`indirect_x1'") < 0)
di as txt "  indirect[x1] = `indirect_x1' (theta=-1, expect < 0)  " ///
   as `=cond(`pass_neg1', "txt", "err")' cond(`pass_neg1', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_neg1'
local NFAIL = `NFAIL' + (1 - `pass_neg1')

local pass_neg2 = (real("`indirect_x2'") < 0)
di as txt "  indirect[x2] = `indirect_x2' (theta=-1, expect < 0)  " ///
   as `=cond(`pass_neg2', "txt", "err")' cond(`pass_neg2', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_neg2'
local NFAIL = `NFAIL' + (1 - `pass_neg2')

// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Validation 3 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
