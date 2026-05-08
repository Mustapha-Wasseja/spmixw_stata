*! Validation 5: SDEM Panel — Stata port of validation_05_sdem.R
*!
*! DGP: u = (I - lambda*Wbig)^{-1} eps;  y = X*beta + WX*theta + FE + u
*! lambda = 0.7,  beta = (1, 1),  theta = (-1, -1),  sige = 0.1
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
    lambda_true = 0.7
    beta_true   = (1 \ 1)
    theta_true  = (-1 \ -1)
    sigma2      = 0.1

    rseed(10203040)
    X  = rnormal(nobs, 2, 0, 1)
    WX = _spmixw_compute_wx(X, W, N, TT)

    sfe = (1::N)  :/ N
    tfe = (1::TT) :/ TT
    fe  = (J(TT,1,1) # sfe) + (tfe # J(N,1,1))

    eps = rnormal(nobs, 1, 0, sqrt(sigma2))
    // u = (I - lambda*W)^{-1} eps, applied per period
    A       = I(N) - lambda_true :* W
    Eps_mat = rowshape(eps, TT)'
    U_mat   = lusolve(A, Eps_mat)
    u       = vec(U_mat)

    y = X * beta_true + WX * theta_true + fe + u

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
print_header Validation 5: SDEM Panel - Homoscedastic, Two-way FE

spmixw y x1 x2, model(sdem) w(W) effects(twoway) rval(0) ///
    ndraw(5000) nomit(1500) seed(10203040)

check_param, name("beta_1")  estimate(`=_b[x1]')   truth(1.0)  tol(0.10) reference("0.9885")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2")  estimate(`=_b[x2]')   truth(1.0)  tol(0.10) reference("1.0094")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("theta_1") estimate(`=_b[W_x1]') truth(-1.0) tol(0.15) reference("-1.0459")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("theta_2") estimate(`=_b[W_x2]') truth(-1.0) tol(0.15) reference("-0.9455")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("lambda")  estimate(`=e(rho)')   truth(0.7)  tol(0.08) reference("0.6880")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// =============================================================
// Effects checks (Phase 2.B.3): for SDEM the formula is exact —
// direct = beta, indirect = theta, total = beta + theta. So the
// effects table just rebroadcasts the coefficient table.
// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Effects checks via estat impact (SDEM)"
di as txt "{hline 70}"

estat impact

matrix dir_t = r(direct)
matrix ind_t = r(indirect)

local d1 = dir_t[1, 1]
local d2 = dir_t[2, 1]
local i1 = ind_t[1, 1]
local i2 = ind_t[2, 1]

local pass_d1 = (abs(`d1' - 1.0) < 0.10)
di as txt "  direct[x1] = `: di %7.4f `d1''  (β1=1.0)  " ///
   as `=cond(`pass_d1', "txt", "err")' cond(`pass_d1', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_d1'
local NFAIL = `NFAIL' + (1 - `pass_d1')

local pass_i1 = (abs(`i1' - (-1.0)) < 0.15)
di as txt "  indirect[x1] = `: di %7.4f `i1''  (θ1=-1.0)  " ///
   as `=cond(`pass_i1', "txt", "err")' cond(`pass_i1', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_i1'
local NFAIL = `NFAIL' + (1 - `pass_i1')

local pass_d2 = (abs(`d2' - 1.0) < 0.10)
di as txt "  direct[x2] = `: di %7.4f `d2''  (β2=1.0)  " ///
   as `=cond(`pass_d2', "txt", "err")' cond(`pass_d2', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_d2'
local NFAIL = `NFAIL' + (1 - `pass_d2')

local pass_i2 = (abs(`i2' - (-1.0)) < 0.15)
di as txt "  indirect[x2] = `: di %7.4f `i2''  (θ2=-1.0)  " ///
   as `=cond(`pass_i2', "txt", "err")' cond(`pass_i2', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_i2'
local NFAIL = `NFAIL' + (1 - `pass_i2')

// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Validation 5 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
