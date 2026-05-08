*! Validation 9: SDM Convex Panel
*!
*! DGP: y = (I - rho * Wc(gamma))^{-1} (X*beta + sum_m W_m*X*theta_m + FE + eps)
*! rho = 0.5, beta = (1, 1), theta_1 = (-0.5, -0.5), theta_2 = (-0.5, -0.5),
*! gamma = (0.4, 0.6), sige = 1
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

    function build_banded3(real scalar N, real scalar b) {
        real matrix    W
        real colvector rs
        real scalar    i, k

        W = J(N, N, 0)
        for (i = 1; i <= N; i++) {
            for (k = 1; k <= b; k++) {
                if (i - k >= 1) W[i, i-k] = 1
                if (i + k <= N) W[i, i+k] = 1
            }
        }
        rs = rowsum(W)
        for (i = 1; i <= N; i++) {
            if (rs[i] > 0) W[i, .] = W[i, .] :/ rs[i]
        }
        return(W)
    }

    W1 = build_banded3(N, 4)
    W2 = build_banded3(N, 6)

    st_matrix("W1", W1)
    st_matrix("W2", W2)
end

mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    nobs = N * TT
    rho_true    = 0.5
    beta_true   = (1 \ 1)
    theta1_true = (-0.5 \ -0.5)
    theta2_true = (-0.5 \ -0.5)
    gamma_true  = (0.4 \ 0.6)
    sigma2      = 1
    M = 2
    Wstack_2 = (st_matrix("W1"), st_matrix("W2"))

    rseed(10203040)
    X   = rnormal(nobs, 2, 0, 1)
    W1X = _spmixw_compute_wx(X, st_matrix("W1"), N, TT)
    W2X = _spmixw_compute_wx(X, st_matrix("W2"), N, TT)

    sfe = (1::N)  :/ N
    tfe = (1::TT) :/ TT
    fe  = (J(TT,1,1) # sfe) + (tfe # J(N,1,1))

    eps = rnormal(nobs, 1, 0, sqrt(sigma2))
    mu  = X * beta_true + W1X * theta1_true + W2X * theta2_true + fe + eps

    Wc = gamma_true[1] :* st_matrix("W1") + gamma_true[2] :* st_matrix("W2")
    A     = I(N) - rho_true :* Wc
    M_mat = rowshape(mu, TT)'
    Y_mat = lusolve(A, M_mat)
    y     = vec(Y_mat)

    region = J(TT, 1, 1) # (1::N)
    year   = (1::TT) # J(N, 1, 1)
    M_dat  = (region, year, y, X)

    st_addvar("int",    "region")
    st_addvar("int",    "year")
    st_addvar("double", "y")
    st_addvar("double", "x1")
    st_addvar("double", "x2")
    st_addobs(rows(M_dat))
    st_store(., ., M_dat)
end

xtset region year

local NPASS = 0
local NFAIL = 0

print_header Validation 9: SDM Conv Panel - M=2, true gamma=(0.4, 0.6)

spmixw y x1 x2, model(sdm_conv) wmats(W1 W2) effects(twoway) ///
    ndraw(15000) nomit(5000) seed(10203040) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

check_param, name("beta_1") estimate(`=_b[x1]') truth(1.0) tol(0.20)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2") estimate(`=_b[x2]') truth(1.0) tol(0.20)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("rho")    estimate(`=e(rho)')  truth(0.5) tol(0.20)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

matrix gam_pm = e(gam)
local g1 = gam_pm[1, 1]
local g2 = gam_pm[1, 2]

check_param, name("gamma_1") estimate(`g1') truth(0.4) tol(0.20)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("gamma_2") estimate(`g2') truth(0.6) tol(0.20)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

di as txt "  gamma posterior = (`: di %5.3f `g1'', `: di %5.3f `g2'')"
di as txt "  acceptance: rho=`: di %5.3f `=e(rho_acc_rate)''" ///
   "  gamma=`: di %5.3f `=e(gam_acc_rate)''"

di as txt _newline "{hline 70}"
di as txt "  Validation 9 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
