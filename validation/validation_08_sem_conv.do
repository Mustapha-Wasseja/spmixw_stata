*! Validation 8: SEM Convex Panel
*!
*! DGP: u = (I - lambda * sum_m gamma_m W_m)^{-1} eps; y = X*beta + FE + u
*! lambda = 0.5, beta = (1, 1), gamma = (0.4, 0.6), sige = 1
*! Two banded W matrices.
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

    function build_banded2(real scalar N, real scalar b) {
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

    W1 = build_banded2(N, 4)
    W2 = build_banded2(N, 6)

    st_matrix("W1", W1)
    st_matrix("W2", W2)
end

mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    lambda_true = 0.5
    beta_true   = (1 \ 1)
    gamma_true  = (0.4 \ 0.6)
    sigma2      = 1
    M = 2
    Wstack_2 = (st_matrix("W1"), st_matrix("W2"))
    y = .; X = .
    _spmixw_simulate_sem_conv(N, TT, Wstack_2, M,
                              lambda_true, gamma_true, beta_true,
                              sigma2, "twoway", y, X)

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

print_header Validation 8: SEM Conv Panel - M=2, true gamma=(0.4, 0.6)

spmixw y x1 x2, model(sem_conv) wmats(W1 W2) effects(twoway) ///
    ndraw(15000) nomit(5000) seed(10203040) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

check_param, name("beta_1") estimate(`=_b[x1]') truth(1.0) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2") estimate(`=_b[x2]') truth(1.0) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("lambda") estimate(`=e(rho)') truth(0.5) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

matrix gam_pm = e(gam)
local g1 = gam_pm[1, 1]
local g2 = gam_pm[1, 2]

check_param, name("gamma_1") estimate(`g1') truth(0.4) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("gamma_2") estimate(`g2') truth(0.6) tol(0.15)
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

di as txt "  gamma posterior = (`: di %5.3f `g1'', `: di %5.3f `g2'')"
di as txt "  acceptance: lam=`: di %5.3f `=e(rho_acc_rate)''" ///
   "  gamma=`: di %5.3f `=e(gam_acc_rate)''"

local pass_acc_l = (e(rho_acc_rate) > 0.10) & (e(rho_acc_rate) < 0.80)
di as txt "  lambda acc in (0.10, 0.80): " ///
   as `=cond(`pass_acc_l', "txt", "err")' cond(`pass_acc_l', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_acc_l'
local NFAIL = `NFAIL' + (1 - `pass_acc_l')

local pass_acc_g = (e(gam_acc_rate) > 0.03)
di as txt "  gamma acc > 0.03: " ///
   as `=cond(`pass_acc_g', "txt", "err")' cond(`pass_acc_g', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_acc_g'
local NFAIL = `NFAIL' + (1 - `pass_acc_g')

di as txt _newline "{hline 70}"
di as txt "  Validation 8 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
