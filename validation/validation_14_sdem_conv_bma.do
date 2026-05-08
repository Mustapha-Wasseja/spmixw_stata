*! Validation 14: SDEM Convex BMA (smoke test, M=2)
*!
*! DGP: y = X*beta + sum_m W_m*X*theta_m + FE + u
*!      u = (I - lambda * Wc(gamma))^{-1} eps
*!      lambda = 0.5, beta = (1, 1), theta = (-0.5, -0.5) on each W
*!      gamma_true = (0.3, 0.7) on (W1, W2)
*! Two banded W matrices => 3 non-empty subsets.
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

    function build_banded14(real scalar N, real scalar b) {
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

    W1 = build_banded14(N, 4)
    W2 = build_banded14(N, 6)

    st_matrix("W1", W1)
    st_matrix("W2", W2)
end

mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    nobs = N * TT
    lambda_true = 0.5
    beta_true   = (1 \ 1)
    theta1_true = (-0.5 \ -0.5)
    theta2_true = (-0.5 \ -0.5)
    gamma_true  = (0.3 \ 0.7)
    sigma2      = 1
    M = 2

    rseed(10203040)
    X   = rnormal(nobs, 2, 0, 1)
    W1X = _spmixw_compute_wx(X, st_matrix("W1"), N, TT)
    W2X = _spmixw_compute_wx(X, st_matrix("W2"), N, TT)

    sfe = (1::N)  :/ N
    tfe = (1::TT) :/ TT
    fe  = (J(TT,1,1) # sfe) + (tfe # J(N,1,1))

    Wc = gamma_true[1] :* st_matrix("W1") + gamma_true[2] :* st_matrix("W2")
    A     = I(N) - lambda_true :* Wc
    eps_v = rnormal(nobs, 1, 0, sqrt(sigma2))
    Eps_M = rowshape(eps_v, TT)'
    U_M   = lusolve(A, Eps_M)
    u     = vec(U_M)

    y = X * beta_true + W1X * theta1_true + W2X * theta2_true + fe + u

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

print_header Validation 14: SDEM Conv BMA -- 2 W matrices, true gamma=(0.3, 0.7)

spmixw_bma y x1 x2, model(sdem_conv) wmats(W1 W2) effects(twoway) ///
    ndraw(8000) nomit(2000) seed(10203040) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

matrix bma_g = e(bma_gamma)
local g1 = bma_g[1, 1]
local g2 = bma_g[1, 2]
di as txt "  BMA-weighted gamma = (`: di %5.3f `g1'', `: di %5.3f `g2'')"

local pass_g1 = (abs(`g1' - 0.3) < 0.25)
di as txt "  BMA gamma_1 ~ 0.3: " ///
   as `=cond(`pass_g1', "txt", "err")' cond(`pass_g1', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_g1'
local NFAIL = `NFAIL' + (1 - `pass_g1')

local pass_g2 = (abs(`g2' - 0.7) < 0.25)
di as txt "  BMA gamma_2 ~ 0.7: " ///
   as `=cond(`pass_g2', "txt", "err")' cond(`pass_g2', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_g2'
local NFAIL = `NFAIL' + (1 - `pass_g2')

local pass_lam = (abs(e(bma_rho) - 0.5) < 0.25)
di as txt "  BMA lambda = " as res %6.4f e(bma_rho) ///
   as txt " vs truth 0.5: " ///
   as `=cond(`pass_lam', "txt", "err")' cond(`pass_lam', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_lam'
local NFAIL = `NFAIL' + (1 - `pass_lam')

matrix bma_t = e(bma_theta)
local nc = colsof(bma_t)
local pass_theta_shape = (`nc' == 4)
di as txt "  e(bma_theta) has 4 columns (M*k = 2*2): " ///
   as `=cond(`pass_theta_shape', "txt", "err")' ///
   cond(`pass_theta_shape', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_theta_shape'
local NFAIL = `NFAIL' + (1 - `pass_theta_shape')

di as txt _newline "{hline 70}"
di as txt "  Validation 14 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
