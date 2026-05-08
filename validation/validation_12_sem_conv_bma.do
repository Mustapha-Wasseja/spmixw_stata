*! Validation 12: SEM Convex BMA
*!
*! Mirror of validation_11 (SAR conv BMA) but for SEM convex.
*! DGP: u = (I - lambda * sum_m gamma_m W_m)^{-1} eps; y = X*beta + FE + u
*!     lambda_true = 0.5, beta = (1, 1)
*!     gamma_true = (0.0, 0.4, 0.6)  -- W1 has zero weight
*! Three banded W matrices (3 W => 7 non-empty subsets).
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

    function build_banded12(real scalar N, real scalar b) {
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

    W1 = build_banded12(N, 4)
    W2 = build_banded12(N, 6)
    W3 = build_banded12(N, 8)

    st_matrix("W1", W1)
    st_matrix("W2", W2)
    st_matrix("W3", W3)
end

mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    lambda_true     = 0.5
    beta_true       = (1 \ 1)
    gamma_true_full = (0 \ 0.4 \ 0.6)
    sigma2          = 1
    M = 3
    Wstack_3 = (st_matrix("W1"), st_matrix("W2"), st_matrix("W3"))
    y = .; X = .
    _spmixw_simulate_sem_conv(N, TT, Wstack_3, M,
                              lambda_true, gamma_true_full, beta_true,
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

print_header Validation 12: SEM Conv BMA -- 3 W matrices, true gamma=(0, 0.4, 0.6)

spmixw_bma y x1 x2, model(sem_conv) wmats(W1 W2 W3) effects(twoway) ///
    ndraw(8000) nomit(2000) seed(10203040) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

matrix probs   = e(probs)
matrix subsets = e(subsets)
local nm = e(nmodels)

local best_i = 1
local best_p = probs[1, 1]
forvalues i = 2/`nm' {
    if (probs[`i', 1] > `best_p') {
        local best_p = probs[`i', 1]
        local best_i = `i'
    }
}
di as txt _newline "Best model: `best_i' with posterior probability " ///
   as res %6.4f `best_p'

local has_W2 = subsets[`best_i', 2]
local has_W3 = subsets[`best_i', 3]

local pass_overlap = (`has_W2' == 1 | `has_W3' == 1)
di as txt "  Best model overlaps truth {W2, W3}: " ///
   as `=cond(`pass_overlap', "txt", "err")' cond(`pass_overlap', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_overlap'
local NFAIL = `NFAIL' + (1 - `pass_overlap')

local mass_excl_W1 = 0
forvalues i = 1/`nm' {
    if (subsets[`i', 1] == 0) {
        local mass_excl_W1 = `mass_excl_W1' + probs[`i', 1]
    }
}
di as txt "  Posterior mass on subsets without W1: " as res %6.4f `mass_excl_W1'

local pass_excl = (`mass_excl_W1' > 0.10)
di as txt "  P(W1 excluded) > 0.10: " ///
   as `=cond(`pass_excl', "txt", "err")' cond(`pass_excl', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_excl'
local NFAIL = `NFAIL' + (1 - `pass_excl')

matrix bma_g = e(bma_gamma)
local g1 = bma_g[1, 1]
local g2 = bma_g[1, 2]
local g3 = bma_g[1, 3]
di as txt "  BMA-weighted gamma = (`: di %5.3f `g1'', `: di %5.3f `g2'', " ///
    "`: di %5.3f `g3'')"

local pass_g1 = (`g1' < 0.30)
di as txt "  BMA gamma_1 < 0.30: " ///
   as `=cond(`pass_g1', "txt", "err")' cond(`pass_g1', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_g1'
local NFAIL = `NFAIL' + (1 - `pass_g1')

local pass_gsum = (`g2' + `g3' > 0.50)
di as txt "  BMA gamma_2 + gamma_3 > 0.50: " ///
   as `=cond(`pass_gsum', "txt", "err")' cond(`pass_gsum', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_gsum'
local NFAIL = `NFAIL' + (1 - `pass_gsum')

local pass_lam = (abs(e(bma_rho) - 0.5) < 0.20)
di as txt "  BMA lambda = " as res %6.4f e(bma_rho) ///
   as txt " vs truth 0.5: " ///
   as `=cond(`pass_lam', "txt", "err")' cond(`pass_lam', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_lam'
local NFAIL = `NFAIL' + (1 - `pass_lam')

di as txt _newline "{hline 70}"
di as txt "  Validation 12 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
