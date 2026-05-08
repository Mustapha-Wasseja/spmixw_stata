*! Validation 15: estat impact for sar_conv / sdm_conv / sdem_conv
*!
*! For each of the three convex models, fit, then call estat impact and
*! check basic invariants (identity direct+indirect=total) plus model-
*! specific sanity bounds.
*!
*! DGP: M=2 banded W's, true gamma=(0.4, 0.6).
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

    function build_banded15(real scalar N, real scalar b) {
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

    W1 = build_banded15(N, 4)
    W2 = build_banded15(N, 6)

    st_matrix("W1", W1)
    st_matrix("W2", W2)
end

local NPASS = 0
local NFAIL = 0

// ============================================================
// (a) SAR convex impacts
// ============================================================
mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    rho_true   = 0.6
    beta_true  = (1 \ 1)
    gamma_true = (0.4 \ 0.6)
    sigma2     = 1
    M = 2
    Wstack_2 = (st_matrix("W1"), st_matrix("W2"))
    y = .; X = .
    _spmixw_simulate_sar_conv(N, TT, Wstack_2, M,
                               rho_true, gamma_true, beta_true,
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

print_header Validation 15a: SAR conv estat impact

spmixw y x1 x2, model(sar_conv) wmats(W1 W2) effects(twoway) ///
    ndraw(10000) nomit(2500) seed(10203040) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

estat impact, seed(20260508)

matrix dir_t = r(direct)
matrix ind_t = r(indirect)
matrix tot_t = r(total)

mata:
    dir_pm = st_matrix("dir_t")[., 1]
    ind_pm = st_matrix("ind_t")[., 1]
    tot_pm = st_matrix("tot_t")[., 1]
    max_id = max(abs(dir_pm + ind_pm - tot_pm))
    st_local("identity_err", strofreal(max_id))
    st_local("ratio_x1",     strofreal(tot_pm[1] / dir_pm[1]))
    st_local("dir_x1",       strofreal(dir_pm[1]))
    st_local("ind_x1",       strofreal(ind_pm[1]))
end

local pass_id = (real("`identity_err'") < 1e-8)
di as txt "  SAR conv: direct + indirect = total: max err = `identity_err'  " ///
   as `=cond(`pass_id', "txt", "err")' cond(`pass_id', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_id'
local NFAIL = `NFAIL' + (1 - `pass_id')

local rho_pm = e(rho)
local upper = 1 / (1 - `rho_pm')
local r_x1  = real("`ratio_x1'")
local pass_ratio = (`r_x1' > 1.0) & (`r_x1' < `upper' + 0.05)
di as txt "  SAR conv: total/direct ratio for x1 = `: di %6.3f `r_x1'' " ///
   as txt "(in (1, " %5.3f `upper' "]): " ///
   as `=cond(`pass_ratio', "txt", "err")' cond(`pass_ratio', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_ratio'
local NFAIL = `NFAIL' + (1 - `pass_ratio')

di as txt "  SAR conv: direct(x1) = `: di %6.3f `dir_x1''  indirect(x1) = `: di %6.3f `ind_x1''"

// ============================================================
// (b) SDEM convex impacts
// DGP: y = X*beta + sum_m W_m*X*theta_m + FE + u
//      u = (I - lambda * Wc(gamma))^{-1} eps
// theta_1 = theta_2 = (-0.5, -0.5).  Indirect_j (truth) = sum_m theta_{m,j} = -1.
// ============================================================
clear all
quietly do "${SPMIXW_PKG}/tests/_setup.do"
quietly do "${SPMIXW_PKG}/validation/validation_helpers.do"

set seed 10203040

mata:
    function build_banded15b(real scalar N, real scalar b) {
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

    N = 200
    W1 = build_banded15b(N, 4)
    W2 = build_banded15b(N, 6)
    st_matrix("W1", W1)
    st_matrix("W2", W2)
end

mata:
    N    = 200
    TT   = 10
    nobs = N * TT
    lambda_true = 0.5
    beta_true   = (1 \ 1)
    theta1_true = (-0.5 \ -0.5)
    theta2_true = (-0.5 \ -0.5)
    gamma_true  = (0.4 \ 0.6)
    sigma2      = 1

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

print_header Validation 15b: SDEM conv estat impact

spmixw y x1 x2, model(sdem_conv) wmats(W1 W2) effects(twoway) ///
    ndraw(10000) nomit(2500) seed(10203040) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

estat impact, seed(20260508)

matrix dir_t = r(direct)
matrix ind_t = r(indirect)
matrix tot_t = r(total)

mata:
    dir_pm = st_matrix("dir_t")[., 1]
    ind_pm = st_matrix("ind_t")[., 1]
    tot_pm = st_matrix("tot_t")[., 1]
    max_id = max(abs(dir_pm + ind_pm - tot_pm))
    st_local("identity_err", strofreal(max_id))
    st_local("dir_x1",       strofreal(dir_pm[1]))
    st_local("ind_x1",       strofreal(ind_pm[1]))
    st_local("dir_x2",       strofreal(dir_pm[2]))
    st_local("ind_x2",       strofreal(ind_pm[2]))
end

local pass_id_b = (real("`identity_err'") < 1e-8)
di as txt "  SDEM conv: direct + indirect = total: max err = `identity_err'  " ///
   as `=cond(`pass_id_b', "txt", "err")' cond(`pass_id_b', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_id_b'
local NFAIL = `NFAIL' + (1 - `pass_id_b')

local pass_dir_b = (abs(real("`dir_x1'") - 1.0) < 0.20) & ///
                   (abs(real("`dir_x2'") - 1.0) < 0.20)
di as txt "  SDEM conv: direct(x1)=`: di %6.3f real("`dir_x1'")' " ///
   "direct(x2)=`: di %6.3f real("`dir_x2'")' (truth=1.0): " ///
   as `=cond(`pass_dir_b', "txt", "err")' cond(`pass_dir_b', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_dir_b'
local NFAIL = `NFAIL' + (1 - `pass_dir_b')

local pass_ind_b = (abs(real("`ind_x1'") - (-1.0)) < 0.30) & ///
                   (abs(real("`ind_x2'") - (-1.0)) < 0.30)
di as txt "  SDEM conv: indirect(x1)=`: di %6.3f real("`ind_x1'")' " ///
   "indirect(x2)=`: di %6.3f real("`ind_x2'")' (truth=-1.0): " ///
   as `=cond(`pass_ind_b', "txt", "err")' cond(`pass_ind_b', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_ind_b'
local NFAIL = `NFAIL' + (1 - `pass_ind_b')

// ============================================================
// (c) SDM convex impacts
// DGP: y = (I - rho * Wc(gamma))^{-1} (X*beta + sum_m W_m*X*theta_m + FE + eps)
// rho_true = 0.5, theta_1 = theta_2 = (-0.5, -0.5).
// ============================================================
clear all
quietly do "${SPMIXW_PKG}/tests/_setup.do"
quietly do "${SPMIXW_PKG}/validation/validation_helpers.do"

set seed 10203040

mata:
    function build_banded15c(real scalar N, real scalar b) {
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

    N = 200
    W1 = build_banded15c(N, 4)
    W2 = build_banded15c(N, 6)
    st_matrix("W1", W1)
    st_matrix("W2", W2)
end

mata:
    N    = 200
    TT   = 10
    nobs = N * TT
    rho_true    = 0.5
    beta_true   = (1 \ 1)
    theta1_true = (-0.5 \ -0.5)
    theta2_true = (-0.5 \ -0.5)
    gamma_true  = (0.4 \ 0.6)
    sigma2      = 1

    rseed(10203040)
    X   = rnormal(nobs, 2, 0, 1)
    W1X = _spmixw_compute_wx(X, st_matrix("W1"), N, TT)
    W2X = _spmixw_compute_wx(X, st_matrix("W2"), N, TT)

    sfe = (1::N)  :/ N
    tfe = (1::TT) :/ TT
    fe  = (J(TT,1,1) # sfe) + (tfe # J(N,1,1))

    eps = rnormal(nobs, 1, 0, sqrt(sigma2))
    mu  = X * beta_true + W1X * theta1_true + W2X * theta2_true + fe + eps

    Wc    = gamma_true[1] :* st_matrix("W1") + gamma_true[2] :* st_matrix("W2")
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

print_header Validation 15c: SDM conv estat impact

spmixw y x1 x2, model(sdm_conv) wmats(W1 W2) effects(twoway) ///
    ndraw(10000) nomit(2500) seed(10203040) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

estat impact, seed(20260508)

matrix dir_t = r(direct)
matrix ind_t = r(indirect)
matrix tot_t = r(total)

mata:
    dir_pm = st_matrix("dir_t")[., 1]
    ind_pm = st_matrix("ind_t")[., 1]
    tot_pm = st_matrix("tot_t")[., 1]
    max_id = max(abs(dir_pm + ind_pm - tot_pm))
    st_local("identity_err", strofreal(max_id))
    st_local("tot_x1",       strofreal(tot_pm[1]))
    st_local("dir_x1",       strofreal(dir_pm[1]))
end

local pass_id_c = (real("`identity_err'") < 1e-8)
di as txt "  SDM conv: direct + indirect = total: max err = `identity_err'  " ///
   as `=cond(`pass_id_c', "txt", "err")' cond(`pass_id_c', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_id_c'
local NFAIL = `NFAIL' + (1 - `pass_id_c')

// Structural check: total_j must equal (beta_j + sum_m theta_{m,j}) /
// (1 - rho), evaluated at posterior means (modulo small Jensen-type
// discrepancy from per-draw averaging). This tests the impact formula, not
// the MCMC's recovery of the underlying coefficients.
//
// e(b) layout for sdm_conv (after sar_conv routing): [beta (k=2),
// W1_x1, W1_x2, W2_x1, W2_x2 (theta blocks), Wy (rho)].
matrix b_pm = e(b)
local beta1   = el(b_pm, 1, 1)
local th11    = el(b_pm, 1, 3)   // W1_x1
local th21    = el(b_pm, 1, 5)   // W2_x1
local rho_pm  = e(rho)
local pred_tot_x1 = (`beta1' + `th11' + `th21') / (1 - `rho_pm')
local tot_x1      = real("`tot_x1'")

local pass_struct = (abs(`tot_x1' - `pred_tot_x1') < 0.10 * max(abs(`pred_tot_x1'), 0.5))
di as txt "  SDM conv: total(x1)=`: di %6.3f `tot_x1'' vs " ///
   "(beta+sum_theta)/(1-rho)=`: di %6.3f `pred_tot_x1'' " ///
   "(structural check): " ///
   as `=cond(`pass_struct', "txt", "err")' cond(`pass_struct', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_struct'
local NFAIL = `NFAIL' + (1 - `pass_struct')

di as txt "  SDM conv: direct(x1) = `: di %6.3f real("`dir_x1'")'  " ///
   "indirect(x1) = `: di %6.3f `tot_x1' - real("`dir_x1'")'"

di as txt _newline "{hline 70}"
di as txt "  Validation 15 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
