*! Validation 2: SAR Panel — Stata port of validation_02_sar.R
*!
*! DGP: y = (I - rho*Wbig)^{-1} (X*beta + FE + eps), rho=0.6, beta=(1,1), sige=1
*! W: banded contiguity (each region has neighbours {i-1, i+1}), row-normalised.
*! Effects checks (direct/indirect/total identity, total/direct ratio) are
*! deferred to Phase 2.B when `estat impact` is added.
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

// -- Generate SAR DGP via the package simulator ----------------------------
mata:
    N    = strtoreal(st_local("N"))
    TT   = strtoreal(st_local("T"))
    rho_true   = 0.6
    beta_true  = (1 \ 1)
    sigma2     = 1
    W = st_matrix("W")
    y = .; X = .
    _spmixw_simulate_sar(N, TT, W, rho_true, beta_true, sigma2, "twoway", y, X)

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
print_header Validation 2: SAR Panel - Homoscedastic, Two-way FE

spmixw y x1 x2, model(sar) w(W) effects(twoway) rval(0) ///
    ndraw(5000) nomit(1500) seed(10203040)

check_param, name("beta_1") estimate(`=_b[x1]') truth(1.0) tol(0.10) reference("0.9845")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("beta_2") estimate(`=_b[x2]') truth(1.0) tol(0.10) reference("1.0100")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("rho") estimate(`=e(rho)') truth(0.6) tol(0.05) reference("0.5990")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

check_param, name("sigma^2") estimate(`=e(sige)') truth(1.0) tol(0.30) reference("0.8694")
local NPASS = `NPASS' + r(pass)
local NFAIL = `NFAIL' + (1 - r(pass))

// =============================================================
// Effects checks (Phase 2.B.3): direct + indirect = total identity,
// total/direct ≈ 1/(1-rho) ≈ 2.5 for rho=0.6.
// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Effects checks via estat impact (SAR)"
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
    st_local("identity_err", strofreal(max_id))
    st_local("ratio_x1",     strofreal(tot_pm[1] / dir_pm[1]))
end

local pass_id = (real("`identity_err'") < 1e-8)
di as txt "  direct + indirect = total identity: max err = `identity_err'  " ///
   as `=cond(`pass_id', "txt", "err")' cond(`pass_id', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_id'
local NFAIL = `NFAIL' + (1 - `pass_id')

// For SAR, total/direct ∈ (1, 1/(1-ρ)] — the lower bound holds whenever
// indirect > 0 (always for positive β and 0 < ρ < 1); the upper bound is
// approached as tr(W^r)/N → 0. Our banded W has tr(W^2)/N ≈ 0.5, so the
// observed ratio is well below 1/(1-ρ); a W-specific tolerance like
// |ratio - 1/(1-ρ)| < 0.5 (which the R port uses with KNN W) does NOT
// hold. Use the W-agnostic bound instead.
local rho_pm = e(rho)
local upper = 1 / (1 - `rho_pm')
local r_x1  = real("`ratio_x1'")
local pass_ratio = (`r_x1' > 1.0) & (`r_x1' < `upper' + 0.05)
di as txt "  total/direct ratio for x1 = `: di %5.3f `r_x1''  " ///
   "(must be in (1.0, `: di %5.3f `upper''])  " ///
   as `=cond(`pass_ratio', "txt", "err")' cond(`pass_ratio', "PASS", "FAIL")
local NPASS = `NPASS' + `pass_ratio'
local NFAIL = `NFAIL' + (1 - `pass_ratio')

// =============================================================
di as txt _newline "{hline 70}"
di as txt "  Validation 2 Summary: `NPASS' passed, `NFAIL' failed"
di as txt "{hline 70}"
if `NFAIL' > 0 exit 9
