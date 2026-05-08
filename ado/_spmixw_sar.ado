*! version 0.1.0  2026-05-08  spmixw — SAR panel subroutine
*! Called by spmixw.ado; not intended for direct use.

program _spmixw_sar, eclass
    version 13.0

    syntax varlist(min=2 numeric) [if] [in],                ///
        ID(varname numeric)                                 ///
        TIME(varname numeric)                               ///
        MOD_(integer)                                       ///
        W(string)                                           ///
        NDraw(integer)                                      ///
        NOmit(integer)                                      ///
        Thin(integer)                                       ///
        Rval(real)                                          ///
        Nu(real)                                            ///
        D0(real)                                            ///
        RMIn(real)                                          ///
        RMAx(real)                                          ///
        GRIDStep(real)                                      ///
        [                                                   ///
            BPrior(string)                                  ///
            BVar(string)                                    ///
            SAving(string)                                  ///
            Level(cilevel)                                  ///
            NODISP                                          ///
        ]

    marksample touse
    markout `touse' `id' `time'

    gettoken depvar covars : varlist

    // -- Establish panel dimensions on the active sample ---------------------
    tempvar tagid tagt
    quietly egen `tagid' = tag(`id')   if `touse'
    quietly egen `tagt'  = tag(`time') if `touse'

    quietly count if `tagid' == 1
    local N = r(N)

    quietly count if `tagt' == 1
    local TT = r(N)

    quietly count if `touse'
    local nobs = r(N)
    if `nobs' != `N' * `TT' {
        di as err "spmixw v0.2 requires a balanced panel (got N=`N', T=`TT', " ///
            "but `nobs' obs); unbalanced support is deferred."
        exit 459
    }

    // -- Validate W dimensions -----------------------------------------------
    if rowsof(`w') != `N' | colsof(`w') != `N' {
        di as err "w(`w') must be `N' x `N' (matches the number of panels)"
        exit 198
    }

    // -- Sort to (time, id) so observation (i, t) sits at row i + (t-1)*N ----
    preserve
    quietly keep if `touse'
    sort `time' `id'

    local k : word count `covars'

    // -- Optional priors as Stata matrices (else use diffuse defaults) -------
    tempname c_beta C_beta
    if "`bprior'" != "" {
        capture confirm matrix `bprior'
        if _rc {
            di as err "bprior(`bprior') is not an existing matrix"
            exit 198
        }
        matrix `c_beta' = `bprior'
        if rowsof(`c_beta') != `k' | colsof(`c_beta') != 1 {
            di as err "bprior must be `k' x 1"
            exit 198
        }
    }
    else {
        matrix `c_beta' = J(`k', 1, 0)
    }

    if "`bvar'" != "" {
        capture confirm matrix `bvar'
        if _rc {
            di as err "bvar(`bvar') is not an existing matrix"
            exit 198
        }
        matrix `C_beta' = `bvar'
        if rowsof(`C_beta') != `k' | colsof(`C_beta') != `k' {
            di as err "bvar must be `k' x `k'"
            exit 198
        }
    }
    else {
        matrix `C_beta' = I(`k') * 1e12
    }

    // -- Hand off to Mata (timed) --------------------------------------------
    tempname b_pm V_pm draws_b draws_s draws_p draws_full vmean detval
    scalar sige_pm = .
    scalar rho_pm  = .

    timer clear 1
    timer on 1

    mata: _spmixw_sar_caller(                          ///
        "`depvar'", "`covars'",                        ///
        `N', `TT',                                     ///
        "`w'",                                         ///
        `ndraw', `nomit', `thin', `mod_',              ///
        `rval', `nu', `d0',                            ///
        `rmin', `rmax', `gridstep',                    ///
        "`c_beta'", "`C_beta'",                        ///
        "`b_pm'", "`V_pm'",                            ///
        "`draws_b'", "`draws_s'", "`draws_p'",         ///
        "`draws_full'",                                ///
        "`vmean'", "`detval'",                         ///
        "sige_pm", "rho_pm")

    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    timer clear 1

    // -- Build the joint (beta, rho) coefficient vector + covariance --------
    // ρ enters the SAR specification as the coefficient on Wy, so it belongs
    // in the coefficient table. The Mata bridge built `b_pm` and `V_pm`
    // already enlarged to (k+1) — last entry is Wy. Set names accordingly.
    matrix colnames `b_pm' = `covars' Wy
    matrix colnames `V_pm' = `covars' Wy
    matrix rownames `V_pm' = `covars' Wy

    restore

    ereturn post `b_pm' `V_pm', depname(`depvar') esample(`touse')

    // FE label
    local effname "none"
    if `mod_' == 1 local effname "region"
    if `mod_' == 2 local effname "time"
    if `mod_' == 3 local effname "twoway"

    ereturn scalar N        = `nobs'
    ereturn scalar N_g      = `N'
    ereturn scalar T        = `TT'
    ereturn scalar k        = `k'
    ereturn scalar ndraw    = `ndraw'
    ereturn scalar nomit    = `nomit'
    ereturn scalar thin     = `thin'
    ereturn scalar rval     = `rval'
    ereturn scalar nu       = `nu'
    ereturn scalar d0       = `d0'
    ereturn scalar mod_     = `mod_'
    ereturn scalar sige     = sige_pm
    ereturn scalar rho      = rho_pm
    ereturn scalar rmin     = `rmin'
    ereturn scalar rmax     = `rmax'
    ereturn scalar gridstep = `gridstep'
    ereturn scalar time     = `elapsed'

    ereturn local depvar    `depvar'
    ereturn local covars    `covars'
    ereturn local id        `id'
    ereturn local timevar   `time'
    ereturn local cmd       "spmixw"
    ereturn local model     "sar"
    ereturn local effects   "`effname'"
    ereturn local wmatrix   "`w'"
    ereturn local estat_cmd "spmixw_estat"

    ereturn matrix draws_b     = `draws_b'
    ereturn matrix draws_s     = `draws_s'
    ereturn matrix draws_p     = `draws_p'
    ereturn matrix draws_full  = `draws_full'
    ereturn matrix vmean       = `vmean'
    ereturn matrix lndet       = `detval'

    // -- Optional saving() of raw draws --------------------------------------
    if `"`saving'"' != "" {
        preserve
        clear
        quietly svmat `draws_b', names(b)
        quietly svmat `draws_s', names(sigma2)
        quietly svmat `draws_p', names(rho)
        quietly save `"`saving'"', replace
        restore
    }

    // -- Display -------------------------------------------------------------
    if "`nodisp'" == "" {
        _spmixw_sar_display, level(`level')
    }
end


// ----------------------------------------------------------------------
// Mata bridge — pulls Stata data, builds detval, calls _spmixw_sar
// ----------------------------------------------------------------------
version 13.0
mata:
mata set matastrict on

void _spmixw_sar_caller(string scalar  yvar,
                        string scalar  xvars,
                        real   scalar  N,
                        real   scalar  TT,
                        string scalar  W_name,
                        real   scalar  ndraw,
                        real   scalar  nomit,
                        real   scalar  thin,
                        real   scalar  mod_,
                        real   scalar  rval,
                        real   scalar  nu,
                        real   scalar  d0,
                        real   scalar  rmin,
                        real   scalar  rmax,
                        real   scalar  gridstep,
                        string scalar  c_beta_name,
                        string scalar  C_beta_name,
                        string scalar  b_pm_name,
                        string scalar  V_pm_name,
                        string scalar  draws_b_name,
                        string scalar  draws_s_name,
                        string scalar  draws_p_name,
                        string scalar  draws_full_name,
                        string scalar  vmean_name,
                        string scalar  detval_name,
                        string scalar  sige_pm_scalar,
                        string scalar  rho_pm_scalar)
{
    real colvector y, beta_pm, sdraw, pdraw, vmean, c_beta_col, b_full
    real matrix    X, W, bdraw, V_full, c_beta_mat, C_beta_mat, draws_full
    real matrix    detval, detval_panel
    real scalar    sige_pm, rho_pm

    y = st_data(., yvar)
    X = st_data(., tokens(xvars))
    W = st_matrix(W_name)

    c_beta_mat = st_matrix(c_beta_name)
    C_beta_mat = st_matrix(C_beta_name)
    c_beta_col = c_beta_mat[., 1]

    // Build the per-unit log-determinant grid, then panel-scale by T.
    detval = _spmixw_logdet_exact(W, rmin, rmax, gridstep)
    detval_panel = detval
    detval_panel[., 2] = TT :* detval[., 2]

    _spmixw_sar(y, X, W, N, TT, ndraw, nomit, thin, mod_,
                rval, nu, d0, c_beta_col, C_beta_mat,
                detval_panel,
                bdraw, sdraw, pdraw, vmean, beta_pm, sige_pm, rho_pm)

    // Joint (beta, rho) draws so the coefficient table can show ρ as Wy.
    draws_full = (bdraw, pdraw)
    b_full     = (beta_pm \ rho_pm)
    V_full     = variance(draws_full)

    st_matrix(b_pm_name,    b_full')
    st_matrix(V_pm_name,    V_full)
    st_matrix(draws_b_name, bdraw)
    st_matrix(draws_s_name, sdraw)
    st_matrix(draws_p_name, pdraw)
    st_matrix(draws_full_name, draws_full)
    st_matrix(vmean_name,   vmean)
    st_matrix(detval_name,  detval)

    st_numscalar(sige_pm_scalar, sige_pm)
    st_numscalar(rho_pm_scalar,  rho_pm)
}

end


// ----------------------------------------------------------------------
// Display helper — Bayesian-style table:
//   Coef. -> Post. Mean
//   Std. Err. -> Post. SD
//   [95% Conf. Interval] -> [95% Cred. Interval] (quantile-based)
//   Drops the misleading z and P>|z| columns.
//   ρ appears in the table as the row labelled "Wy".
// ----------------------------------------------------------------------
program _spmixw_sar_display
    version 13.0
    syntax, [Level(cilevel)]

    tempname draws_full
    matrix `draws_full' = e(draws_full)

    di
    di as txt "Bayesian SAR panel (spmixw)" _col(58) "N obs    = " ///
        as res %9.0f e(N)
    di as txt "Model:    sar" _col(58) "N groups = " ///
        as res %9.0f e(N_g)
    di as txt "Effects:  " as res "`e(effects)'" _col(58) "T        = " ///
        as res %9.0f e(T)
    di as txt "MCMC:     " as res "ndraw=`=e(ndraw)' nomit=`=e(nomit)' " ///
        "thin=`=e(thin)'  rval=`=e(rval)'"
    di as txt "sigma^2 (post. mean) = " as res %9.4f e(sige) ///
        as txt "    elapsed = " as res %6.2f e(time) as txt " s"
    di

    _spmixw_print_table, draws(`draws_full') ///
        names(`e(covars)' Wy) depvar(`e(depvar)') level(`level')
end
