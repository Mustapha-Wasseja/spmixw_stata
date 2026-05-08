*! version 0.1.0  2026-05-08  spmixw -- OLS panel subroutine
*! Called by spmixw.ado; not intended for direct use.

program _spmixw_ols, eclass
    version 13.0

    syntax varlist(min=2 numeric) [if] [in],                ///
        ID(varname numeric)                                 ///
        TIME(varname numeric)                               ///
        MOD_(integer)                                       ///
        NDraw(integer)                                      ///
        NOmit(integer)                                      ///
        Thin(integer)                                       ///
        Rval(real)                                          ///
        Nu(real)                                            ///
        D0(real)                                            ///
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
        di as err "spmixw v0.1 requires a balanced panel (got N=`N', T=`TT', " ///
            "but `nobs' obs); unbalanced support is deferred."
        exit 459
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
    tempname b_pm V_pm draws_b draws_s vmean
    scalar sige_pm = .

    timer clear 1
    timer on 1

    mata: _spmixw_ols_caller(                          ///
        "`depvar'", "`covars'",                        ///
        `N', `TT',                                     ///
        `ndraw', `nomit', `thin', `mod_',              ///
        `rval', `nu', `d0',                            ///
        "`c_beta'", "`C_beta'",                        ///
        "`b_pm'", "`V_pm'",                            ///
        "`draws_b'", "`draws_s'", "`vmean'",           ///
        "sige_pm")

    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    timer clear 1

    // -- ereturn post --------------------------------------------------------
    matrix colnames `b_pm' = `covars'
    matrix colnames `V_pm' = `covars'
    matrix rownames `V_pm' = `covars'

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
    ereturn scalar time     = `elapsed'

    ereturn local depvar    `depvar'
    ereturn local covars    `covars'
    ereturn local id        `id'
    ereturn local timevar   `time'
    ereturn local cmd       "spmixw"
    ereturn local model     "ols"
    ereturn local effects   "`effname'"
    ereturn local estat_cmd "spmixw_estat"

    ereturn matrix draws_b  = `draws_b'
    ereturn matrix draws_s  = `draws_s'
    ereturn matrix vmean    = `vmean'

    // -- Optional saving() of raw draws --------------------------------------
    if `"`saving'"' != "" {
        preserve
        clear
        quietly svmat `draws_b', names(b)
        quietly svmat `draws_s', names(sigma2)
        quietly save `"`saving'"', replace
        restore
    }

    // -- Display -------------------------------------------------------------
    if "`nodisp'" == "" {
        _spmixw_ols_display, level(`level')
    }
end


// ----------------------------------------------------------------------
// Mata bridge — pulls Stata data, calls _spmixw_ols, posts results back
// ----------------------------------------------------------------------
version 13.0
mata:
mata set matastrict on

void _spmixw_ols_caller(string scalar  yvar,
                        string scalar  xvars,
                        real   scalar  N,
                        real   scalar  TT,
                        real   scalar  ndraw,
                        real   scalar  nomit,
                        real   scalar  thin,
                        real   scalar  mod_,
                        real   scalar  rval,
                        real   scalar  nu,
                        real   scalar  d0,
                        string scalar  c_beta_name,
                        string scalar  C_beta_name,
                        string scalar  b_pm_name,
                        string scalar  V_pm_name,
                        string scalar  draws_b_name,
                        string scalar  draws_s_name,
                        string scalar  vmean_name,
                        string scalar  sige_pm_scalar)
{
    real colvector y, beta_pm, sdraw, vmean, c_beta_col
    real matrix    X, bdraw, V_pm, c_beta_mat, C_beta_mat
    real scalar    sige_pm

    y = st_data(., yvar)
    X = st_data(., tokens(xvars))

    c_beta_mat = st_matrix(c_beta_name)
    C_beta_mat = st_matrix(C_beta_name)
    c_beta_col = c_beta_mat[., 1]

    _spmixw_ols(y, X, N, TT, ndraw, nomit, thin, mod_,
                rval, nu, d0, c_beta_col, C_beta_mat,
                bdraw, sdraw, vmean, beta_pm, sige_pm)

    // Posterior covariance from the retained draws (variance() takes columns)
    V_pm = variance(bdraw)

    st_matrix(b_pm_name,    beta_pm')
    st_matrix(V_pm_name,    V_pm)
    st_matrix(draws_b_name, bdraw)
    st_matrix(draws_s_name, sdraw)
    st_matrix(vmean_name,   vmean)

    st_numscalar(sige_pm_scalar, sige_pm)
}

end


// ----------------------------------------------------------------------
// Display helper — Bayesian-style table:
//   Post. Mean / Post. SD / [<level>% Cred. Interval] (quantile-based).
// ----------------------------------------------------------------------
program _spmixw_ols_display
    version 13.0
    syntax, [Level(cilevel)]

    tempname draws_b
    matrix `draws_b' = e(draws_b)

    di
    di as txt "Bayesian OLS panel (spmixw)" _col(58) "N obs    = " ///
        as res %9.0f e(N)
    di as txt "Model:    ols" _col(58) "N groups = " ///
        as res %9.0f e(N_g)
    di as txt "Effects:  " as res "`e(effects)'" _col(58) "T        = " ///
        as res %9.0f e(T)
    di as txt "MCMC:     " as res "ndraw=`=e(ndraw)' nomit=`=e(nomit)' " ///
        "thin=`=e(thin)'  rval=`=e(rval)'"
    di as txt "sigma^2 (post. mean) = " as res %9.4f e(sige) ///
        as txt "    elapsed = " as res %6.2f e(time) as txt " s"
    di

    _spmixw_print_table, draws(`draws_b') ///
        names(`e(covars)') depvar(`e(depvar)') level(`level')
end
