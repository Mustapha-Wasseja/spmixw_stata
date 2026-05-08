*! version 0.1.0  2026-05-08  spmixw -- SEM convex panel subroutine
*!
*! Bayesian SEM with W_c(gamma) = sum_m gamma_m * W_m. Multi-W via wmats().
*!
*! Calls the same Mata MCMC kernel (_spmixw_sar_conv_caller, in lspmixw.mlib)
*! as the SAR convex ado -- the math is identical, only the user-facing
*! labels differ (row label "We" instead of "Wy", header "SEM (convex W)",
*! e(model) = "sem_conv").

program _spmixw_sem_conv, eclass
    version 13.0

    syntax varlist(min=2 numeric) [if] [in],                ///
        ID(varname numeric)                                 ///
        TIME(varname numeric)                               ///
        MOD_(integer)                                       ///
        WMats(string)                                       ///
        NDraw(integer)                                      ///
        NOmit(integer)                                      ///
        Thin(integer)                                       ///
        RMIn(real)                                          ///
        RMAx(real)                                          ///
        TAYLORorder(integer)                                ///
        [                                                   ///
            SAving(string)                                  ///
            Level(cilevel)                                  ///
            NODISP                                          ///
        ]

    marksample touse
    markout `touse' `id' `time'

    gettoken depvar covars : varlist

    local M : word count `wmats'
    if `M' < 1 {
        di as err "model(sem_conv) requires at least one matrix in wmats()"
        exit 198
    }
    foreach Wm of local wmats {
        capture confirm matrix `Wm'
        if _rc {
            di as err "wmats: '`Wm'' is not an existing Stata matrix"
            exit 198
        }
    }

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
        di as err "spmixw v0.3 requires a balanced panel (got N=`N', T=`TT', " ///
            "but `nobs' obs)."
        exit 459
    }

    foreach Wm of local wmats {
        if rowsof(`Wm') != `N' | colsof(`Wm') != `N' {
            di as err "wmats: `Wm' must be `N' x `N' (matches xtset panel size)"
            exit 198
        }
    }

    preserve
    quietly keep if `touse'
    sort `time' `id'

    local k : word count `covars'

    tempname b_pm V_pm draws_b draws_s draws_p draws_g draws_full vmean
    tempname gam_pm
    scalar sige_pm = .
    scalar rho_pm  = .
    scalar rho_acc = .
    scalar gam_acc = .
    scalar log_marg = .

    timer clear 1
    timer on 1

    mata: _spmixw_sar_conv_caller(                  ///
        "`depvar'", "`covars'",                     ///
        `M', `N', `TT',                             ///
        "`wmats'",                                  ///
        `ndraw', `nomit', `thin', `mod_',           ///
        `rmin', `rmax', `taylororder',              ///
        "`b_pm'", "`V_pm'",                         ///
        "`draws_b'", "`draws_s'", "`draws_p'",      ///
        "`draws_g'", "`draws_full'",                ///
        "`gam_pm'",                                 ///
        "sige_pm", "rho_pm", "rho_acc", "gam_acc", "log_marg")

    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    timer clear 1

    matrix colnames `b_pm' = `covars' We
    matrix colnames `V_pm' = `covars' We
    matrix rownames `V_pm' = `covars' We

    restore

    ereturn post `b_pm' `V_pm', depname(`depvar') esample(`touse')

    local effname "none"
    if `mod_' == 1 local effname "region"
    if `mod_' == 2 local effname "time"
    if `mod_' == 3 local effname "twoway"

    ereturn scalar N        = `nobs'
    ereturn scalar N_g      = `N'
    ereturn scalar T        = `TT'
    ereturn scalar k        = `k'
    ereturn scalar M        = `M'
    ereturn scalar ndraw    = `ndraw'
    ereturn scalar nomit    = `nomit'
    ereturn scalar thin     = `thin'
    ereturn scalar mod_     = `mod_'
    ereturn scalar sige     = sige_pm
    ereturn scalar rho      = rho_pm
    ereturn scalar rmin     = `rmin'
    ereturn scalar rmax     = `rmax'
    ereturn scalar taylor_order = `taylororder'
    ereturn scalar rho_acc_rate = rho_acc
    ereturn scalar gam_acc_rate = gam_acc
    ereturn scalar logmarg      = log_marg
    ereturn scalar time     = `elapsed'

    ereturn local depvar    `depvar'
    ereturn local covars    `covars'
    ereturn local id        `id'
    ereturn local timevar   `time'
    ereturn local cmd       "spmixw"
    ereturn local model     "sem_conv"
    ereturn local effects   "`effname'"
    ereturn local wmats     "`wmats'"
    ereturn local estat_cmd "spmixw_estat"

    ereturn matrix draws_b     = `draws_b'
    ereturn matrix draws_s     = `draws_s'
    ereturn matrix draws_p     = `draws_p'
    ereturn matrix draws_g     = `draws_g'
    ereturn matrix draws_full  = `draws_full'
    ereturn matrix gam         = `gam_pm'

    if `"`saving'"' != "" {
        preserve
        clear
        quietly svmat `draws_b', names(b)
        quietly svmat `draws_s', names(sigma2)
        quietly svmat `draws_p', names(lambda)
        quietly svmat `draws_g', names(gam)
        quietly save `"`saving'"', replace
        restore
    }

    if "`nodisp'" == "" {
        _spmixw_sem_conv_display, level(`level')
    }
end


// ----------------------------------------------------------------------
// Display helper -- SEM-flavoured header, We row, gamma block.
// ----------------------------------------------------------------------
program _spmixw_sem_conv_display
    version 13.0
    syntax, [Level(cilevel)]

    tempname draws_full
    matrix `draws_full' = e(draws_full)

    di
    di as txt "Bayesian SEM (convex W) panel (spmixw)" _col(58) "N obs    = " ///
        as res %9.0f e(N)
    di as txt "Model:    sem_conv" _col(58) "N groups = " ///
        as res %9.0f e(N_g)
    di as txt "Effects:  " as res "`e(effects)'" _col(58) "T        = " ///
        as res %9.0f e(T)
    di as txt "MCMC:     " as res "ndraw=`=e(ndraw)' nomit=`=e(nomit)' " ///
        "thin=`=e(thin)'  taylor_order=`=e(taylor_order)'"
    di as txt "MH acceptance: " as res "lam=" %5.3f e(rho_acc_rate) ///
        "  gam=" %5.3f e(gam_acc_rate)
    di as txt "sigma^2 (post. mean) = " as res %9.4f e(sige) ///
        as txt "    elapsed = " as res %6.2f e(time) as txt " s"
    di

    _spmixw_print_table, draws(`draws_full') ///
        names(`e(covars)' We) depvar(`e(depvar)') level(`level')

    tempname draws_g
    matrix `draws_g' = e(draws_g)

    local M = e(M)
    local gnames
    forvalues m = 1/`M' {
        local gnames `gnames' gam`m'
    }
    di
    di as txt "Convex weights gamma (M = `M'):"
    _spmixw_print_table, draws(`draws_g') names(`gnames') ///
        depvar("gamma") level(`level')
end
