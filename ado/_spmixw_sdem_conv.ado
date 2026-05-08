*! version 0.1.0  2026-05-08  spmixw -- SDEM convex panel wrapper
*!
*! Augments X with [W_1*X, W_2*X, ..., W_M*X] then routes to sem_conv.

program _spmixw_sdem_conv, eclass
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

    tempvar tagid tagt
    quietly egen `tagid' = tag(`id')   if `touse'
    quietly egen `tagt'  = tag(`time') if `touse'
    quietly count if `tagid' == 1
    local N = r(N)
    quietly count if `tagt' == 1
    local TT = r(N)

    local M : word count `wmats'
    if `M' < 1 {
        di as err "model(sdem_conv) requires at least one matrix in wmats()"
        exit 198
    }

    preserve
    quietly keep if `touse'
    sort `time' `id'

    local wxvars
    local m_idx = 0
    foreach Wm of local wmats {
        local m_idx = `m_idx' + 1
        local block_vars
        foreach x of local covars {
            local newname "W`m_idx'_`x'"
            capture confirm new variable `newname'
            if _rc {
                di as err "spmixw, model(sdem_conv): cannot create variable " ///
                    "`newname' -- already exists in dataset. Rename and retry."
                exit 110
            }
            quietly gen double `newname' = .
            local block_vars `block_vars' `newname'
        }
        mata: _spmixw_fill_wx_in_data("`Wm'", "`covars'", "`block_vars'", `N', `TT')
        local wxvars `wxvars' `block_vars'
    }

    _spmixw_sem_conv `depvar' `covars' `wxvars',                       ///
        id(`id') time(`time') mod_(`mod_')                             ///
        wmats(`wmats')                                                 ///
        ndraw(`ndraw') nomit(`nomit') thin(`thin')                     ///
        rmin(`rmin') rmax(`rmax')                                      ///
        taylororder(`taylororder')                                     ///
        saving(`"`saving'"') level(`level')                            ///
        nodisp

    ereturn local model       "sdem_conv"
    ereturn local covars_orig "`covars'"
    ereturn local covars_wx   "`wxvars'"

    restore

    if "`nodisp'" == "" {
        _spmixw_sdem_conv_display, level(`level')
    }
end


program _spmixw_sdem_conv_display
    version 13.0
    syntax, [Level(cilevel)]

    tempname draws_full
    matrix `draws_full' = e(draws_full)

    di
    di as txt "Bayesian SDEM (convex W) panel (spmixw)" _col(58) "N obs    = " ///
        as res %9.0f e(N)
    di as txt "Model:    sdem_conv" _col(58) "N groups = " ///
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
