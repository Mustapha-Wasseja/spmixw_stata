*! version 0.1.0  2026-05-08  spmixw — SLX panel wrapper
*! Augments X with WX, then routes to `_spmixw_ols` (no spatial parameter).
*! Coefficient table contains [original X | W_<X>] only — no Wy / We.

program _spmixw_slx, eclass
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
        [                                                   ///
            BPrior(string)                                  ///
            BVar(string)                                    ///
            SAving(string)                                  ///
            Level(cilevel)                                  ///
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

    preserve
    quietly keep if `touse'
    sort `time' `id'

    local k : word count `covars'
    local wxvars
    foreach x of local covars {
        capture confirm new variable W_`x'
        if _rc {
            di as err "spmixw, model(slx): cannot create variable W_`x'" ///
                " — already exists in dataset. Rename and retry."
            exit 110
        }
        quietly gen double W_`x' = .
        local wxvars `wxvars' W_`x'
    }

    mata: _spmixw_fill_wx_in_data("`w'", "`covars'", "`wxvars'", `N', `TT')

    _spmixw_ols `depvar' `covars' `wxvars',                ///
        id(`id') time(`time') mod_(`mod_')                 ///
        ndraw(`ndraw') nomit(`nomit') thin(`thin')         ///
        rval(`rval') nu(`nu') d0(`d0')                     ///
        bprior(`bprior') bvar(`bvar')                      ///
        saving(`"`saving'"') level(`level')                ///
        nodisp

    ereturn local model       "slx"
    ereturn local covars_orig "`covars'"
    ereturn local covars_wx   "`wxvars'"
    ereturn local wmatrix     "`w'"

    restore

    _spmixw_slx_display, level(`level')
end


program _spmixw_slx_display
    version 13.0
    syntax, [Level(cilevel)]

    tempname draws_b
    matrix `draws_b' = e(draws_b)

    di
    di as txt "Bayesian SLX panel (spmixw)" _col(58) "N obs    = " ///
        as res %9.0f e(N)
    di as txt "Model:    slx" _col(58) "N groups = " ///
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
