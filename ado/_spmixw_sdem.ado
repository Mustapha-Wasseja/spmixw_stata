*! version 0.1.0  2026-05-08  spmixw — SDEM panel wrapper
*! Augments X with WX in a `preserve`d copy of the data, then routes to
*! `_spmixw_sem` (which estimates the augmented SEM). After the call, e()
*! is patched so e(model) = "sdem" and the augmented column list is
*! displayed with the SDEM header.

program _spmixw_sdem, eclass
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
            di as err "spmixw, model(sdem): cannot create variable W_`x'" ///
                " — already exists in dataset. Rename and retry."
            exit 110
        }
        quietly gen double W_`x' = .
        local wxvars `wxvars' W_`x'
    }

    mata: _spmixw_fill_wx_in_data("`w'", "`covars'", "`wxvars'", `N', `TT')

    _spmixw_sem `depvar' `covars' `wxvars',                ///
        id(`id') time(`time') mod_(`mod_')                 ///
        w(`w')                                             ///
        ndraw(`ndraw') nomit(`nomit') thin(`thin')         ///
        rval(`rval') nu(`nu') d0(`d0')                     ///
        bprior(`bprior') bvar(`bvar')                      ///
        rmin(`rmin') rmax(`rmax') gridstep(`gridstep')     ///
        saving(`"`saving'"') level(`level')                ///
        nodisp

    ereturn local model       "sdem"
    ereturn local covars_orig "`covars'"
    ereturn local covars_wx   "`wxvars'"

    restore

    _spmixw_sdem_display, level(`level')
end


program _spmixw_sdem_display
    version 13.0
    syntax, [Level(cilevel)]

    tempname draws_full
    matrix `draws_full' = e(draws_full)

    di
    di as txt "Bayesian SDEM panel (spmixw)" _col(58) "N obs    = " ///
        as res %9.0f e(N)
    di as txt "Model:    sdem" _col(58) "N groups = " ///
        as res %9.0f e(N_g)
    di as txt "Effects:  " as res "`e(effects)'" _col(58) "T        = " ///
        as res %9.0f e(T)
    di as txt "MCMC:     " as res "ndraw=`=e(ndraw)' nomit=`=e(nomit)' " ///
        "thin=`=e(thin)'  rval=`=e(rval)'"
    di as txt "sigma^2 (post. mean) = " as res %9.4f e(sige) ///
        as txt "    elapsed = " as res %6.2f e(time) as txt " s"
    di

    _spmixw_print_table, draws(`draws_full') ///
        names(`e(covars)' We) depvar(`e(depvar)') level(`level')
end
