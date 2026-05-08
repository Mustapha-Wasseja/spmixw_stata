*! version 0.1.0  2026-05-08  spmixw — SDM panel wrapper
*! Augments X with WX in a `preserve`d copy of the data, then routes to
*! `_spmixw_sar` (which estimates the augmented SAR). After the call, e()
*! is patched so e(model) = "sdm" and the augmented column list is
*! displayed with the SDM header. The added W_<x> variables exist only
*! within the preserve scope.

program _spmixw_sdm, eclass
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

    // -- Establish panel dimensions on the active sample ---------------------
    tempvar tagid tagt
    quietly egen `tagid' = tag(`id')   if `touse'
    quietly egen `tagt'  = tag(`time') if `touse'
    quietly count if `tagid' == 1
    local N = r(N)
    quietly count if `tagt' == 1
    local TT = r(N)

    // -- Augment X with WX inside a preserve scope ---------------------------
    preserve
    quietly keep if `touse'
    sort `time' `id'

    local k : word count `covars'
    local wxvars
    foreach x of local covars {
        capture confirm new variable W_`x'
        if _rc {
            di as err "spmixw, model(sdm): cannot create variable W_`x'" ///
                " — already exists in dataset. Rename and retry."
            exit 110
        }
        quietly gen double W_`x' = .
        local wxvars `wxvars' W_`x'
    }

    // Fill W_x* via Mata
    mata: _spmixw_fill_wx_in_data("`w'", "`covars'", "`wxvars'", `N', `TT')

    // -- Call SAR with augmented X (suppress its display) -------------------
    _spmixw_sar `depvar' `covars' `wxvars',                            ///
        id(`id') time(`time') mod_(`mod_')                             ///
        w(`w')                                                         ///
        ndraw(`ndraw') nomit(`nomit') thin(`thin')                     ///
        rval(`rval') nu(`nu') d0(`d0')                                 ///
        bprior(`bprior') bvar(`bvar')                                  ///
        rmin(`rmin') rmax(`rmax') gridstep(`gridstep')                 ///
        saving(`"`saving'"') level(`level')                            ///
        nodisp

    // -- Patch metadata: model label, original-vs-augmented covar lists ----
    ereturn local model       "sdm"
    ereturn local covars_orig "`covars'"
    ereturn local covars_wx   "`wxvars'"

    restore

    // -- Display with SDM header -------------------------------------------
    _spmixw_sdm_display, level(`level')
end


// ----------------------------------------------------------------------
// Display helper — same Bayesian table as SAR, with header "SDM panel".
// The coefficient table includes original X, W_<X>, and Wy rows.
// ----------------------------------------------------------------------
program _spmixw_sdm_display
    version 13.0
    syntax, [Level(cilevel)]

    tempname draws_full
    matrix `draws_full' = e(draws_full)

    di
    di as txt "Bayesian SDM panel (spmixw)" _col(58) "N obs    = " ///
        as res %9.0f e(N)
    di as txt "Model:    sdm" _col(58) "N groups = " ///
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
