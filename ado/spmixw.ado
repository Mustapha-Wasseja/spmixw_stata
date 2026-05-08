*! version 0.1.0  2026-05-08  spmixw -- Bayesian spatial panel data models
*! Author: Mustapha Wasseja Mohammed <muswaseja@gmail.com>
*! Top-level dispatcher for OLS / SAR / SEM / SDM / SDEM / SLX and their
*! _conv (convex W) counterparts. See spmixw_bma for BMA over W subsets.

program spmixw, eclass
    version 13.0

    syntax varlist(min=2 numeric) [if] [in],                ///
        Model(string)                                       ///
        [                                                   ///
            W(string)                                       ///
            WMats(string)                                   ///
            EFfects(string)                                 ///
            NDraw(integer 5500)                             ///
            NOmit(integer 1500)                             ///
            Thin(integer 1)                                 ///
            Rval(real 4)                                    ///
            Nu(real 0)                                      ///
            D0(real 0)                                      ///
            BPrior(string)                                  ///
            BVar(string)                                    ///
            RMIn(real -1)                                   ///
            RMAx(real 1)                                    ///
            GRIDStep(real 0.001)                            ///
            TAYLORorder(integer 6)                          ///
            SEED(integer -1)                                ///
            SAving(string)                                  ///
            Level(cilevel)                                  ///
        ]

    // -- Validate model ------------------------------------------------------
    local model = lower("`model'")
    if !inlist("`model'", "ols", "sar", "sem", "sdm", "sdem", "slx") &       ///
       !inlist("`model'", "sar_conv", "sem_conv", "sdm_conv", "sdem_conv") {
        di as err "model(`model') not supported in v0.3.1"
        di as err "  Available standard: ols, sar, sem, sdm, sdem, slx."
        di as err "  Available convex:   sar_conv, sem_conv, sdm_conv, sdem_conv."
        di as err "  BMA / slx_conv arrive in later phases."
        exit 198
    }

    // -- Single-W models require w(), convex models require wmats() ---------
    if inlist("`model'", "sar", "sem", "sdm", "sdem", "slx") & "`w'" == "" {
        di as err "model(`model') requires the {bf:w(}{it:matname}{bf:)} option"
        exit 198
    }
    if inlist("`model'", "sar_conv", "sem_conv", "sdm_conv", "sdem_conv") &  ///
       "`wmats'" == "" {
        di as err "model(`model') requires the {bf:wmats(}{it:matnames}{bf:)} option"
        exit 198
    }
    if "`w'" != "" {
        capture confirm matrix `w'
        if _rc {
            di as err "w(`w') is not an existing Stata matrix"
            exit 198
        }
    }

    // -- Verify panel --------------------------------------------------------
    capture xtset
    if _rc {
        di as err "spmixw requires xtset; run {bf:xtset} {it:panelvar} {it:timevar} first."
        exit 459
    }
    local id   `r(panelvar)'
    local time `r(timevar)'
    if "`time'" == "" {
        di as err "spmixw requires a time variable in xtset."
        exit 459
    }

    // -- Effects -------------------------------------------------------------
    if "`effects'" == "" local effects "none"
    local effects = lower("`effects'")
    local mod_ = .
    if "`effects'" == "none"   local mod_ 0
    if "`effects'" == "region" local mod_ 1
    if "`effects'" == "time"   local mod_ 2
    if "`effects'" == "twoway" local mod_ 3
    if `mod_' == . {
        di as err "effects() must be one of: none, region, time, twoway"
        exit 198
    }

    // -- MCMC sanity ---------------------------------------------------------
    if `ndraw' <= `nomit' {
        di as err "ndraw must exceed nomit"
        exit 198
    }
    if `thin' < 1 {
        di as err "thin must be >= 1"
        exit 198
    }
    if `rval' < 0 {
        di as err "rval must be >= 0 (0 = homoscedastic)"
        exit 198
    }

    // -- Spatial grid sanity -------------------------------------------------
    if `rmin' >= `rmax' {
        di as err "rmin must be < rmax"
        exit 198
    }
    if `gridstep' <= 0 {
        di as err "gridstep must be > 0"
        exit 198
    }

    // -- Seed ----------------------------------------------------------------
    if `seed' >= 0 set seed `seed'

    // -- Touse ---------------------------------------------------------------
    marksample touse
    markout `touse' `id' `time'

    // -- y, X ----------------------------------------------------------------
    gettoken depvar covars : varlist
    if "`covars'" == "" {
        di as err "spmixw requires at least one covariate"
        exit 102
    }

    // -- Dispatch ------------------------------------------------------------
    if "`model'" == "ols" {
        _spmixw_ols `depvar' `covars' if `touse',           ///
            id(`id') time(`time') mod_(`mod_')              ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rval(`rval') nu(`nu') d0(`d0')                  ///
            bprior(`bprior') bvar(`bvar')                   ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sar" {
        _spmixw_sar `depvar' `covars' if `touse',           ///
            id(`id') time(`time') mod_(`mod_')              ///
            w(`w')                                          ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rval(`rval') nu(`nu') d0(`d0')                  ///
            bprior(`bprior') bvar(`bvar')                   ///
            rmin(`rmin') rmax(`rmax') gridstep(`gridstep')  ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sem" {
        _spmixw_sem `depvar' `covars' if `touse',           ///
            id(`id') time(`time') mod_(`mod_')              ///
            w(`w')                                          ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rval(`rval') nu(`nu') d0(`d0')                  ///
            bprior(`bprior') bvar(`bvar')                   ///
            rmin(`rmin') rmax(`rmax') gridstep(`gridstep')  ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sdm" {
        _spmixw_sdm `depvar' `covars' if `touse',           ///
            id(`id') time(`time') mod_(`mod_')              ///
            w(`w')                                          ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rval(`rval') nu(`nu') d0(`d0')                  ///
            bprior(`bprior') bvar(`bvar')                   ///
            rmin(`rmin') rmax(`rmax') gridstep(`gridstep')  ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sdem" {
        _spmixw_sdem `depvar' `covars' if `touse',          ///
            id(`id') time(`time') mod_(`mod_')              ///
            w(`w')                                          ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rval(`rval') nu(`nu') d0(`d0')                  ///
            bprior(`bprior') bvar(`bvar')                   ///
            rmin(`rmin') rmax(`rmax') gridstep(`gridstep')  ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "slx" {
        _spmixw_slx `depvar' `covars' if `touse',           ///
            id(`id') time(`time') mod_(`mod_')              ///
            w(`w')                                          ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rval(`rval') nu(`nu') d0(`d0')                  ///
            bprior(`bprior') bvar(`bvar')                   ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sar_conv" {
        _spmixw_sar_conv `depvar' `covars' if `touse',      ///
            id(`id') time(`time') mod_(`mod_')              ///
            wmats(`wmats')                                  ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rmin(`rmin') rmax(`rmax')                       ///
            taylororder(`taylororder')                      ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sem_conv" {
        _spmixw_sem_conv `depvar' `covars' if `touse',      ///
            id(`id') time(`time') mod_(`mod_')              ///
            wmats(`wmats')                                  ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rmin(`rmin') rmax(`rmax')                       ///
            taylororder(`taylororder')                      ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sdm_conv" {
        _spmixw_sdm_conv `depvar' `covars' if `touse',      ///
            id(`id') time(`time') mod_(`mod_')              ///
            wmats(`wmats')                                  ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rmin(`rmin') rmax(`rmax')                       ///
            taylororder(`taylororder')                      ///
            saving(`"`saving'"') level(`level')
    }
    else if "`model'" == "sdem_conv" {
        _spmixw_sdem_conv `depvar' `covars' if `touse',     ///
            id(`id') time(`time') mod_(`mod_')              ///
            wmats(`wmats')                                  ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')      ///
            rmin(`rmin') rmax(`rmax')                       ///
            taylororder(`taylororder')                      ///
            saving(`"`saving'"') level(`level')
    }
end
