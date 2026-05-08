*! version 0.1.0  2026-05-08  spmixw_bma -- BMA over W subsets, all four convex models
*!
*! Enumerates all 2^M - 1 non-empty subsets of M weight matrices, fits the
*! corresponding convex-W estimator (SAR / SEM / SDM / SDEM convex) on each
*! subset, and BMA-averages point estimates by posterior model probabilities
*! (softmax of log-marginal likelihoods).
*!
*! Reference: Debarsy, N. and LeSage, J. P. (2021). "Bayesian model averaging
*! for spatial autoregressive models based on convex combinations of
*! different types of connectivity matrices." JBES 40(2), 547-558.

program spmixw_bma, eclass
    version 13.0

    syntax varlist(min=2 numeric) [if] [in],                ///
        Model(string)                                       ///
        WMats(string)                                       ///
        [                                                   ///
            EFfects(string)                                 ///
            NDraw(integer 15000)                            ///
            NOmit(integer 5000)                             ///
            Thin(integer 1)                                 ///
            RMIn(real -0.99)                                ///
            RMAx(real 0.99)                                 ///
            TAYLORorder(integer 6)                          ///
            SEED(integer -1)                                ///
            Level(cilevel)                                  ///
        ]

    // -- Validate model ----------------------------------------------------
    local model = lower("`model'")
    if !inlist("`model'", "sar_conv", "sem_conv", "sdm_conv", "sdem_conv") {
        di as err "model() must be one of: sar_conv, sem_conv, sdm_conv, sdem_conv"
        exit 198
    }
    local has_wx = inlist("`model'", "sdm_conv", "sdem_conv")
    local subprog "_spmixw_`model'"

    // Display-friendly labels (rho vs lambda, Wy vs We)
    local sp_label "rho"
    local row_label "Wy"
    if inlist("`model'", "sem_conv", "sdem_conv") {
        local sp_label  "lambda"
        local row_label "We"
    }

    // -- Parse wmats list --------------------------------------------------
    local M : word count `wmats'
    if `M' < 2 {
        di as err "BMA requires at least 2 W matrices in wmats()"
        exit 198
    }
    foreach Wm of local wmats {
        capture confirm matrix `Wm'
        if _rc {
            di as err "wmats: '`Wm'' is not an existing Stata matrix"
            exit 198
        }
    }

    // -- Verify panel ------------------------------------------------------
    capture xtset
    if _rc {
        di as err "spmixw_bma requires xtset; run {bf:xtset} {it:panelvar} {it:timevar} first."
        exit 459
    }
    local id   `r(panelvar)'
    local time `r(timevar)'

    // -- FE encoding -------------------------------------------------------
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

    // -- Touse, depvar/covars ---------------------------------------------
    marksample touse
    markout `touse' `id' `time'

    gettoken depvar covars : varlist
    local k : word count `covars'

    if `seed' >= 0 set seed `seed'

    // -- Enumerate non-empty subsets ---------------------------------------
    local n_models = 2^`M' - 1
    tempname subsets
    matrix `subsets' = J(`n_models', `M', 0)

    forvalues i = 1/`n_models' {
        local mask = `i'
        forvalues m = 1/`M' {
            local bit = mod(`mask', 2)
            matrix `subsets'[`i', `m'] = `bit'
            local mask = floor(`mask' / 2)
        }
    }

    // -- Per-model storage -------------------------------------------------
    tempname logm beta_all rho_all sige_all gam_all
    matrix `logm'     = J(`n_models', 1,     0)
    matrix `beta_all' = J(`n_models', `k',   0)
    matrix `rho_all'  = J(`n_models', 1,     0)
    matrix `sige_all' = J(`n_models', 1,     0)
    matrix `gam_all'  = J(`n_models', `M',   0)

    // For SDM/SDEM: theta_all stores W*X coefficients on the full (M x k) grid.
    // Column index for (W_m, x_j) is (m-1)*k + j. Models that do not include
    // W_m have zeros in those columns (and contribute zero posterior weight to
    // those (m, j) entries via probs * theta_all).
    if `has_wx' {
        tempname theta_all
        matrix `theta_all' = J(`n_models', `=`M' * `k'', 0)
    }

    di as txt _newline "{hline 70}"
    di as txt "  spmixw_bma: " upper("`model'") " BMA over `n_models' W subsets"
    di as txt "{hline 70}"

    timer clear 2
    timer on 2

    forvalues i = 1/`n_models' {
        // Build subset wmats list and remember member positions
        local sub_wmats
        local sub_positions   // m index (1..M) of each member, in order
        forvalues m = 1/`M' {
            if `subsets'[`i', `m'] == 1 {
                local sub_wmats `sub_wmats' `: word `m' of `wmats''
                local sub_positions `sub_positions' `m'
            }
        }
        local nw : word count `sub_wmats'

        di as txt "  model `i'/`n_models': W = {" as res "`sub_wmats'" ///
            as txt "} (`nw' matrices)" _continue

        quietly `subprog' `varlist' if `touse',           ///
            id(`id') time(`time') mod_(`mod_')             ///
            wmats(`sub_wmats')                             ///
            ndraw(`ndraw') nomit(`nomit') thin(`thin')     ///
            rmin(`rmin') rmax(`rmax')                      ///
            taylororder(`taylororder')                     ///
            saving("") level(`level')                      ///
            nodisp

        di as txt "  logm = " as res %9.3f e(logmarg)

        matrix `logm'[`i', 1]     = e(logmarg)
        matrix `rho_all'[`i', 1]  = e(rho)
        matrix `sige_all'[`i', 1] = e(sige)

        // β_orig: e(b) cols 1..k. (For SDM/SDEM, cols k+1..k+k*nw are W*X
        // blocks, last col is rho/lambda. For SAR/SEM, last col is rho/lambda.)
        forvalues j = 1/`k' {
            matrix `beta_all'[`i', `j'] = el(e(b), 1, `j')
        }

        // For SDM/SDEM: extract the W*X block coefficients and map them back
        // to the full (M, k) grid using sub_positions.
        if `has_wx' {
            local mp = 0
            foreach m_abs of local sub_positions {
                local mp = `mp' + 1
                forvalues j = 1/`k' {
                    local src_col = `k' + (`mp' - 1) * `k' + `j'
                    local dst_col = (`m_abs' - 1) * `k' + `j'
                    matrix `theta_all'[`i', `dst_col'] = el(e(b), 1, `src_col')
                }
            }
        }

        // γ: e(gam) is 1 × nw, ordered by sub_positions.
        local idx = 0
        forvalues m = 1/`M' {
            if `subsets'[`i', `m'] == 1 {
                local idx = `idx' + 1
                matrix `gam_all'[`i', `m'] = el(e(gam), 1, `idx')
            }
        }
    }

    timer off 2
    quietly timer list 2
    local elapsed = r(t2)
    timer clear 2

    // -- Posterior model probabilities -------------------------------------
    tempname probs bma_beta bma_gamma
    matrix `probs' = J(`n_models', 1, 0)
    mata: st_matrix("`probs'", _spmixw_model_probs(st_matrix("`logm'")))

    // -- BMA-weighted point estimates --------------------------------------
    matrix `bma_beta'  = `probs'' * `beta_all'      // 1 × k
    matrix `bma_gamma' = `probs'' * `gam_all'       // 1 × M

    if `has_wx' {
        tempname bma_theta
        matrix `bma_theta' = `probs'' * `theta_all'  // 1 × (M*k)
    }

    tempname bma_rho_m bma_sige_m
    matrix `bma_rho_m'  = `probs'' * `rho_all'      // 1 × 1
    matrix `bma_sige_m' = `probs'' * `sige_all'     // 1 × 1
    local bma_rho  = el(`bma_rho_m',  1, 1)
    local bma_sige = el(`bma_sige_m', 1, 1)

    // -- Display -----------------------------------------------------------
    di as txt _newline "{hline 70}"
    di as txt "  BMA results"
    di as txt "{hline 70}"
    di as txt "  Model = `model',  M = `M' weight matrices, " ///
        "n_models = `n_models'"
    di as txt "  elapsed = " as res %5.1f `elapsed' as txt " s"
    di
    di as txt "Posterior model probabilities:"
    di as txt "{hline 70}"
    di as txt "  Model" _col(12) "W subset" _col(40) "logmarg" _col(55) "p(M|y)"
    di as txt "{hline 70}"
    forvalues i = 1/`n_models' {
        local lbl
        forvalues m = 1/`M' {
            if `subsets'[`i', `m'] == 1 local lbl `lbl' `: word `m' of `wmats''
        }
        di as txt %5.0f `i' _col(12) "`lbl'" _col(36) ///
            as res %12.4f el(`logm', `i', 1) ///
            "  " %8.4f el(`probs', `i', 1)
    }
    di as txt "{hline 70}"

    di
    di as txt "BMA-weighted point estimates:"
    di as txt "  `sp_label'" _col(12) "= " as res %8.4f `bma_rho'
    di as txt "  sigma^2" _col(12) "= " as res %8.4f `bma_sige'
    di as txt "  beta:"
    forvalues j = 1/`k' {
        local nm : word `j' of `covars'
        di as txt "    " %-15s "`nm'" "= " as res %8.4f el(`bma_beta', 1, `j')
    }
    di as txt "  gamma:"
    forvalues m = 1/`M' {
        local nm : word `m' of `wmats'
        di as txt "    " %-15s "`nm'" "= " as res %8.4f el(`bma_gamma', 1, `m')
    }
    if `has_wx' {
        di as txt "  W*X coefficients (zeros = W not in any model with mass):"
        forvalues m = 1/`M' {
            local Wnm : word `m' of `wmats'
            forvalues j = 1/`k' {
                local xnm : word `j' of `covars'
                local col = (`m' - 1) * `k' + `j'
                di as txt "    " %-15s "`Wnm'*`xnm'" "= " ///
                    as res %8.4f el(`bma_theta', 1, `col')
            }
        }
    }

    // -- ereturn -----------------------------------------------------------
    ereturn clear
    ereturn local cmd       "spmixw_bma"
    ereturn local model     "`model'"
    ereturn local depvar    "`depvar'"
    ereturn local covars    "`covars'"
    ereturn local wmats     "`wmats'"
    ereturn local effects   "`effects'"

    ereturn scalar M        = `M'
    ereturn scalar nmodels  = `n_models'
    ereturn scalar bma_rho  = `bma_rho'
    ereturn scalar bma_sige = `bma_sige'
    ereturn scalar time     = `elapsed'

    matrix colnames `bma_beta'  = `covars'
    matrix colnames `bma_gamma' = `wmats'
    matrix colnames `gam_all'   = `wmats'
    matrix colnames `subsets'   = `wmats'

    ereturn matrix bma_beta  = `bma_beta'
    ereturn matrix bma_gamma = `bma_gamma'
    ereturn matrix probs     = `probs'
    ereturn matrix logm      = `logm'
    ereturn matrix subsets   = `subsets'
    ereturn matrix beta_all  = `beta_all'
    ereturn matrix rho_all   = `rho_all'
    ereturn matrix sige_all  = `sige_all'
    ereturn matrix gam_all   = `gam_all'

    if `has_wx' {
        // Build (M*k)-length name list: W1_x1 W1_x2 ... W1_xk W2_x1 ...
        local theta_names
        forvalues m = 1/`M' {
            local Wnm : word `m' of `wmats'
            forvalues j = 1/`k' {
                local xnm : word `j' of `covars'
                local theta_names `theta_names' `Wnm'_`xnm'
            }
        }
        matrix colnames `bma_theta'  = `theta_names'
        matrix colnames `theta_all'  = `theta_names'
        ereturn matrix bma_theta  = `bma_theta'
        ereturn matrix theta_all  = `theta_all'
    }
end
