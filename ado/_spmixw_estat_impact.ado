*! version 0.1.0  2026-05-08  spmixw - estat impact (LeSage-Pace effects)
*!
*! Computes direct, indirect, and total spatial effects for:
*!   sar / sdm  : stochastic-trace tr(W^j)/N, with W = e(wmatrix)
*!   sdem / slx : exact identity direct=beta, indirect=theta, total=beta+theta
*!   sar_conv  : same SAR formula with W replaced by W_c(gamma_pm)
*!   sdm_conv  : SDM formula with W_c, plus per-W cross-traces tr(W_c^r W_m)/N
*!   sdem_conv : trivial direct=beta, indirect=sum_m theta_m, total=beta+indir

program _spmixw_estat_impact, rclass
    version 13.0

    syntax [, Level(cilevel) ITER(integer 50) ORDer(integer 100) ///
              SEED(integer -1)]

    local model `e(model)'

    local is_std  = inlist("`model'", "sar", "sdm", "sdem", "slx")
    local is_conv = inlist("`model'", "sar_conv", "sdm_conv", "sdem_conv")

    if !`is_std' & !`is_conv' {
        di as err "estat impact: not applicable for model(`model')"
        di as err "  Effects available for: sar, sdm, sdem, slx, " ///
                  "sar_conv, sdm_conv, sdem_conv."
        exit 198
    }

    // ------------------------------------------------------------------
    // Determine `p` (number of original X variables) and the effect names.
    //   SAR        : p = e(k); names = e(covars)
    //   SDM/SDEM/SLX (std) : p = #(covars_orig); names = covars_orig
    //   SAR conv   : p = e(k); names = e(covars)  (no augmentation)
    //   SDM/SDEM conv      : p = #(covars_orig); names = covars_orig
    // ------------------------------------------------------------------
    local effvars
    local p
    if `is_std' {
        if inlist("`model'", "sdm", "sdem", "slx") {
            local effvars `e(covars_orig)'
            local p : word count `effvars'
        }
        else {
            local effvars `e(covars)'
            local p = e(k)
        }
    }
    else {
        // Convex models
        if "`model'" == "sar_conv" {
            local effvars `e(covars)'
            local p = e(k)
        }
        else {
            local effvars `e(covars_orig)'
            local p : word count `effvars'
        }
    }

    if `p' < 1 {
        di as err "estat impact: cannot determine number of variables"
        exit 198
    }

    tempname direct_t indirect_t total_t effects_full

    // Stata's `set seed` shares state with Mata's RNG.
    if `seed' >= 0 set seed `seed'

    if `is_std' {
        // ------------------------------------------------------------------
        // Standard models: existing caller, single W in e(wmatrix).
        // ------------------------------------------------------------------
        local W `e(wmatrix)'
        local has_pdraw = inlist("`model'", "sar", "sdm")

        mata: _spmixw_estat_impact_caller(             ///
            "`model'", "`W'", `has_pdraw',             ///
            `p', `iter', `order',                      ///
            "`effects_full'",                          ///
            "`direct_t'", "`indirect_t'", "`total_t'", ///
            `level' / 100)
    }
    else {
        // ------------------------------------------------------------------
        // Convex models: read wmats list and gamma posterior mean, build
        // W_c, compute tracew and (for SDM conv) per-W cross-traces.
        // ------------------------------------------------------------------
        local wmats `e(wmats)'
        local M = e(M)
        if `M' < 1 {
            di as err "estat impact: e(M) is missing or zero"
            exit 198
        }

        mata: _spmixw_estat_impact_conv_caller(        ///
            "`model'", "`wmats'", `M',                 ///
            `p', `iter', `order',                      ///
            "`effects_full'",                          ///
            "`direct_t'", "`indirect_t'", "`total_t'", ///
            `level' / 100)
    }

    matrix rownames `direct_t'   = `effvars'
    matrix rownames `indirect_t' = `effvars'
    matrix rownames `total_t'    = `effvars'
    matrix colnames `direct_t'   = mean sd lo hi
    matrix colnames `indirect_t' = mean sd lo hi
    matrix colnames `total_t'    = mean sd lo hi

    // ------------------------------------------------------------------
    // Display
    // ------------------------------------------------------------------
    di
    di as txt "Spatial effects (LeSage-Pace scalar summary, " ///
        as txt "level=`level'%)"
    di as txt "  model = " as res "`model'" ///
        as txt "    iter = " as res `iter' ///
        as txt "    order = " as res `order'

    di
    di as txt "Direct effects:"
    _spmixw_print_effects_subtable, table(`direct_t')   names(`effvars') level(`level')

    di
    di as txt "Indirect effects:"
    _spmixw_print_effects_subtable, table(`indirect_t') names(`effvars') level(`level')

    di
    di as txt "Total effects:"
    _spmixw_print_effects_subtable, table(`total_t')    names(`effvars') level(`level')

    return matrix direct   = `direct_t'
    return matrix indirect = `indirect_t'
    return matrix total    = `total_t'
    return scalar p        = `p'
    return local model     "`model'"

    capture confirm matrix `effects_full'
    if !_rc {
        return matrix effects_draws = `effects_full'
    }
end


// ----------------------------------------------------------------------
// Mata bridge - standard (non-convex) models.
// ----------------------------------------------------------------------
version 13.0
mata:
mata set matastrict on

void _spmixw_estat_impact_caller(string scalar model_name,
                                 string scalar W_name,
                                 real   scalar has_pdraw,
                                 real   scalar p,
                                 real   scalar uiter,
                                 real   scalar maxorder,
                                 string scalar effects_full_name,
                                 string scalar direct_t_name,
                                 string scalar indirect_t_name,
                                 string scalar total_t_name,
                                 real   scalar lvl)
{
    real matrix    bdraw, W, all_draws, direct_d, indirect_d, total_d
    real colvector pdraw, tracew

    bdraw = st_matrix("e(draws_b)")
    if (has_pdraw) {
        pdraw = st_matrix("e(draws_p)")[., 1]
    }
    else {
        pdraw = J(rows(bdraw), 1, 0)
    }

    if (model_name == "sar") {
        W         = st_matrix(W_name)
        tracew    = _spmixw_trace_mc(W, maxorder, uiter)
        all_draws = _spmixw_effects_sar(bdraw, pdraw, tracew, p)
    }
    else if (model_name == "sdm") {
        W         = st_matrix(W_name)
        tracew    = _spmixw_trace_mc(W, maxorder, uiter)
        all_draws = _spmixw_effects_sdm(bdraw, pdraw, tracew, p)
    }
    else if (model_name == "sdem" | model_name == "slx") {
        all_draws = _spmixw_effects_simple(bdraw, p)
    }
    else {
        _error(3300, "_spmixw_estat_impact_caller: unsupported model " + model_name)
    }

    direct_d   = all_draws[., (1::p)]
    indirect_d = all_draws[., ((p+1)::(2*p))]
    total_d    = all_draws[., ((2*p+1)::(3*p))]

    st_matrix(direct_t_name,   _spmixw_summary(direct_d,   lvl))
    st_matrix(indirect_t_name, _spmixw_summary(indirect_d, lvl))
    st_matrix(total_t_name,    _spmixw_summary(total_d,    lvl))

    if (max((rows(all_draws), cols(all_draws))) <= st_numscalar("c(matsize)")) {
        st_matrix(effects_full_name, all_draws)
    }
}


// ----------------------------------------------------------------------
// Mata bridge - convex models (sar_conv, sdm_conv, sdem_conv).
// Reads e(wmats) (space-sep list of W matrix names), e(gam) (1 x M
// posterior mean weights), and e(draws_b) / e(draws_p) for the per-draw
// betas / rhos. Builds W_c at posterior mean gamma and routes to the
// matching effects function.
// ----------------------------------------------------------------------
void _spmixw_estat_impact_conv_caller(string scalar model_name,
                                      string scalar wmats_list,
                                      real   scalar M,
                                      real   scalar p,
                                      real   scalar uiter,
                                      real   scalar maxorder,
                                      string scalar effects_full_name,
                                      string scalar direct_t_name,
                                      string scalar indirect_t_name,
                                      string scalar total_t_name,
                                      real   scalar lvl)
{
    real matrix    bdraw, W_c, all_draws, direct_d, indirect_d, total_d
    real matrix    tracew_cm, W_m
    real colvector pdraw, tracew, gam_pm
    real scalar    N, m
    string rowvector W_names

    bdraw  = st_matrix("e(draws_b)")
    gam_pm = st_matrix("e(gam)")[1, .]'      // M x 1
    if (rows(gam_pm) != M) {
        _error(3300, "_spmixw_estat_impact_conv_caller: e(gam) length != M")
    }

    W_names = tokens(wmats_list)
    if (cols(W_names) != M) {
        _error(3300, "_spmixw_estat_impact_conv_caller: wmats count != M")
    }

    // Build W_c at posterior-mean gamma.
    N   = rows(st_matrix(W_names[1]))
    W_c = J(N, N, 0)
    for (m = 1; m <= M; m++) {
        W_c = W_c + gam_pm[m] :* st_matrix(W_names[m])
    }

    if (model_name == "sar_conv") {
        pdraw     = st_matrix("e(draws_p)")[., 1]
        tracew    = _spmixw_trace_mc(W_c, maxorder, uiter)
        all_draws = _spmixw_effects_sar(bdraw, pdraw, tracew, p)
    }
    else if (model_name == "sdm_conv") {
        pdraw  = st_matrix("e(draws_p)")[., 1]
        tracew = _spmixw_trace_mc(W_c, maxorder, uiter)

        // Cross-traces: tracew_cm[r+1, m] = tr(W_c^r W_m) / N for r = 0..maxorder.
        tracew_cm = J(maxorder + 1, M, 0)
        for (m = 1; m <= M; m++) {
            W_m = st_matrix(W_names[m])
            tracew_cm[., m] = _spmixw_trace_cross_mc(W_c, W_m, maxorder, uiter)
        }
        all_draws = _spmixw_effects_sdm_conv(bdraw, pdraw, tracew, tracew_cm, p, M)
    }
    else if (model_name == "sdem_conv") {
        all_draws = _spmixw_effects_sdem_conv(bdraw, p, M)
    }
    else {
        _error(3300, "_spmixw_estat_impact_conv_caller: unsupported model " + model_name)
    }

    direct_d   = all_draws[., (1::p)]
    indirect_d = all_draws[., ((p+1)::(2*p))]
    total_d    = all_draws[., ((2*p+1)::(3*p))]

    st_matrix(direct_t_name,   _spmixw_summary(direct_d,   lvl))
    st_matrix(indirect_t_name, _spmixw_summary(indirect_d, lvl))
    st_matrix(total_t_name,    _spmixw_summary(total_d,    lvl))

    if (max((rows(all_draws), cols(all_draws))) <= st_numscalar("c(matsize)")) {
        st_matrix(effects_full_name, all_draws)
    }
}

end


// ----------------------------------------------------------------------
// Pretty-printer for an effects sub-table.
// `table` is a (p x 4) matrix with columns [mean, sd, lo, hi].
// ----------------------------------------------------------------------
program _spmixw_print_effects_subtable
    version 13.0
    syntax , TABLE(name) NAMES(string) [Level(cilevel)]

    local k = rowsof(`table')
    local nwords : word count `names'
    if `k' != `nwords' {
        di as err "_spmixw_print_effects_subtable: row count mismatch"
        exit 198
    }

    di as txt "{hline 78}"
    di as txt %14s "Variable" "  |" ///
        as txt "  Post. Mean   Post. SD   [`level'% Cred. Interval]"
    di as txt "{hline 14}+{hline 63}"

    forvalues j = 1/`k' {
        local nm : word `j' of `names'
        local m  = `table'[`j', 1]
        local sd = `table'[`j', 2]
        local lo = `table'[`j', 3]
        local hi = `table'[`j', 4]
        di as txt %14s "`nm'" "  |" ///
            as res %12.4f `m' %11.4f `sd' ///
            "    " %9.4f `lo' "  " %9.4f `hi'
    }
    di as txt "{hline 78}"
end
