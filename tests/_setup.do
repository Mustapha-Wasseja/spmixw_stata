*! Sourced by every test/validation script.
*! Adds spmixw_stata/ado to the adopath so spmixw and lspmixw.mlib resolve.

version 13.0

// Locate package root from this script's directory.
// Stata's `c(filename)` after `do _setup` won't help, but tests use a known
// relative path. Tests/validation files set $SPMIXW_PKG before sourcing.

if "${SPMIXW_PKG}" == "" {
    di as err "_setup.do: please set global SPMIXW_PKG to spmixw_stata/ before sourcing"
    exit 198
}

capture confirm file "${SPMIXW_PKG}/ado/spmixw.ado"
if _rc {
    di as err "spmixw.ado not found at \${SPMIXW_PKG}/ado/spmixw.ado"
    di as err "  current value: ${SPMIXW_PKG}"
    exit 601
}

capture confirm file "${SPMIXW_PKG}/ado/lspmixw.mlib"
if _rc {
    di as txt "Note: lspmixw.mlib not yet built. Run build/build_mlib.do first."
}

quietly adopath ++ "${SPMIXW_PKG}/ado"

// Refresh Mata library index (no-op if no .mlib present yet).
capture mata mlib index

// Defensive: source each .mata directly into the session. This guarantees
// the functions exist regardless of whether lspmixw.mlib was built or is
// resolvable on the current adopath. The .mlib is the production path; this
// is the dev / test path.
//
// `mata: void f() {...}` errors with r(3000) if f already exists in this
// session, so drop first (capture-wrapped because the first run won't have
// the functions defined yet).
capture mata: mata drop _spmixw_demean()
capture mata: mata drop _spmixw_ols()
capture mata: mata drop _spmixw_simulate_ols()
capture mata: mata drop _spmixw_logdet_exact()
capture mata: mata drop _spmixw_griddy_rho_sar()
capture mata: mata drop _spmixw_sar()
capture mata: mata drop _spmixw_simulate_sar()
capture mata: mata drop _spmixw_summary()
capture mata: mata drop _spmixw_linear_interp()
capture mata: mata drop _spmixw_griddy_rho_sem()
capture mata: mata drop _spmixw_sem()
capture mata: mata drop _spmixw_simulate_sem()
capture mata: mata drop _spmixw_compute_wx()
capture mata: mata drop _spmixw_fill_wx_in_data()
capture mata: mata drop _spmixw_trace_mc()
capture mata: mata drop _spmixw_effects_sar()
capture mata: mata drop _spmixw_effects_sdm()
capture mata: mata drop _spmixw_effects_simple()
capture mata: mata drop _spmixw_logdet_taylor()
capture mata: mata drop _spmixw_eval_taylor_lndet()
capture mata: mata drop _spmixw_gamma_proposal_uniform()
capture mata: mata drop _spmixw_gamma_proposal_adapted()
capture mata: mata drop _spmixw_eval_cond_sar_conv()
capture mata: mata drop _spmixw_sar_conv()
capture mata: mata drop _spmixw_simulate_sar_conv()
capture mata: mata drop _spmixw_simulate_sem_conv()
capture mata: mata drop _spmixw_sar_conv_caller()
capture mata: mata drop _spmixw_model_probs()
capture mata: mata drop _spmixw_trace_cross_mc()
capture mata: mata drop _spmixw_effects_sdm_conv()
capture mata: mata drop _spmixw_effects_sdem_conv()

quietly do "${SPMIXW_PKG}/mata/_spmixw_demean.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_ols.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_simulate.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_logdet_exact.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_griddy_rho_sar.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_sar.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_simulate_sar.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_summary.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_linear_interp.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_griddy_rho_sem.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_sem.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_simulate_sem.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_compute_wx.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_trace_mc.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_effects.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_logdet_taylor.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_sar_conv.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_simulate_sar_conv.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_simulate_sem_conv.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_sar_conv_caller.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_model_probs.mata"
quietly do "${SPMIXW_PKG}/mata/_spmixw_effects_conv.mata"

di as txt "spmixw test setup OK (adopath: ${SPMIXW_PKG}/ado)"
