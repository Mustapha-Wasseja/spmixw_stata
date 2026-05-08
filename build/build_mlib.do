*! version 0.1.0  2026-05-08
*! Build lspmixw.mlib from the Mata source files in mata/.
*!
*! Run from spmixw_stata/build/   :   do build_mlib.do
*! or from spmixw_stata/          :   do build/build_mlib.do
*! or with full path              :   do /path/to/build_mlib.do

version 13.0

// -- Locate package root --------------------------------------------------
// Tries (a) global SPMIXW_PKG, (b) c(filename) of this do-file when run by
// full path, then (c) cwd-relative fallbacks for "do build_mlib.do" from
// build/ or "do build/build_mlib.do" from the package root.

local PKG ""

if "${SPMIXW_PKG}" != "" {
    capture confirm file "${SPMIXW_PKG}/mata/_spmixw_demean.mata"
    if !_rc local PKG "${SPMIXW_PKG}"
}

if "`PKG'" == "" {
    local me = c(filename)
    if "`me'" != "" {
        local me_unix = subinstr("`me'", "\", "/", .)
        if regexm("`me_unix'", "^(.+)/build/build_mlib\.do$") {
            local cand = regexs(1)
            capture confirm file "`cand'/mata/_spmixw_demean.mata"
            if !_rc local PKG "`cand'"
        }
    }
}

if "`PKG'" == "" {
    capture confirm file "mata/_spmixw_demean.mata"
    if !_rc local PKG "."
}
if "`PKG'" == "" {
    capture confirm file "../mata/_spmixw_demean.mata"
    if !_rc local PKG ".."
}

if "`PKG'" == "" {
    di as err "build_mlib.do: cannot locate spmixw_stata/mata/."
    di as err "  Set global SPMIXW_PKG or cd to the package root before running."
    exit 601
}

local MATA "`PKG'/mata"
// Drop the .mlib next to the ados so a single `adopath ++ "ado/"` resolves
// both. A copy is left in build/ for archival.
local OUT  "`PKG'/ado"
local ARCH "`PKG'/build"

// -- Source files (order matters: callees must be compiled before callers) -
local FILES                          ///
    _spmixw_demean.mata              ///
    _spmixw_ols.mata                 ///
    _spmixw_simulate.mata            ///
    _spmixw_logdet_exact.mata        ///
    _spmixw_griddy_rho_sar.mata      ///
    _spmixw_sar.mata                 ///
    _spmixw_simulate_sar.mata        ///
    _spmixw_summary.mata             ///
    _spmixw_linear_interp.mata       ///
    _spmixw_griddy_rho_sem.mata      ///
    _spmixw_sem.mata                 ///
    _spmixw_simulate_sem.mata        ///
    _spmixw_compute_wx.mata          ///
    _spmixw_trace_mc.mata            ///
    _spmixw_effects.mata             ///
    _spmixw_logdet_taylor.mata       ///
    _spmixw_sar_conv.mata            ///
    _spmixw_sar_conv_caller.mata     ///
    _spmixw_simulate_sar_conv.mata   ///
    _spmixw_simulate_sem_conv.mata   ///
    _spmixw_model_probs.mata         ///
    _spmixw_effects_conv.mata

// -- Compile ---------------------------------------------------------------
di as txt "Building lspmixw.mlib from `MATA'/"

// Some Stata 13 builds mis-parse `mata clear` as `mata` + `clear` (entering
// Mata mode and evaluating clear as an expression). Use explicit per-function
// drops instead -- same defensive pattern as tests/_setup.do.
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
capture mata: mata drop _spmixw_sar_conv_caller()
capture mata: mata drop _spmixw_simulate_sar_conv()
capture mata: mata drop _spmixw_simulate_sem_conv()
capture mata: mata drop _spmixw_model_probs()
capture mata: mata drop _spmixw_trace_cross_mc()
capture mata: mata drop _spmixw_effects_sdm_conv()
capture mata: mata drop _spmixw_effects_sdem_conv()

foreach f of local FILES {
    di as txt "  compiling `f' ..."
    quietly do "`MATA'/`f'"
}

di as txt "Packaging into `OUT'/lspmixw.mlib ..."

// Some Stata 13 builds reject the Stata-level form `mata mlib create/add`
// outright (parser returns "invalid expression r(3000)"). Run the same
// commands from inside a `mata: ... end` block instead -- inside Mata,
// `mlib create` / `mlib add` are recognised subcommands and the flaky
// outer parser is bypassed entirely.
local olddir = c(pwd)
cd "`OUT'"

mata:
    mlib create lspmixw, replace
    mlib add lspmixw _spmixw_demean()
    mlib add lspmixw _spmixw_ols()
    mlib add lspmixw _spmixw_simulate_ols()
    mlib add lspmixw _spmixw_logdet_exact()
    mlib add lspmixw _spmixw_griddy_rho_sar()
    mlib add lspmixw _spmixw_sar()
    mlib add lspmixw _spmixw_simulate_sar()
    mlib add lspmixw _spmixw_summary()
    mlib add lspmixw _spmixw_linear_interp()
    mlib add lspmixw _spmixw_griddy_rho_sem()
    mlib add lspmixw _spmixw_sem()
    mlib add lspmixw _spmixw_simulate_sem()
    mlib add lspmixw _spmixw_compute_wx()
    mlib add lspmixw _spmixw_fill_wx_in_data()
    mlib add lspmixw _spmixw_trace_mc()
    mlib add lspmixw _spmixw_effects_sar()
    mlib add lspmixw _spmixw_effects_sdm()
    mlib add lspmixw _spmixw_effects_simple()
    mlib add lspmixw _spmixw_logdet_taylor()
    mlib add lspmixw _spmixw_eval_taylor_lndet()
    mlib add lspmixw _spmixw_gamma_proposal_uniform()
    mlib add lspmixw _spmixw_gamma_proposal_adapted()
    mlib add lspmixw _spmixw_eval_cond_sar_conv()
    mlib add lspmixw _spmixw_sar_conv()
    mlib add lspmixw _spmixw_sar_conv_caller()
    mlib add lspmixw _spmixw_simulate_sar_conv()
    mlib add lspmixw _spmixw_simulate_sem_conv()
    mlib add lspmixw _spmixw_model_probs()
    mlib add lspmixw _spmixw_trace_cross_mc()
    mlib add lspmixw _spmixw_effects_sdm_conv()
    mlib add lspmixw _spmixw_effects_sdem_conv()
    mlib index
end

cd "`olddir'"

// -- Verify the .mlib actually exists -----------------------------------
capture confirm file "`OUT'/lspmixw.mlib"
if _rc {
    di as err "build_mlib.do: lspmixw.mlib was NOT created at `OUT'/."
    di as err "  inspect any Mata errors above."
    exit 601
}

// Archival copy
capture mkdir "`ARCH'"
quietly copy "`OUT'/lspmixw.mlib" "`ARCH'/lspmixw.mlib", replace

di as txt "{hline 60}"
di as txt "BUILD OK: lspmixw.mlib written to:"
di as res "  `OUT'/lspmixw.mlib"
di as txt "  (archival copy: `ARCH'/lspmixw.mlib)"
di as txt "{hline 60}"
