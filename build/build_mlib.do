*! version 0.1.0  2026-05-07
*! Build lspmixw.mlib from spmixw_stata/mata/*.mata
*!
*! Run from spmixw_stata/build/   :   do build_mlib.do
*! or from spmixw_stata/          :   do build/build_mlib.do
*! or with full path              :   do /path/to/build_mlib.do

version 13.0

// -- Locate package root --------------------------------------------------
// Works whether this do-file is run from build/, from package root, or by
// full path: try the cwd's mata/ first, fall back to ../mata/.

local PKG ""
capture confirm file "mata/_spmixw_demean.mata"
if !_rc {
    local PKG "."
}
else {
    capture confirm file "../mata/_spmixw_demean.mata"
    if !_rc {
        local PKG ".."
    }
}
if "`PKG'" == "" {
    di as err "build_mlib.do: cannot locate spmixw_stata/mata/. " ///
        "Run from spmixw_stata/ or spmixw_stata/build/."
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

mata clear

foreach f of local FILES {
    di as txt "  compiling `f' ..."
    quietly do "`MATA'/`f'"
}

di as txt "Packaging into `OUT'/lspmixw.mlib ..."

// `mata mlib *` are Stata commands (not Mata code) — no `mata:` prefix.
mata mlib create lspmixw, dir("`OUT'") replace
mata mlib add lspmixw _spmixw_demean()           , dir("`OUT'")
mata mlib add lspmixw _spmixw_ols()              , dir("`OUT'")
mata mlib add lspmixw _spmixw_simulate_ols()     , dir("`OUT'")
mata mlib add lspmixw _spmixw_logdet_exact()     , dir("`OUT'")
mata mlib add lspmixw _spmixw_griddy_rho_sar()   , dir("`OUT'")
mata mlib add lspmixw _spmixw_sar()              , dir("`OUT'")
mata mlib add lspmixw _spmixw_simulate_sar()     , dir("`OUT'")
mata mlib add lspmixw _spmixw_summary()          , dir("`OUT'")
mata mlib add lspmixw _spmixw_linear_interp()    , dir("`OUT'")
mata mlib add lspmixw _spmixw_griddy_rho_sem()   , dir("`OUT'")
mata mlib add lspmixw _spmixw_sem()              , dir("`OUT'")
mata mlib add lspmixw _spmixw_simulate_sem()     , dir("`OUT'")
mata mlib add lspmixw _spmixw_compute_wx()       , dir("`OUT'")
mata mlib add lspmixw _spmixw_fill_wx_in_data()  , dir("`OUT'")
mata mlib add lspmixw _spmixw_trace_mc()         , dir("`OUT'")
mata mlib add lspmixw _spmixw_effects_sar()      , dir("`OUT'")
mata mlib add lspmixw _spmixw_effects_sdm()      , dir("`OUT'")
mata mlib add lspmixw _spmixw_effects_simple()   , dir("`OUT'")
mata mlib add lspmixw _spmixw_logdet_taylor()    , dir("`OUT'")
mata mlib add lspmixw _spmixw_eval_taylor_lndet(), dir("`OUT'")
mata mlib add lspmixw _spmixw_gamma_proposal_uniform()  , dir("`OUT'")
mata mlib add lspmixw _spmixw_gamma_proposal_adapted()  , dir("`OUT'")
mata mlib add lspmixw _spmixw_eval_cond_sar_conv()      , dir("`OUT'")
mata mlib add lspmixw _spmixw_sar_conv()         , dir("`OUT'")
mata mlib add lspmixw _spmixw_sar_conv_caller()  , dir("`OUT'")
mata mlib add lspmixw _spmixw_simulate_sar_conv(), dir("`OUT'")
mata mlib add lspmixw _spmixw_simulate_sem_conv(), dir("`OUT'")
mata mlib add lspmixw _spmixw_model_probs()      , dir("`OUT'")
mata mlib add lspmixw _spmixw_trace_cross_mc()   , dir("`OUT'")
mata mlib add lspmixw _spmixw_effects_sdm_conv() , dir("`OUT'")
mata mlib add lspmixw _spmixw_effects_sdem_conv(), dir("`OUT'")

mata mlib index

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
