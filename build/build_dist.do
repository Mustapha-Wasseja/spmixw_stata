*! version 0.1.0  2026-05-08
*! Build a flat distribution folder ready for SSC submission.
*!
*! Produces build/dist/ containing every .ado, .sthlp, .mlib, plus a flat
*! spmixw.pkg (no subdirectory prefixes) and stata.toc. Once built, test with
*!     net install spmixw, from("path-to-build/dist")
*! and zip build/dist/ for emailing to kit.baum@bc.edu.
*!
*! Usage (after setting global SPMIXW_PKG):
*!     do "${SPMIXW_PKG}/build/build_dist.do"

version 13.0

// -- Locate package root (same logic as build_mlib.do) -------------------
local PKG ""

if "${SPMIXW_PKG}" != "" {
    capture confirm file "${SPMIXW_PKG}/ado/spmixw.ado"
    if !_rc local PKG "${SPMIXW_PKG}"
}

if "`PKG'" == "" {
    local me = c(filename)
    if "`me'" != "" {
        local me_unix = subinstr("`me'", "\", "/", .)
        if regexm("`me_unix'", "^(.+)/build/build_dist\.do$") {
            local cand = regexs(1)
            capture confirm file "`cand'/ado/spmixw.ado"
            if !_rc local PKG "`cand'"
        }
    }
}

if "`PKG'" == "" {
    capture confirm file "ado/spmixw.ado"
    if !_rc local PKG "."
}
if "`PKG'" == "" {
    capture confirm file "../ado/spmixw.ado"
    if !_rc local PKG ".."
}

if "`PKG'" == "" {
    di as err "build_dist.do: cannot locate spmixw_stata/ado/."
    di as err "  Set global SPMIXW_PKG or cd to the package root before running."
    exit 601
}

local DIST "`PKG'/build/dist"

// -- Step 1: prepare empty dist/ folder ----------------------------------
di as txt "{hline 70}"
di as txt "Step 1/2: assembling flat dist/ folder at `DIST'/"
di as txt "{hline 70}"

capture mkdir "`PKG'/build"
capture mkdir "`DIST'"

// Wipe any prior contents so dist/ reflects exactly the current package.
local stale : dir "`DIST'" files "*"
foreach f of local stale {
    capture erase "`DIST'/`f'"
}

// -- Step 2: copy every shipped file flat into dist/ ---------------------
//
// The package ships .mata source files alongside the .ado files; the
// _spmixw_init.ado helper sources them at first invocation. This avoids a
// dependency on `mata mlib *` commands, which are broken on some Stata 13
// builds. If a pre-built lspmixw.mlib exists at PKG/ado/, it is bundled
// too so installs that *can* auto-load .mlib functions are slightly faster.
local ADO_FILES                  ///
    spmixw.ado                   ///
    spmixw_bma.ado               ///
    spmixw_estat.ado             ///
    _spmixw_init.ado             ///
    _spmixw_ols.ado              ///
    _spmixw_sar.ado              ///
    _spmixw_sem.ado              ///
    _spmixw_sdm.ado              ///
    _spmixw_sdem.ado             ///
    _spmixw_slx.ado              ///
    _spmixw_sar_conv.ado         ///
    _spmixw_sem_conv.ado         ///
    _spmixw_sdm_conv.ado         ///
    _spmixw_sdem_conv.ado        ///
    _spmixw_estat_impact.ado     ///
    _spmixw_print_table.ado

local MATA_FILES                       ///
    _spmixw_demean.mata                ///
    _spmixw_ols.mata                   ///
    _spmixw_simulate.mata              ///
    _spmixw_logdet_exact.mata          ///
    _spmixw_griddy_rho_sar.mata        ///
    _spmixw_sar.mata                   ///
    _spmixw_simulate_sar.mata          ///
    _spmixw_summary.mata               ///
    _spmixw_linear_interp.mata         ///
    _spmixw_griddy_rho_sem.mata        ///
    _spmixw_sem.mata                   ///
    _spmixw_simulate_sem.mata          ///
    _spmixw_compute_wx.mata            ///
    _spmixw_trace_mc.mata              ///
    _spmixw_effects.mata               ///
    _spmixw_logdet_taylor.mata         ///
    _spmixw_sar_conv.mata              ///
    _spmixw_sar_conv_caller.mata       ///
    _spmixw_simulate_sar_conv.mata     ///
    _spmixw_simulate_sem_conv.mata     ///
    _spmixw_model_probs.mata           ///
    _spmixw_effects_conv.mata

local HELP_FILES                 ///
    spmixw.sthlp                 ///
    spmixw_bma.sthlp             ///
    spmixw_postestimation.sthlp

foreach f of local ADO_FILES {
    quietly copy "`PKG'/ado/`f'" "`DIST'/`f'", replace
    di as txt "  ado:    `f'"
}

foreach f of local MATA_FILES {
    quietly copy "`PKG'/mata/`f'" "`DIST'/`f'", replace
    di as txt "  mata:   `f'"
}

foreach f of local HELP_FILES {
    quietly copy "`PKG'/help/`f'" "`DIST'/`f'", replace
    di as txt "  help:   `f'"
}

// Bundle a pre-built .mlib if one happens to exist; otherwise skip.
local has_mlib 0
capture confirm file "`PKG'/ado/lspmixw.mlib"
if !_rc {
    quietly copy "`PKG'/ado/lspmixw.mlib" "`DIST'/lspmixw.mlib", replace
    di as txt "  mlib:   lspmixw.mlib (bundled)"
    local has_mlib 1
}
else {
    di as txt "  mlib:   (none -- install will source .mata files at runtime)"
}

// -- Step 3: write flat spmixw.pkg ---------------------------------------
tempname pkg
file open `pkg' using "`DIST'/spmixw.pkg", write replace
file write `pkg' "* spmixw.pkg" _n
file write `pkg' "v 3" _n
file write `pkg' "d 'SPMIXW': Bayesian spatial panel data models with convex combinations of weight matrices" _n
file write `pkg' "d" _n
file write `pkg' "d {bf:spmixw} fits Bayesian fixed-effects panel data models for SAR, SEM," _n
file write `pkg' "d SDM, SDEM, SLX and their convex-W counterparts (Debarsy & LeSage 2021)." _n
file write `pkg' "d Companion command {bf:spmixw_bma} runs Bayesian Model Averaging over" _n
file write `pkg' "d non-empty subsets of M weight matrices. {bf:estat impact} reports" _n
file write `pkg' "d LeSage-Pace direct, indirect, and total effects." _n
file write `pkg' "d" _n
file write `pkg' "d Author: Mustapha Wasseja Mohammed" _n
file write `pkg' "d Support: muswaseja@gmail.com" _n
file write `pkg' "d" _n
file write `pkg' "d Distribution-Date: 20260508" _n
file write `pkg' "d" _n
file write `pkg' "d Stata version 13.0 or later" _n
file write `pkg' "" _n
foreach f of local ADO_FILES {
    file write `pkg' "f `f'" _n
}
foreach f of local MATA_FILES {
    file write `pkg' "f `f'" _n
}
if `has_mlib' {
    file write `pkg' "f lspmixw.mlib" _n
}
foreach f of local HELP_FILES {
    file write `pkg' "f `f'" _n
}
file close `pkg'
di as txt "  pkg:    spmixw.pkg (flat layout)"

// -- Step 4: write stata.toc ---------------------------------------------
tempname toc
file open `toc' using "`DIST'/stata.toc", write replace
file write `toc' "* stata.toc for spmixw" _n
file write `toc' "v 3" _n
file write `toc' "d spmixw -- Bayesian spatial panel data models with convex W" _n
file write `toc' "d Mustapha Wasseja Mohammed <muswaseja@gmail.com>" _n
file write `toc' "p spmixw" _n
file close `toc'
di as txt "  toc:    stata.toc"

// -- Done ----------------------------------------------------------------
di as txt _newline "{hline 70}"
di as txt "DIST OK: flat distribution written to:"
di as res "  `DIST'"
di as txt _newline "Test the install on a clean Stata session with:"
di as res `"  net install spmixw, from("`DIST'") replace"'
di as txt "{hline 70}"
