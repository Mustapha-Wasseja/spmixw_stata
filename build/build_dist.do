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

// -- Step 1: ensure lspmixw.mlib is built --------------------------------
di as txt "{hline 70}"
di as txt "Step 1/3: building lspmixw.mlib via build_mlib.do ..."
di as txt "{hline 70}"
do "`PKG'/build/build_mlib.do"

capture confirm file "`PKG'/ado/lspmixw.mlib"
if _rc {
    di as err "build_dist.do: lspmixw.mlib not present after build_mlib.do."
    exit 601
}

// -- Step 2: prepare empty dist/ folder ----------------------------------
di as txt _newline "{hline 70}"
di as txt "Step 2/3: assembling flat dist/ folder at `DIST'/"
di as txt "{hline 70}"

capture mkdir "`PKG'/build"
capture mkdir "`DIST'"

// Wipe any prior contents so dist/ reflects exactly the current package.
local stale : dir "`DIST'" files "*"
foreach f of local stale {
    capture erase "`DIST'/`f'"
}

// -- Step 3: copy every shipped file flat into dist/ ---------------------
local ADO_FILES                  ///
    spmixw.ado                   ///
    spmixw_bma.ado               ///
    spmixw_estat.ado             ///
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

local HELP_FILES                 ///
    spmixw.sthlp                 ///
    spmixw_bma.sthlp             ///
    spmixw_postestimation.sthlp

foreach f of local ADO_FILES {
    quietly copy "`PKG'/ado/`f'" "`DIST'/`f'", replace
    di as txt "  ado:    `f'"
}

foreach f of local HELP_FILES {
    quietly copy "`PKG'/help/`f'" "`DIST'/`f'", replace
    di as txt "  help:   `f'"
}

quietly copy "`PKG'/ado/lspmixw.mlib" "`DIST'/lspmixw.mlib", replace
di as txt "  mlib:   lspmixw.mlib"

// -- Step 4: write flat spmixw.pkg ---------------------------------------
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
file write `pkg' "f lspmixw.mlib" _n
foreach f of local HELP_FILES {
    file write `pkg' "f `f'" _n
}
file close `pkg'
di as txt "  pkg:    spmixw.pkg (flat layout)"

// -- Step 5: write stata.toc ---------------------------------------------
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
