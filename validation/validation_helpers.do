*! Shared helpers for validation_*.do scripts (Stata side).
*! Mirrors R_package/spmixW/inst/validation/validation_helpers.R.

version 13.0

// Drop existing definitions so this file is idempotent across re-sourcing.
capture program drop check_param
capture program drop print_header

// ---- check_param --------------------------------------------------------
// Compares an estimate against a known truth, prints a row with PASS/FAIL,
// and returns r(pass) ∈ {0, 1}.
program check_param, rclass
    syntax , Name(string) Estimate(real) Truth(real) Tol(real) ///
        [Reference(string)]

    local err = abs(`estimate' - `truth')
    local pass = (`err' <= `tol')
    local status = cond(`pass', "PASS", "FAIL")
    local refstr = cond("`reference'" == "", "N/A", "`reference'")

    if `pass' {
        di as txt "  " %-20s "`name'" ///
            " | True="  %-8.4f `truth' ///
            " | Stata=" %-8.4f `estimate' ///
            " | R="     %-8s "`refstr'" ///
            " | tol="   %-5.2f `tol' ///
            " | " as txt "PASS"
    }
    else {
        di as txt "  " %-20s "`name'" ///
            " | True="  %-8.4f `truth' ///
            " | Stata=" %-8.4f `estimate' ///
            " | R="     %-8s "`refstr'" ///
            " | tol="   %-5.2f `tol' ///
            " | " as err "FAIL"
    }

    return scalar pass = `pass'
end

// ---- print_header -------------------------------------------------------
// Uses `0' (raw command-line text) instead of `syntax` so commas in the
// title don't get mis-parsed as option separators.
program print_header
    di
    di as txt "{hline 70}"
    di as txt `"  `0'"'
    di as txt "{hline 70}"
    di as txt "  " %-20s "Parameter" " | " %-8s "True" " | " ///
        %-8s "Stata" " | " %-8s "R ref" " | " %-5s "Tol" " | Status"
    di as txt "{hline 70}"
end
