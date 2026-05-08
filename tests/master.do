*! Run the full spmixw test suite. Exits non-zero on any failure.
version 13.0

if "${SPMIXW_PKG}" == "" {
    di as err "Set: global SPMIXW_PKG /path/to/spmixw_stata"
    di as err "before running master.do"
    exit 198
}

// Capture into a local so child do-files cannot clobber it.
local PKG `"${SPMIXW_PKG}"'

di as txt _newline "===== spmixw test suite ====="
di as txt `"  PKG = `PKG'"'

local TESTS test_demean test_ols
local NPASS = 0
local NFAIL = 0

foreach t of local TESTS {
    di as txt _newline "----- `t' -----"
    // Re-establish the global before every child run, in case a prior child
    // (or `clear all` inside one) wiped it.
    global SPMIXW_PKG `"`PKG'"'
    capture noisily do `"`PKG'/tests/`t'.do"'
    if _rc {
        local NFAIL = `NFAIL' + 1
        di as err "FAIL: `t' (rc=`_rc')"
    }
    else {
        local NPASS = `NPASS' + 1
        di as txt "PASS: `t'"
    }
}

di as txt _newline "===== Summary: `NPASS' passed, `NFAIL' failed ====="
if `NFAIL' > 0 exit 9
