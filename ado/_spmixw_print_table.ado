*! version 0.1.0  2026-05-08  spmixw — shared Bayesian posterior table
*! Internal helper. Prints a posterior-summary table for the given
*! MCMC draws and column names.

program _spmixw_print_table
    version 13.0
    syntax, DRAWS(name) NAMES(string) DEPvar(string) [Level(cilevel)]

    local lvl = `level' / 100

    // Compute mean / SD / credible interval from the draws via Mata.
    tempname summ
    mata: st_matrix("`summ'", ///
        _spmixw_summary(st_matrix("`draws'"), `lvl'))

    local k = rowsof(`summ')
    local nwords : word count `names'
    if `k' != `nwords' {
        di as err "_spmixw_print_table: row count of draws (`k') != words in names() (`nwords')"
        exit 198
    }

    di as txt "{hline 78}"
    di as txt %14s "`depvar'" "  |" ///
        as txt "  Post. Mean   Post. SD   [`level'% Cred. Interval]"
    di as txt "{hline 14}+{hline 63}"

    forvalues j = 1/`k' {
        local nm : word `j' of `names'
        local m  = `summ'[`j', 1]
        local sd = `summ'[`j', 2]
        local lo = `summ'[`j', 3]
        local hi = `summ'[`j', 4]
        di as txt %14s "`nm'" "  |" ///
            as res %12.4f `m' %11.4f `sd' ///
            "    " %9.4f `lo' "  " %9.4f `hi'
    }

    di as txt "{hline 78}"
end
