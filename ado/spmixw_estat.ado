*! version 0.1.0  2026-05-08  spmixw — estat dispatcher
*!
*! Wired in via `ereturn local estat_cmd "spmixw_estat"` in each spmixw
*! estimator. Routes `estat <subcommand>` after a spmixw fit.

program spmixw_estat, rclass
    version 13.0

    if "`e(cmd)'" != "spmixw" {
        di as err "spmixw_estat: no spmixw estimation results found"
        exit 301
    }

    gettoken subcmd 0 : 0, parse(" ,")

    if "`subcmd'" == "impact" {
        _spmixw_estat_impact `0'
    }
    else {
        di as err "estat `subcmd' is not implemented for spmixw"
        di as err "  Available subcommands: impact"
        exit 198
    }

    // Promote r() from the subcommand back up
    return add
end
