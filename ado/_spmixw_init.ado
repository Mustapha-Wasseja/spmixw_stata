*! version 0.1.0  2026-05-08  spmixw -- runtime Mata initialiser
*!
*! Loads the package's Mata functions into memory by sourcing the shipped
*! _spmixw_*.mata files via Stata's adopath. Used as a fallback for Stata
*! builds where the lspmixw.mlib library cannot be loaded (e.g. Stata 13
*! installs where mata mlib * commands are broken).
*!
*! Idempotent: probes for one of the package's Mata functions and exits
*! immediately if it is already loaded (e.g. via lspmixw.mlib auto-load).

program _spmixw_init
    version 13.0

    // Probe: is one of the package's Mata functions already loaded?
    capture mata: __spmixw_probe = &_spmixw_demean()
    if _rc == 0 {
        capture mata: mata drop __spmixw_probe
        exit 0
    }

    // Source each .mata file in dependency order.
    local files                          ///
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

    foreach f of local files {
        local path : copy local f
        local found : findfile "`f'"
        if "`found'" == "" {
            di as err "spmixw: cannot locate `f' on Stata's adopath."
            di as err "  The package ships .mata source files alongside .ado;"
            di as err "  ensure the install directory is on the adopath."
            exit 601
        }
        quietly do "`found'"
    }

    // Final verify
    capture mata: __spmixw_probe = &_spmixw_demean()
    if _rc {
        di as err "spmixw: Mata functions failed to load after sourcing."
        exit 601
    }
    capture mata: mata drop __spmixw_probe
end
