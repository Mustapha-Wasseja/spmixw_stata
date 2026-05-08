*! version 0.1.0  2026-05-08  spmixw: Mata bridge for SAR/SEM convex
*!
*! Pulls Stata data, runs `_spmixw_sar_conv()` (the convex MCMC kernel,
*! shared between SAR conv and SEM conv since the math is identical),
*! and posts results back via st_matrix / st_numscalar.
*!
*! Lives in its own .mata file so it ends up in lspmixw.mlib and is
*! callable from both _spmixw_sar_conv.ado and _spmixw_sem_conv.ado
*! regardless of which ado is loaded first.
version 13.0

mata:
mata set matastrict on

void _spmixw_sar_conv_caller(string scalar  yvar,
                             string scalar  xvars,
                             real   scalar  M,
                             real   scalar  N,
                             real   scalar  TT,
                             string scalar  wmats_list,
                             real   scalar  ndraw,
                             real   scalar  nomit,
                             real   scalar  thin,
                             real   scalar  mod_,
                             real   scalar  rmin,
                             real   scalar  rmax,
                             real   scalar  taylor_order,
                             string scalar  b_pm_name,
                             string scalar  V_pm_name,
                             string scalar  draws_b_name,
                             string scalar  draws_s_name,
                             string scalar  draws_p_name,
                             string scalar  draws_g_name,
                             string scalar  draws_full_name,
                             string scalar  gam_pm_name,
                             string scalar  sige_pm_scalar,
                             string scalar  rho_pm_scalar,
                             string scalar  rho_acc_scalar,
                             string scalar  gam_acc_scalar,
                             string scalar  log_marg_scalar)
{
    real colvector y, beta_pm, sdraw, pdraw, b_full, gam_pm
    real matrix    X, W_stack, bdraw, gdraw, V_full, draws_full
    real scalar    sige_pm, rho_pm, rho_acc_rate, gam_acc_rate, log_marg
    real scalar    m
    string rowvector W_names

    y = st_data(., yvar)
    X = st_data(., tokens(xvars))

    // Build the stacked W matrix inside Mata. Each individual W is N x N
    // (well within Stata's matsize cap); the stacked N x (N*M) matrix
    // would otherwise exceed matsize for moderate (N, M).
    W_names = tokens(wmats_list)
    if (cols(W_names) != M) {
        _error(3300, "_spmixw_sar_conv_caller: M doesn't match wmats list")
    }
    W_stack = J(N, N * M, 0)
    for (m = 1; m <= M; m++) {
        W_stack[., ((m - 1) * N + 1) :: (m * N)] = st_matrix(W_names[m])
    }

    _spmixw_sar_conv(y, X, W_stack, M, N, TT,
                     ndraw, nomit, thin, mod_,
                     rmin, rmax, taylor_order,
                     bdraw, sdraw, pdraw, gdraw,
                     beta_pm, sige_pm, rho_pm, gam_pm,
                     rho_acc_rate, gam_acc_rate, log_marg)

    draws_full = (bdraw, pdraw)
    b_full     = (beta_pm \ rho_pm)
    V_full     = variance(draws_full)

    st_matrix(b_pm_name,    b_full')
    st_matrix(V_pm_name,    V_full)
    st_matrix(draws_b_name, bdraw)
    st_matrix(draws_s_name, sdraw)
    st_matrix(draws_p_name, pdraw)
    st_matrix(draws_g_name, gdraw)
    st_matrix(draws_full_name, draws_full)
    st_matrix(gam_pm_name,  gam_pm')

    st_numscalar(sige_pm_scalar, sige_pm)
    st_numscalar(rho_pm_scalar,  rho_pm)
    st_numscalar(rho_acc_scalar, rho_acc_rate)
    st_numscalar(gam_acc_scalar, gam_acc_rate)
    st_numscalar(log_marg_scalar, log_marg)
}

end
