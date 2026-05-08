*! version 0.1.0  2026-05-08  spmixw: SAR convex DGP for tests/validation
*!
*! Generates a SAR panel with W_c(gamma) = sum_m gamma_m W_m, applied
*! per period via lusolve(I_N - rho*W_c) on each time-period vector.
version 13.0

mata:
mata set matastrict on

void _spmixw_simulate_sar_conv(real scalar    N,
                               real scalar    TT,
                               real matrix    W_stack,
                               real scalar    M,
                               real scalar    rho,
                               real colvector gamma,
                               real colvector beta,
                               real scalar    sigma2,
                               string scalar  effects,
                               real colvector y,
                               real matrix    X)
{
    real scalar    nobs, k, m
    real colvector fe_vec, eps, mu, sfe, tfe
    real matrix    Wc, A, M_mat, Y_mat

    nobs = N * TT
    k    = rows(beta)

    if (rows(W_stack) != N | cols(W_stack) != N * M) {
        _error(3200, "_spmixw_simulate_sar_conv: W_stack must be N x (N*M)")
    }
    if (rows(gamma) != M) {
        _error(3200, "_spmixw_simulate_sar_conv: gamma must have length M")
    }

    X = rnormal(nobs, k, 0, 1)

    fe_vec = J(nobs, 1, 0)
    if (effects == "region" | effects == "twoway") {
        sfe    = (1::N) :/ N
        fe_vec = fe_vec :+ (J(TT, 1, 1) # sfe)
    }
    if (effects == "time" | effects == "twoway") {
        tfe    = (1::TT) :/ TT
        fe_vec = fe_vec :+ (tfe # J(N, 1, 1))
    }

    eps = rnormal(nobs, 1, 0, sqrt(sigma2))
    mu  = X * beta + fe_vec + eps

    // Build the convex W_c (N x N), then (I - rho W_c) y_t = mu_t per period
    Wc = J(N, N, 0)
    for (m = 1; m <= M; m++) {
        Wc = Wc + gamma[m] :* W_stack[., ((m - 1) * N + 1) :: (m * N)]
    }
    A     = I(N) - rho :* Wc
    M_mat = rowshape(mu, TT)'
    Y_mat = lusolve(A, M_mat)
    y     = vec(Y_mat)
}

end
