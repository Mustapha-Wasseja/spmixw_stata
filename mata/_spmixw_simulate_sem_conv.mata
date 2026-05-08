*! version 0.1.0  2026-05-08  spmixw: SEM convex DGP for tests/validation
*!
*! u = (I - lambda * Wc(gamma))^{-1} eps; y = X*beta + FE + u
*! Per-period via lusolve(I_N - lambda * Wc) on each time-period eps slice.
version 13.0

mata:
mata set matastrict on

void _spmixw_simulate_sem_conv(real scalar    N,
                               real scalar    TT,
                               real matrix    W_stack,
                               real scalar    M,
                               real scalar    lambda,
                               real colvector gamma,
                               real colvector beta,
                               real scalar    sigma2,
                               string scalar  effects,
                               real colvector y,
                               real matrix    X)
{
    real scalar    nobs, k, m
    real colvector fe_vec, eps, u, sfe, tfe
    real matrix    Wc, A, Eps_mat, U_mat

    nobs = N * TT
    k    = rows(beta)

    if (rows(W_stack) != N | cols(W_stack) != N * M) {
        _error(3200, "_spmixw_simulate_sem_conv: W_stack must be N x (N*M)")
    }
    if (rows(gamma) != M) {
        _error(3200, "_spmixw_simulate_sem_conv: gamma must have length M")
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

    // Build Wc, then u = (I - lambda*Wc)^{-1} eps applied per period
    Wc = J(N, N, 0)
    for (m = 1; m <= M; m++) {
        Wc = Wc + gamma[m] :* W_stack[., ((m - 1) * N + 1) :: (m * N)]
    }
    A       = I(N) - lambda :* Wc
    Eps_mat = rowshape(eps, TT)'
    U_mat   = lusolve(A, Eps_mat)
    u       = vec(U_mat)

    y = X * beta + fe_vec + u
}

end
