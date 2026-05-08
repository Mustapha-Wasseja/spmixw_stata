*! version 0.1.0  2026-05-08  spmixw: SEM panel DGP for tests/validation
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_simulate_sem()
//
// Generates a balanced spatial panel from the SEM DGP:
//
//     u = (I_NT - lambda * I_T (x) W)^{-1} eps
//     y = X*beta + FE + u
//
// FE: mu_i = i/N (region), nu_t = t/T (time).
//
// The spatial inverse is applied time-period-by-time-period via lusolve.
//
// Inputs / outputs identical to _spmixw_simulate_sar but the spatial
// parameter acts on the disturbance, not the response.
// ----------------------------------------------------------------------
void _spmixw_simulate_sem(real scalar    N,
                          real scalar    TT,
                          real matrix    W,
                          real scalar    lambda,
                          real colvector beta,
                          real scalar    sigma2,
                          string scalar  effects,
                          real colvector y,
                          real matrix    X)
{
    real scalar    nobs, k
    real colvector fe_vec, eps, u, sfe, tfe
    real matrix    A, Eps_mat, U_mat

    nobs = N * TT
    k    = rows(beta)

    if (rows(W) != N | cols(W) != N) {
        _error(3200, "_spmixw_simulate_sem: W must be N x N")
    }
    if (k < 1) _error(3300, "_spmixw_simulate_sem: beta must have length >= 1")

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

    // u = (I_N - lambda*W)^{-1} eps, applied per period
    A       = I(N) - lambda :* W
    Eps_mat = rowshape(eps, TT)'
    U_mat   = lusolve(A, Eps_mat)
    u       = vec(U_mat)

    y = X * beta + fe_vec + u
}

end
