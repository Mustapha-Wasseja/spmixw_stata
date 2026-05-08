*! version 0.1.0  2026-05-08  spmixw: SAR panel DGP for tests/validation
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_simulate_sar()
//
// Generates a balanced spatial panel from the SAR DGP:
//
//     y = (I_NT - rho * I_T (x) W)^{-1} ( X*beta + FE + eps )
//
// FE: mu_i = i/N (region), nu_t = t/T (time).
//
// To avoid forming the NT x NT operator, the spatial inverse is applied
// time-period-by-time-period via lusolve on (I_N - rho*W).
//
// Inputs:
//   N         scalar    cross-sectional units
//   TT        scalar    time periods
//   W         N x N     spatial weight matrix
//   rho       scalar    SAR parameter
//   beta      k x 1     true coefficients
//   sigma2    scalar    error variance
//   effects   string    "none", "region", "time", "twoway"
//
// Outputs (allocated by callee):
//   y         NT x 1
//   X         NT x k
// ----------------------------------------------------------------------
void _spmixw_simulate_sar(real scalar    N,
                          real scalar    TT,
                          real matrix    W,
                          real scalar    rho,
                          real colvector beta,
                          real scalar    sigma2,
                          string scalar  effects,
                          real colvector y,
                          real matrix    X)
{
    real scalar    nobs, k
    real colvector fe_vec, eps, mu, sfe, tfe
    real matrix    A, M_mat, Y_mat

    nobs = N * TT
    k    = rows(beta)

    if (rows(W) != N | cols(W) != N) {
        _error(3200, "_spmixw_simulate_sar: W must be N x N")
    }
    if (k < 1) _error(3300, "_spmixw_simulate_sar: beta must have length >= 1")

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

    // (I_N - rho*W) y_t = mu_t for each t
    A     = I(N) - rho :* W
    M_mat = rowshape(mu, TT)'
    Y_mat = lusolve(A, M_mat)
    y     = vec(Y_mat)
}

end
