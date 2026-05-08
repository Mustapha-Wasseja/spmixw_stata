*! version 0.1.0  2026-05-07  spmixw: synthetic panel DGP
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_simulate_ols()
//
// Phase-1 DGP. Generates a non-spatial panel (rho = 0, theta = NULL):
//
//     y_it = X_it * beta + mu_i + nu_t + eps_it
//
// where eps ~ N(0, sigma2), mu_i = i/N, nu_t = t/T.
//
// Output stored at position i + (t-1)*N. Uses the Stata-side seed via
// rseed() set by the caller — no separate seed argument here.
//
// Inputs:
//   N         scalar    cross-sectional units
//   TT        scalar    time periods
//   beta      k x 1     true coefficients
//   sigma2    scalar    error variance
//   effects   string    "none", "region", "time", "twoway"
//
// Outputs (allocated by callee):
//   y         NT x 1    response
//   X         NT x k    covariates
// ----------------------------------------------------------------------
void _spmixw_simulate_ols(real scalar    N,
                          real scalar    TT,
                          real colvector beta,
                          real scalar    sigma2,
                          string scalar  effects,
                          real colvector y,
                          real matrix    X)
{
    real scalar    nobs, k
    real colvector fe_vec, eps, sfe, tfe

    nobs = N * TT
    k    = rows(beta)

    if (k < 1) _error(3300, "_spmixw_simulate_ols: beta must have length >= 1")

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
    y   = X * beta + fe_vec + eps
}

end
