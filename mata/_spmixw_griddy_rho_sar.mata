*! version 0.1.0  2026-05-08  spmixw: griddy Gibbs sampler for SAR rho
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_griddy_rho_sar()
//
// Draws rho from its conditional posterior in the SAR panel model using
// the griddy Gibbs scheme of LeSage & Pace (2009, Ch. 5).
//
// Given the pre-computed (rho, lndet) grid (panel-scaled by T) and the
// three concentrated-likelihood quantities {epe0, eped, epe0d}, the
// kernel of the conditional log-posterior at grid point rho_g is
//
//     den(rho_g) = T * lndet(rho_g)
//                 - ((NT - k) / 2) * log(epe0 - 2*rho_g*epe0d + rho_g^2*eped)
//
// where epe0 = e0'e0, eped = ed'ed, epe0d = ed'e0 with
//     e0 = ywith - X * b0g,    b0g = (X'X + sige*Q)^{-1}(X'ywith + sige*Qpc)
//     ed = Wywith - X * bdg,   bdg = (X'X + sige*Q)^{-1}(X'Wywith + sige*Qpc)
//
// Sampling is by inverse-CDF on the discrete grid using the trapezoid rule
// for the normalising integral. The grid step in `detval` controls
// resolution.
//
// Inputs:
//   detval   ng x 2   [rho_grid, lndet_panel] (lndet already multiplied by T)
//   epe0     scalar   e0' e0
//   eped     scalar   ed' ed
//   epe0d    scalar   ed' e0
//   nt       scalar   NT (panel total observations)
//   k        scalar   number of regressors
//
// Returns: a single rho draw (real scalar).
// ----------------------------------------------------------------------
real scalar _spmixw_griddy_rho_sar(real matrix detval,
                                   real scalar epe0,
                                   real scalar eped,
                                   real scalar epe0d,
                                   real scalar nt,
                                   real scalar k)
{
    real scalar    nmk, nrho, isum, target, idraw
    real colvector rho_grid, lndet_vals, z, x
    real colvector upper_idx, lower_idx, cdf, mask

    nmk        = (nt - k) / 2
    nrho       = rows(detval)
    rho_grid   = detval[., 1]
    lndet_vals = detval[., 2]

    // Concentrated log-posterior on the grid
    z = epe0 :- 2 :* rho_grid :* epe0d :+ (rho_grid :^ 2) :* eped
    z = -nmk :* log(z)
    z = lndet_vals :+ z

    // Stabilise: subtract max before exponentiating
    z = z :- max(z)
    x = exp(z)

    // Vectorised trapezoid integration
    upper_idx = 2::nrho
    lower_idx = 1::(nrho - 1)
    isum = sum((rho_grid[upper_idx] + rho_grid[lower_idx]) :*
               (x[upper_idx] - x[lower_idx])) / 2

    z = abs(x :/ isum)

    // Vectorised inverse-CDF sample (cdf is monotone non-decreasing)
    target = runiform(1, 1)[1, 1] * sum(z)
    cdf    = runningsum(z)
    mask   = (cdf :<= target)
    idraw  = sum(mask)
    if (idraw < 1)    idraw = 1
    if (idraw > nrho) idraw = nrho

    return(rho_grid[idraw])
}

end
