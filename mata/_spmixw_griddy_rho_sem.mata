*! version 0.3.0  2026-05-08  spmixw: griddy Gibbs sampler for SEM lambda
*! Phase 2.D performance: blocks (V-weighted cross-products) are now built
*! by the caller and passed in. _spmixw_sem builds them once outside the
*! MCMC loop for homo (V is constant) and per-iteration for hetero. The
*! griddy itself is just a tight loop over the coarse rho grid.
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_griddy_rho_sem()
//
//   xs(λ)    = X - λ Wx,          ys(λ) = y - λ Wy
//   XtX(λ)   = A0 - λ A1s + λ² A2
//   Xty(λ)   = c0 - λ c1  + λ² c2
//   ys'ys(λ) = d0 - 2λ d1 + λ² d2
//   SSR(λ)   = ys'ys(λ) - Xty(λ)' · b(λ),    b(λ) = XtX(λ)^{-1} Xty(λ)
//
// where the blocks A*, c*, d* are V-weighted cross-products of (X, Wx, y, Wy):
//   A0 = X'V·X     A1s = X'V·Wx + Wx'V·X     A2  = Wx'V·Wx
//   c0 = X'V·y     c1  = X'V·Wy + Wx'V·y     c2  = Wx'V·Wy
//   d0 = y'V·y     d1  = y'V·Wy              d2  = Wy'V·Wy
// ----------------------------------------------------------------------
real scalar _spmixw_griddy_rho_sem(real matrix    detval,
                                   real matrix    A0,
                                   real matrix    A1s,
                                   real matrix    A2,
                                   real colvector c0,
                                   real colvector c1,
                                   real colvector c2,
                                   real scalar    d0,
                                   real scalar    d1,
                                   real scalar    d2,
                                   real scalar    nt,
                                   real scalar    k,
                                   real scalar    rmin,
                                   real scalar    rmax)
{
    real scalar    nmk, nrho, ng, i, isum, target, idraw, rho_i
    real scalar    ysys_g
    real colvector rho_grid, lndet_vals, b_g, Xty_g
    real colvector upper_idx, lower_idx, cdf, mask
    real matrix    XtX_g
    real colvector epet, coarse_grid, epe_fine, z, x_g

    nmk        = (nt - k) / 2
    nrho       = rows(detval)
    rho_grid   = detval[., 1]
    lndet_vals = detval[., 2]

    ng           = floor((rmax - rmin - 0.02) / 0.01) + 1
    coarse_grid  = (rmin + 0.01) :+ ((0::(ng-1)) :* 0.01)
    epet         = J(ng, 1, 0)

    for (i = 1; i <= ng; i++) {
        rho_i   = coarse_grid[i]
        XtX_g   = A0  :- rho_i :* A1s :+ (rho_i^2) :* A2
        Xty_g   = c0  :- rho_i :* c1  :+ (rho_i^2) :* c2
        ysys_g  = d0 - 2 * rho_i * d1 + (rho_i^2) * d2
        b_g     = lusolve(XtX_g, Xty_g)
        epet[i] = ysys_g - sum(Xty_g :* b_g)
    }

    epe_fine = _spmixw_linear_interp(coarse_grid, epet, rho_grid)

    z = lndet_vals :- nmk :* log(epe_fine)
    z = z :- max(z)
    x_g = exp(z)

    // Vectorised trapezoid integration (replaces the per-element loop —
    // Mata's loop interpretation overhead made this expensive).
    upper_idx = 2::nrho
    lower_idx = 1::(nrho - 1)
    isum = sum((rho_grid[upper_idx] + rho_grid[lower_idx]) :*
               (x_g[upper_idx] - x_g[lower_idx])) / 2

    z = abs(x_g :/ isum)

    // Vectorised inverse-CDF sample. cdf is monotone non-decreasing, so
    // {cdf ≤ target} is a contiguous prefix; sum(mask) gives its length,
    // which is the index of the last grid point at or below target.
    target = runiform(1, 1)[1, 1] * sum(z)
    cdf    = runningsum(z)
    mask   = (cdf :<= target)
    idraw  = sum(mask)
    if (idraw < 1)    idraw = 1
    if (idraw > nrho) idraw = nrho

    return(rho_grid[idraw])
}

end
