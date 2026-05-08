*! version 0.1.0  2026-05-08  spmixw: exact ln|I_N - rho*W| via eigenvalues
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_logdet_exact()
//
// Pre-computes ln|I_N - rho*W| on a grid of rho values, using the
// eigenvalues of W. For each rho:
//
//     lndet(rho) = sum_i log(1 - rho * Re(lambda_i))
//
// Following the LeSage MATLAB toolbox and R port `log_det_exact`, real parts
// of complex eigenvalues are taken (asymmetric W's eigenvalues come in
// complex-conjugate pairs that combine into real log-dets; using Re()
// matches the reference convention).
//
// Cost is O(N^3) for the eigendecomposition. For large N use a Monte Carlo
// approximation (`_spmixw_logdet_mc`, deferred to a later release).
//
// Inputs:
//   W           N x N    spatial weight matrix (dense; sparse handled in ado)
//   rmin        scalar   lower bound of rho grid (typically -1)
//   rmax        scalar   upper bound of rho grid (typically  1)
//   grid_step   scalar   grid spacing (typically 0.001)
//
// Output:
//   detval      ng x 2   columns are [rho, lndet]
// ----------------------------------------------------------------------
real matrix _spmixw_logdet_exact(real matrix W,
                                 real scalar rmin,
                                 real scalar rmax,
                                 real scalar grid_step)
{
    real scalar       N, ng, j
    real colvector    eigvals_re, rho_grid, lndet
    complex colvector eigvals
    complex matrix    Vc

    N = rows(W)
    if (cols(W) != N) _error(3200, "_spmixw_logdet_exact: W must be square")

    // Eigenvalues of W (asymmetric in general). eigensystem() returns the
    // full system; pass `Vc` for eigenvectors even though we don't use them
    // because the 3-arg form is the documented signature in Stata 13.
    eigensystem(W, Vc, eigvals)
    eigvals_re = Re(eigvals)

    ng = floor((rmax - rmin) / grid_step) + 1
    rho_grid = rmin :+ ((0::(ng-1)) :* grid_step)

    lndet = J(ng, 1, 0)
    for (j = 1; j <= ng; j++) {
        lndet[j] = sum(log(1 :- rho_grid[j] :* eigvals_re))
    }

    return((rho_grid, lndet))
}

end
