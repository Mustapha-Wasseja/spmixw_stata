*! version 0.1.0  2026-05-08  spmixw: posterior summary table from MCMC draws
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_summary()
//
// Given an n_save x k matrix of MCMC draws and a credibility level in
// (0, 1), returns a k x 4 matrix with columns:
//
//     [ posterior mean | posterior std dev | lower quantile | upper quantile ]
//
// Quantiles are the empirical (type-1) quantiles of each column of `draws`
// — i.e. the credible interval is read directly off the chain, not built
// from a normal approximation.
// ----------------------------------------------------------------------
real matrix _spmixw_summary(real matrix draws, real scalar lvl)
{
    real scalar    k, j, n, idx_lo, idx_hi
    real scalar    plo, phi
    real matrix    out
    real colvector sv

    plo = (1 - lvl) / 2
    phi = 1 - plo

    k = cols(draws)
    n = rows(draws)
    out = J(k, 4, 0)

    for (j = 1; j <= k; j++) {
        out[j, 1] = mean(draws[., j])
        out[j, 2] = sqrt(variance(draws[., j]))
        sv = sort(draws[., j], 1)
        idx_lo = ceil(plo * n)
        idx_hi = ceil(phi * n)
        if (idx_lo < 1) idx_lo = 1
        if (idx_hi > n) idx_hi = n
        out[j, 3] = sv[idx_lo]
        out[j, 4] = sv[idx_hi]
    }

    return(out)
}

end
