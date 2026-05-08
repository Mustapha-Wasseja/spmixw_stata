*! version 0.1.0  2026-05-08  spmixw: stochastic trace estimator tr(W^j)/N
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_trace_mc()
//
// Estimates tr(W^j) / N for j = 1..maxorder using the Hutchinson stochastic
// trace estimator (Pace & Barry 1999, Barry & Pace 1999):
//
//     tr(W^j) ≈ (1/q) * sum_{l=1..q} u_l' W^j u_l
//
// where u_l are independent N(0, I_N) vectors. Returned vector is divided
// by N so its entries are tr(W^j)/N — the quantity used in the LeSage-Pace
// scalar-summary effects formulas.
//
// Block-diagonal Wbig = I_T ⊗ W has tr(Wbig^j) = T * tr(W^j), and dividing
// by NT yields tr(W^j)/N — so we work directly with the N x N W and avoid
// instantiating Wbig.
//
// Overrides matching the R port:
//   tracew[1] = 0   (row-normalised W has zero diagonal ⇒ exact tr(W) = 0)
//   tracew[2] = sum(W .* W) / N    (exact tr(W'W) = ||W||_F^2)
//
// Inputs:
//   W         N x N    spatial weight matrix
//   maxorder  scalar   highest power to estimate (R port uses 100)
//   uiter     scalar   number of random vectors (R port uses 50)
//
// Caller is responsible for setting Stata's seed (`set seed N` shares
// state with Mata's RNG).
//
// Returns: maxorder x 1 column vector of tr(W^j)/N estimates.
// ----------------------------------------------------------------------
real colvector _spmixw_trace_mc(real matrix W,
                                real scalar maxorder,
                                real scalar uiter)
{
    real scalar    N, j
    real matrix    rv, wjjju
    real colvector tracew

    N = rows(W)
    if (cols(W) != N) _error(3200, "_spmixw_trace_mc: W must be square")
    if (maxorder < 2) _error(3300, "_spmixw_trace_mc: maxorder must be >= 2")

    rv     = rnormal(N, uiter, 0, 1)
    wjjju  = rv
    tracew = J(maxorder, 1, 0)

    for (j = 1; j <= maxorder; j++) {
        wjjju     = W * wjjju
        // sum over l, i of rv[i,l] * (W^j u_l)[i], divided by (N*uiter):
        // estimates tr(W^j)/N. Using sum() (which returns a real scalar)
        // avoids the mean(colsum(...)) composition that returns a 1x1
        // matrix in some Stata 13 builds and fails to assign to tracew[j]
        // under matastrict.
        tracew[j] = sum(rv :* wjjju) / (N * uiter)
    }

    // Exact overrides
    tracew[1] = 0
    tracew[2] = sum(W :* W) / N

    return(tracew)
}

end
