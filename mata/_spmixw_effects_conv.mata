*! version 0.1.0  2026-05-08  spmixw: LeSage-Pace effects for convex W models
*!
*! Adds three pieces:
*!   _spmixw_trace_cross_mc(W_c, W_m, maxorder, uiter)
*!     Hutchinson estimator of tr(W_c^r W_m)/N for r = 0..maxorder. r=0 and
*!     r=1 are computed exactly; r>=2 use stochastic trace.
*!   _spmixw_effects_sdm_conv(bdraw, pdraw, tracew, tracew_cm, p, M)
*!     LeSage-Pace decomposition for SDM convex (per-W theta blocks plus a
*!     convex-W spatial multiplier on y).
*!   _spmixw_effects_sdem_conv(bdraw, p, M)
*!     SDEM convex: direct = beta, indirect = sum_m theta_m, total = beta +
*!     indirect. No spatial multiplier on X, so traces are not needed.
*!
*! SAR convex effects use the existing _spmixw_effects_sar with tracew built
*! from the convex-W posterior mean -- no separate function needed.
version 13.0

mata:
mata set matastrict on


// ----------------------------------------------------------------------
// Cross-trace Hutchinson estimator: tr(W_c^r W_m) / N for r = 0..maxorder
// ----------------------------------------------------------------------
real colvector _spmixw_trace_cross_mc(real matrix W_c,
                                      real matrix W_m,
                                      real scalar maxorder,
                                      real scalar uiter)
{
    real scalar    N, k
    real matrix    rv, vu
    real colvector cross

    N = rows(W_c)
    if (cols(W_c) != N | rows(W_m) != N | cols(W_m) != N) {
        _error(3200, "_spmixw_trace_cross_mc: W_c and W_m must be N x N")
    }
    if (maxorder < 2) {
        _error(3300, "_spmixw_trace_cross_mc: maxorder must be >= 2")
    }

    cross = J(maxorder + 1, 1, 0)

    // r = 0: tr(W_m) / N (exact)
    cross[1] = trace(W_m) / N

    // r = 1: tr(W_c W_m) / N = sum_ij W_c[i,j] * W_m[j,i] / N (exact)
    cross[2] = sum(W_c :* W_m') / N

    // r >= 2: Hutchinson. vu_l starts as W_m u_l, then iteratively W_c vu_l.
    rv = rnormal(N, uiter, 0, 1)
    vu = W_m * rv
    vu = W_c * vu                       // now vu_l = W_c W_m u_l (r=1, skipped)

    for (k = 2; k <= maxorder; k++) {
        vu = W_c * vu                   // vu_l = W_c^k W_m u_l
        cross[k + 1] = sum(rv :* vu) / (N * uiter)
    }

    return(cross)
}


// ----------------------------------------------------------------------
// SDM convex effects.
//
// bdraw column layout: first p cols = beta; then M blocks of p cols each
//   (block m = theta_m, the W_m * X coefficient vector). Total cols =
//   p * (M + 1).
//
// For variable j, the marginal effect on y is
//   dy/dx_j = (I - rho W_c)^{-1} (beta_j I + sum_m theta_{m,j} W_m)
//
// Scalar summary:
//   total_j    = (beta_j + sum_m theta_{m,j}) * sum_r rho^r
//   direct_j   = beta_j * sum_r rho^r tr(W_c^r)/N
//              + sum_m theta_{m,j} * sum_r rho^r tr(W_c^r W_m)/N
//   indirect_j = total_j - direct_j
// ----------------------------------------------------------------------
real matrix _spmixw_effects_sdm_conv(real matrix    bdraw,
                                     real colvector pdraw,
                                     real colvector tracew,
                                     real matrix    tracew_cm,
                                     real scalar    p,
                                     real scalar    M)
{
    real scalar    n_save, ntrs, i, m
    real colvector trs, ree, rmat
    real matrix    direct_d, indirect_d, total_d
    real rowvector beta_i, theta_i_m, total_row, direct_row, theta_eff_i

    n_save = rows(bdraw)
    trs    = (1 \ tracew)               // r = 0..maxorder
    ntrs   = rows(trs)
    ree    = (0::(ntrs - 1))

    if (rows(tracew_cm) != ntrs) {
        _error(3300, "_spmixw_effects_sdm_conv: tracew_cm has wrong row count")
    }
    if (cols(tracew_cm) != M) {
        _error(3300, "_spmixw_effects_sdm_conv: tracew_cm has wrong col count")
    }
    if (cols(bdraw) != p * (M + 1)) {
        _error(3300, "_spmixw_effects_sdm_conv: bdraw cols must be p*(M+1)")
    }

    direct_d   = J(n_save, p, 0)
    indirect_d = J(n_save, p, 0)
    total_d    = J(n_save, p, 0)

    for (i = 1; i <= n_save; i++) {
        rmat   = pdraw[i] :^ ree
        beta_i = bdraw[|i, 1 \ i, p|]

        theta_eff_i = J(1, p, 0)
        for (m = 1; m <= M; m++) {
            theta_eff_i = theta_eff_i + ///
                bdraw[|i, p + (m - 1) * p + 1 \ i, p + m * p|]
        }

        total_row  = (beta_i :+ theta_eff_i) :* sum(rmat)
        direct_row = beta_i :* sum(trs :* rmat)
        for (m = 1; m <= M; m++) {
            theta_i_m = bdraw[|i, p + (m - 1) * p + 1 \ i, p + m * p|]
            direct_row = direct_row + ///
                theta_i_m :* sum(tracew_cm[., m] :* rmat)
        }

        total_d[i, .]    = total_row
        direct_d[i, .]   = direct_row
        indirect_d[i, .] = total_row :- direct_row
    }

    return((direct_d, indirect_d, total_d))
}


// ----------------------------------------------------------------------
// SDEM convex effects.
//
// bdraw column layout: same as SDM conv (p beta + M*p theta blocks).
// No spatial multiplier on X, so:
//   direct_j   = beta_j
//   indirect_j = sum_m theta_{m,j}
//   total_j    = beta_j + indirect_j
//
// (Holds for any row-normalised W_m with zero diagonal: average row sum is
// 1, so the W_m*X spillover contributes 1 * theta_{m,j} per W to indirect.)
// ----------------------------------------------------------------------
real matrix _spmixw_effects_sdem_conv(real matrix bdraw,
                                      real scalar p,
                                      real scalar M)
{
    real scalar n_save, m
    real matrix direct_d, indirect_d, total_d

    n_save = rows(bdraw)

    if (cols(bdraw) != p * (M + 1)) {
        _error(3300, "_spmixw_effects_sdem_conv: bdraw cols must be p*(M+1)")
    }

    direct_d   = bdraw[., (1::p)]
    indirect_d = J(n_save, p, 0)

    for (m = 1; m <= M; m++) {
        indirect_d = indirect_d + ///
            bdraw[., ((p + (m - 1) * p + 1)::(p + m * p))]
    }

    total_d = direct_d + indirect_d

    return((direct_d, indirect_d, total_d))
}

end
