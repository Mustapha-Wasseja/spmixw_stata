*! version 0.1.0  2026-05-08  spmixw: Taylor-series log-det for convex W
*!
*! Implements the order-p Taylor approximation
*!     log|I_N - rho W_c(gamma)|  ~=  - sum_{p>=2} (rho^p / p) tr(W_c^p)
*! where W_c = sum_m gamma_m W_m. Since
*!     tr(W_c^p) = sum_{i_1,...,i_p} gamma_{i_1} ... gamma_{i_p}
*!                 * tr(W_{i_1} ... W_{i_p})
*! we pre-compute the M^p cross-trace terms once (expensive) and then
*! evaluate the log-det rapidly at any (rho, gamma) via Kronecker products.
*!
*! Reference:
*!   Debarsy, N. and LeSage, J. P. (2021). Bayesian model averaging for
*!   spatial autoregressive models based on convex combinations of
*!   different types of connectivity matrices. JBES 40(2), 547-558.
*! Mirrors R port spmixW::log_det_taylor + eval_taylor_lndet.
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_logdet_taylor()
//
// Computes the trace cross-terms tr(W_{i_1} ... W_{i_p}) for orders
// p = 2, ..., max_order. Returns them packed into a single column vector
// of length M^2 + M^3 + ... + M^max_order, in row-major order over
// (i_1, ..., i_p) with i_p innermost.
//
// Inputs:
//   Wstack    N x (N*M)    horizontal stack of M weight matrices,
//                          [W_1, W_2, ..., W_M]
//   M         scalar       number of weight matrices
//   max_order scalar       highest Taylor order (typically 4, up to 8)
//
// Caller is responsible for slicing the packed vector at the right offsets;
// `_spmixw_eval_taylor_lndet()` does that internally.
// ----------------------------------------------------------------------
real colvector _spmixw_logdet_taylor(real matrix W_stack,
                                     real scalar M,
                                     real scalar max_order)
{
    real scalar    N, p, m, a, n_prev, n_new, cnt, total_size, p_offset
    real matrix    accum, new_accum, Aprev, Wm
    real colvector packed

    N = rows(W_stack)
    if (cols(W_stack) != N * M) {
        _error(3200, "_spmixw_logdet_taylor: W_stack must be N x (N*M)")
    }
    if (max_order < 2) {
        _error(3300, "_spmixw_logdet_taylor: max_order must be >= 2")
    }

    // Pre-compute total packed length: M^2 + M^3 + ... + M^max_order
    total_size = 0
    for (p = 2; p <= max_order; p++) total_size = total_size + M^p
    packed = J(total_size, 1, 0)

    // Order 1: accum holds the M base matrices, stacked column-wise.
    accum = W_stack

    p_offset = 0
    for (p = 2; p <= max_order; p++) {
        n_prev = M^(p - 1)
        n_new  = n_prev * M

        // Allocate next-order accum buffer if we'll need it.
        if (p < max_order) {
            new_accum = J(N, N * n_new, 0)
        }

        cnt = 0
        for (a = 1; a <= n_prev; a++) {
            Aprev = accum[., ((a - 1) * N + 1) :: (a * N)]

            for (m = 1; m <= M; m++) {
                Wm = W_stack[., ((m - 1) * N + 1) :: (m * N)]

                // tr(Aprev * Wm) = sum(Aprev :* Wm') -- element-wise dot
                packed[p_offset + cnt + 1] = sum(Aprev :* Wm')
                cnt = cnt + 1

                if (p < max_order) {
                    new_accum[., ((cnt - 1) * N + 1) :: (cnt * N)] = Aprev * Wm
                }
            }
        }

        if (p < max_order) {
            accum = new_accum
        }

        p_offset = p_offset + n_new
    }

    return(packed)
}


// ----------------------------------------------------------------------
// _spmixw_eval_taylor_lndet()
//
// Rapidly evaluates the Taylor log-det at (rho, gamma) using pre-computed
// cross-terms. Builds Kronecker products of gamma incrementally:
//     g_p = kron(g_{p-1}, gamma)   has length M^p
// and computes tr(W_c^p) = g_p' Tp.
//
// Inputs:
//   packed_traces    output of _spmixw_logdet_taylor()
//   M, max_order     same as used to build packed_traces
//   rho              spatial parameter
//   gamma            M x 1 convex weights (caller ensures they sum to 1)
//   use_order        Taylor order to use; 2 <= use_order <= max_order
// ----------------------------------------------------------------------
real scalar _spmixw_eval_taylor_lndet(real colvector packed_traces,
                                      real scalar    M,
                                      real scalar    max_order,
                                      real scalar    rho,
                                      real colvector gamma,
                                      real scalar    use_order)
{
    real scalar    p, n_p, p_offset, lndet, tr_Wcp, rho_p
    real colvector g_prev, g_p, Tp

    if (use_order < 2 | use_order > max_order) {
        _error(3300, "_spmixw_eval_taylor_lndet: use_order outside [2, max_order]")
    }
    if (rows(gamma) != M) {
        _error(3200, "_spmixw_eval_taylor_lndet: gamma length != M")
    }

    p_offset = 0
    g_prev   = gamma
    lndet    = 0
    rho_p    = 1

    for (p = 2; p <= use_order; p++) {
        rho_p = rho_p * rho
        n_p   = M^p

        // g_p = kron(g_{p-1}, gamma), length M^p
        g_p = g_prev # gamma

        Tp     = packed_traces[(p_offset + 1) :: (p_offset + n_p)]
        tr_Wcp = sum(g_p :* Tp)
        lndet  = lndet - rho_p / p * tr_Wcp

        g_prev   = g_p
        p_offset = p_offset + n_p
    }

    return(lndet)
}

end
