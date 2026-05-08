*! version 0.1.0  2026-05-08  spmixw: compute WX (per-period spatial lag of X)
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_compute_wx()
//
// Given an N x N spatial weight matrix and an NT x k design matrix sorted
// by (time, region) — observation (i, t) at row i + (t-1)*N — returns the
// NT x k matrix of W applied per time period:
//
//     [Wx]_{i,j,t} = sum_{l=1..N} W_{il} * X_{l,j,t}
//
// Equivalent to (I_T (x) W) * X but never forms the NT x NT operator.
// ----------------------------------------------------------------------
real matrix _spmixw_compute_wx(real matrix X,
                               real matrix W,
                               real scalar N,
                               real scalar TT)
{
    real scalar    nobs, k, j
    real matrix    Wx, Xj_mat, WXj_mat

    nobs = N * TT
    k    = cols(X)

    if (rows(X) != nobs) _error(3200, "_spmixw_compute_wx: rows(X) != N*T")
    if (rows(W) != N | cols(W) != N) {
        _error(3200, "_spmixw_compute_wx: W must be N x N")
    }

    Wx = J(nobs, k, 0)
    for (j = 1; j <= k; j++) {
        Xj_mat   = rowshape(X[., j], TT)'   // N x T
        WXj_mat  = W * Xj_mat
        Wx[., j] = vec(WXj_mat)
    }
    return(Wx)
}

// ----------------------------------------------------------------------
// _spmixw_fill_wx_in_data()
//
// Caller-side convenience: read X columns from `xvars` in the active
// dataset, apply W, and write the results back to `wxvars` in-place.
// `xvars` and `wxvars` are space-separated lists of equal length k.
// ----------------------------------------------------------------------
void _spmixw_fill_wx_in_data(string scalar W_name,
                             string scalar xvars,
                             string scalar wxvars,
                             real   scalar N,
                             real   scalar TT)
{
    real matrix    X, W, Wx

    W  = st_matrix(W_name)
    X  = st_data(., tokens(xvars))
    Wx = _spmixw_compute_wx(X, W, N, TT)
    st_store(., tokens(wxvars), Wx)
}

end
