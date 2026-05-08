*! version 0.1.0  2026-05-07  spmixw: panel within-transformation
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_demean()
//
// Within-transformation for fixed-effects panel data.
//
// Convention (matches R port spmixW::demean_panel and the LeSage MATLAB
// toolbox): observation (i, t) is stored at position i + (t-1)*N in y and
// in row i + (t-1)*N of X. Equivalently, the data is sorted by time then
// by region.
//
// Inputs (all read):
//   y       NT x 1   dependent variable
//   X       NT x k   covariates (no intercept; intercept is absorbed)
//   N       scalar   number of cross-sectional units
//   TT      scalar   number of time periods
//   mod_    scalar   FE specification: 0=pooled, 1=region, 2=time, 3=both
//
// Outputs (passed by reference, allocated by callee):
//   ywith   NT x 1   demeaned y
//   xwith   NT x k   demeaned X
//   meanny  N x 1    region means of y    (zero if mod_ in {0,2})
//   meannx  N x k    region means of X    (zero if mod_ in {0,2})
//   meanty  T x 1    time means of y      (zero if mod_ in {0,1})
//   meantx  T x k    time means of X      (zero if mod_ in {0,1})
//
// Cross-checked against ../R_package/spmixW/R/demean_panel.R
// ----------------------------------------------------------------------
void _spmixw_demean(real colvector y,
                    real matrix    X,
                    real scalar    N,
                    real scalar    TT,
                    real scalar    mod_,
                    real colvector ywith,
                    real matrix    xwith,
                    real colvector meanny,
                    real matrix    meannx,
                    real colvector meanty,
                    real matrix    meantx)
{
    real scalar    nobs, nvar, j
    real matrix    ymat, xmat_j, repN_meannx, repT_meantx
    real colvector onesT, onesN, grand_mean_x

    nobs = N * TT
    nvar = cols(X)

    if (rows(y) != nobs) _error(3200, "_spmixw_demean: rows(y) != N*T")
    if (rows(X) != nobs) _error(3200, "_spmixw_demean: rows(X) != N*T")
    if (mod_ < 0 | mod_ > 3 | mod_ != trunc(mod_)) {
        _error(3300, "_spmixw_demean: mod_ must be 0, 1, 2, or 3")
    }

    meanny = J(N,  1, 0)
    meannx = J(N,  nvar, 0)
    meanty = J(TT, 1, 0)
    meantx = J(TT, nvar, 0)
    onesT  = J(TT, 1, 1)
    onesN  = J(N,  1, 1)

    // -- Region means ------------------------------------------------------
    // R: matrix(y, nrow=N, ncol=Time) fills column-by-column so result[i,t]
    // = y[i + (t-1)*N]. The Mata equivalent is rowshape(y, TT)' which gives
    // an N x T matrix with the same (i, t) layout.
    if (mod_ == 1 | mod_ == 3) {
        ymat   = rowshape(y, TT)'                  // N x T
        meanny = rowsum(ymat) :/ TT

        for (j = 1; j <= nvar; j++) {
            xmat_j = rowshape(X[., j], TT)'        // N x T
            meannx[., j] = rowsum(xmat_j) :/ TT
        }
    }

    // -- Time means --------------------------------------------------------
    if (mod_ == 2 | mod_ == 3) {
        ymat   = rowshape(y, TT)'                  // N x T
        meanty = (colsum(ymat) :/ N)'

        for (j = 1; j <= nvar; j++) {
            xmat_j = rowshape(X[., j], TT)'        // N x T
            meantx[., j] = (colsum(xmat_j) :/ N)'
        }
    }

    // -- Apply demeaning ---------------------------------------------------
    if (mod_ == 0) {
        ywith = y
        xwith = X
    }
    else if (mod_ == 1) {
        // Subtract region means: kron(ones_T, meanny) -> NT x 1
        ywith = y :- (onesT # meanny)
        xwith = X :- (onesT # meannx)
    }
    else if (mod_ == 2) {
        // Subtract time means: kron(meanty, ones_N) -> NT x 1
        ywith = y :- (meanty # onesN)
        xwith = X :- (meantx # onesN)
    }
    else {
        // Two-way: subtract both, add back the grand mean
        grand_mean_x = mean(X)'                    // k x 1, column vector
        ywith = y :- (onesT # meanny) :- (meanty # onesN) :+ mean(y)
        xwith = X :- (onesT # meannx) :- (meantx # onesN) :+ J(nobs, 1, grand_mean_x')
    }
}

end
