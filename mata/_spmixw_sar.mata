*! version 0.1.0  2026-05-08  spmixw: Bayesian SAR panel MCMC
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_sar()
//
// Bayesian fixed-effects SAR panel via Gibbs sampling. Mirrors R port
// spmixW::sar_panel.
//
//     y_it = rho * (W y)_it + X_it * beta + FE + eps_it
//     eps  ~ N(0, sigma^2 * V),   V = diag(1/v_i)  (Geweke 1993)
//
// W is N x N. The NT x NT operator I_T (x) W is never formed; instead Wy
// is computed by reshaping y to N x T, multiplying by W, and re-vec'ing.
//
// Inputs:
//   y, X, W, N, TT, ndraw, nomit, thin, mod_, rval, nu, d0, c_beta, C_beta:
//     same conventions as _spmixw_ols
//   detval_panel  ng x 2  pre-computed [rho_grid, T*lndet] from
//                          _spmixw_logdet_exact (caller multiplies by T)
//
// Outputs (passed by reference, allocated by callee):
//   bdraw    n_save x k
//   sdraw    n_save x 1
//   pdraw    n_save x 1     retained rho draws
//   vmean    NT x 1         posterior mean of v_i (ones if homo)
//   beta_pm  k x 1          posterior mean of beta
//   sige_pm  scalar         posterior mean of sigma^2
//   rho_pm   scalar         posterior mean of rho
// ----------------------------------------------------------------------
void _spmixw_sar(real colvector y,
                 real matrix    X,
                 real matrix    W,
                 real scalar    N,
                 real scalar    TT,
                 real scalar    ndraw,
                 real scalar    nomit,
                 real scalar    thin,
                 real scalar    mod_,
                 real scalar    rval,
                 real scalar    nu,
                 real scalar    d0,
                 real colvector c_beta,
                 real matrix    C_beta,
                 real matrix    detval_panel,
                 real matrix    bdraw,
                 real colvector sdraw,
                 real colvector pdraw,
                 real colvector vmean,
                 real colvector beta_pm,
                 real scalar    sige_pm,
                 real scalar    rho_pm)
{
    real scalar    nobs, nvar, n_save, save_idx, iter, homo
    real scalar    sige, rho, nu1, d1, chi
    real scalar    epe0, eped, epe0d
    real colvector ywith, wywith, ys
    real matrix    xwith
    real matrix    Y_mat, WY_mat, Q, AI, prec, XtX
    real colvector Wy, Qpc, V, vi, sqrtV, b0, bhat, e, ev, chiv
    real colvector b0g, bdg, e0, ed
    real colvector dummy_meanny, dummy_meanty, dummy_xcol
    real colvector meanny, meanty
    real matrix    meannx, meantx, xs
    real matrix    dummy_xwith, dummy_meannx, dummy_meantx, Wy_X_dummy
    real colvector Wys, ys_h, yss

    nobs = N * TT
    nvar = cols(X)

    if (ndraw <= nomit) _error(3300, "_spmixw_sar: ndraw must exceed nomit")
    if (rows(W) != N | cols(W) != N) {
        _error(3200, "_spmixw_sar: W must be N x N")
    }

    // -- Compute Wy without forming Wbig (Kron) ---------------------------
    // rowshape(y, TT)' produces an N x T matrix where col t = obs at time t
    // (matches R's matrix(y, nrow=N, ncol=Time)). Multiply by W on the left
    // and vec() (column-stacking) to recover the NT layout i + (t-1)*N.
    Y_mat  = rowshape(y, TT)'
    WY_mat = W * Y_mat
    Wy     = vec(WY_mat)

    // -- Demean y, X, Wy --------------------------------------------------
    _spmixw_demean(y,  X,  N, TT, mod_,
                   ywith, xwith, meanny, meannx, meanty, meantx)

    Wy_X_dummy = J(nobs, 1, 0)
    _spmixw_demean(Wy, Wy_X_dummy, N, TT, mod_,
                   wywith, dummy_xwith,
                   dummy_meanny, dummy_meannx,
                   dummy_meanty, dummy_meantx)

    // -- Prior precision --------------------------------------------------
    Q   = invsym(C_beta)
    Qpc = Q * c_beta
    homo = (rval == 0)

    // -- Pre-compute homo XtX (constant across iterations in homo mode) ---
    if (homo) XtX = quadcross(xwith, xwith)

    // -- MCMC storage -----------------------------------------------------
    n_save = floor((ndraw - nomit) / thin)
    bdraw  = J(n_save, nvar, 0)
    sdraw  = J(n_save, 1,    0)
    pdraw  = J(n_save, 1,    0)
    vmean  = J(nobs,   1,    0)

    V    = J(nobs, 1, 1)
    vi   = J(nobs, 1, 1)
    sige = 1
    rho  = 0.5
    save_idx = 0

    // -- MCMC loop --------------------------------------------------------
    for (iter = 1; iter <= ndraw; iter++) {

        if (homo) {
            // ---- Beta draw -------------------------------------------------
            ys = ywith - rho :* wywith
            prec = XtX :+ sige :* Q
            AI   = invsym(prec)
            b0   = AI * (quadcross(xwith, ys) :+ sige :* Qpc)
            bhat = b0 + cholesky(sige :* AI) * rnormal(nvar, 1, 0, 1)

            // ---- Sige draw -------------------------------------------------
            e   = ys - xwith * bhat
            nu1 = nu + nobs
            d1  = d0 + sum(e :^ 2)
            chi  = rchi2(1, 1, nu1)[1, 1]
            sige = d1 / chi

            // ---- Rho draw (griddy Gibbs) -----------------------------------
            // Reuse AI = (X'X + sige*Q)^{-1}
            b0g = AI * (quadcross(xwith, ywith)  :+ sige :* Qpc)
            bdg = AI * (quadcross(xwith, wywith) :+ sige :* Qpc)
            e0  = ywith  - xwith * b0g
            ed  = wywith - xwith * bdg
            epe0  = sum(e0 :^ 2)
            eped  = sum(ed :^ 2)
            epe0d = (ed' * e0)[1, 1]
            rho   = _spmixw_griddy_rho_sar(detval_panel, epe0, eped, epe0d,
                                            nobs, nvar)
        }
        else {
            // ---- Heteroscedastic beta draw --------------------------------
            sqrtV = sqrt(V)
            xs    = xwith :* sqrtV          // row-scale
            ys_h  = ywith :* sqrtV
            Wys   = wywith :* sqrtV
            XtX   = quadcross(xs, xs)
            prec  = XtX :+ sige :* Q
            AI    = invsym(prec)
            yss   = ys_h - rho :* Wys
            b0    = AI * (quadcross(xs, yss) :+ sige :* Qpc)
            bhat  = b0 + cholesky(sige :* AI) * rnormal(nvar, 1, 0, 1)

            // ---- Sige draw ------------------------------------------------
            e   = yss - xs * bhat
            nu1 = nu + nobs
            d1  = d0 + sum(e :^ 2)
            chi  = rchi2(1, 1, nu1)[1, 1]
            sige = d1 / chi

            // ---- v_i draw -------------------------------------------------
            ev   = ywith - rho :* wywith - xwith * bhat
            chiv = rchi2(nobs, 1, rval + 1)
            vi   = ((ev :^ 2) :/ sige :+ rval) :/ chiv
            V    = 1 :/ vi

            // ---- Rho draw -------------------------------------------------
            b0g = AI * (quadcross(xs, ys_h) :+ sige :* Qpc)
            bdg = AI * (quadcross(xs, Wys)  :+ sige :* Qpc)
            e0  = ys_h - xs * b0g
            ed  = Wys  - xs * bdg
            epe0  = sum(e0 :^ 2)
            eped  = sum(ed :^ 2)
            epe0d = (ed' * e0)[1, 1]
            rho   = _spmixw_griddy_rho_sar(detval_panel, epe0, eped, epe0d,
                                            nobs, nvar)
        }

        // ---- Save draws ---------------------------------------------------
        if (iter > nomit & mod(iter - nomit, thin) == 0) {
            save_idx = save_idx + 1
            bdraw[save_idx, .] = bhat'
            sdraw[save_idx]    = sige
            pdraw[save_idx]    = rho
            vmean = vmean :+ vi
        }
    }

    vmean = vmean :/ n_save

    // -- Posterior means --------------------------------------------------
    beta_pm = mean(bdraw)'
    sige_pm = mean(sdraw)
    rho_pm  = mean(pdraw)
}

end
