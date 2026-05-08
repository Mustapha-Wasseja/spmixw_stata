*! version 0.1.0  2026-05-08  spmixw: Bayesian SEM panel MCMC
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_sem()
//
// Bayesian fixed-effects SEM panel via Gibbs sampling. Mirrors R port
// spmixW::sem_panel.
//
//     y = X*beta + FE + u,   u = lambda * W u + eps,   eps ~ N(0, sigma^2 V)
//
// MCMC kernel (homoscedastic):
//   ys = ywith - lambda*Wy;  xs = xwith - lambda*Wx
//   beta | lambda, sige ~ N( (xs'xs + sige Q)^{-1} (xs'ys + sige Q c),
//                            sige (xs'xs + sige Q)^{-1} )
//   sige | beta, lambda ~ IG  (with d1 from filtered residuals)
//   lambda  via griddy Gibbs (recomputes filtered SSR per coarse-grid point)
//
// Heteroscedastic blocks scale rows of (xwith, ywith) by sqrt(V) and
// recompute Wxs/Wys = W * (sqrtV :* {x, y}) every iteration (V changes).
//
// W is N x N. NT-block W is applied via reshape(N, T) trick — never form
// I_T (x) W explicitly.
// ----------------------------------------------------------------------
void _spmixw_sem(real colvector y,
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
                 real scalar    rmin,
                 real scalar    rmax,
                 real matrix    bdraw,
                 real colvector sdraw,
                 real colvector pdraw,
                 real colvector vmean,
                 real colvector beta_pm,
                 real scalar    sige_pm,
                 real scalar    rho_pm)
{
    real scalar    nobs, nvar, n_save, save_idx, iter, homo, j
    real scalar    sige, rho, rho2, nu1, d1, chi, ysys_rho, ssr_filt, xty_b
    real colvector ywith, ys_filt, ys_h, eps_filt, ys_v, xty_rho
    real matrix    xwith, xs_filt, xs_v, XtX_rho
    real matrix    Y_mat, Xj_mat, WY_mat, WXj_mat
    real colvector Wy, Wy_h
    real matrix    Wx, Wx_h
    real matrix    Q, AI, prec, XtX
    real colvector Qpc, V, vi, sqrtV, b0, bhat, e, ev, chiv
    real colvector dummy_meanny, dummy_meanty
    real colvector meanny, meanty, ev_filt
    real matrix    meannx, meantx
    real matrix    dummy_xwith, dummy_meannx, dummy_meantx, Wy_X_dummy
    // Pre-computed cross-product blocks (Phase 2.D)
    real matrix    A0_, A1s_, A2_
    real colvector c0_, c1_, c2_
    real scalar    d0_, d1_, d2_

    nobs = N * TT
    nvar = cols(X)

    if (ndraw <= nomit) _error(3300, "_spmixw_sem: ndraw must exceed nomit")
    if (rows(W) != N | cols(W) != N) {
        _error(3200, "_spmixw_sem: W must be N x N")
    }

    // -- Demean y, X for FE -----------------------------------------------
    _spmixw_demean(y, X, N, TT, mod_,
                   ywith, xwith, meanny, meannx, meanty, meantx)

    // -- Compute Wy, Wx on the demeaned data (no Wbig) --------------------
    Y_mat  = rowshape(ywith, TT)'                 // N x T
    WY_mat = W * Y_mat
    Wy     = vec(WY_mat)

    Wx = J(nobs, nvar, 0)
    for (j = 1; j <= nvar; j++) {
        Xj_mat  = rowshape(xwith[., j], TT)'      // N x T
        WXj_mat = W * Xj_mat
        Wx[., j] = vec(WXj_mat)
    }

    // -- Prior precision --------------------------------------------------
    Q   = invsym(C_beta)
    Qpc = Q * c_beta
    homo = (rval == 0)

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

    // ------------------------------------------------------------------
    // Build the V-weighted cross-product blocks used by the griddy.
    // For homo (V = 1 forever) the blocks are constants of (X, y, Wx, Wy)
    // and we build them ONCE outside the MCMC loop. For hetero we rebuild
    // them after each V update (inside the hetero branch, before the
    // lambda draw). This eliminates 8+ quadcross calls per iteration in
    // homo, which dominated Stata 13's runtime.
    // ------------------------------------------------------------------
    if (homo) {
        A0_  = quadcross(xwith, xwith)
        A1s_ = quadcross(xwith, Wx) + quadcross(Wx, xwith)
        A2_  = quadcross(Wx, Wx)
        c0_  = quadcross(xwith, ywith)
        c1_  = quadcross(xwith, Wy) + quadcross(Wx, ywith)
        c2_  = quadcross(Wx, Wy)
        d0_  = sum(ywith :^ 2)
        d1_  = sum(ywith :* Wy)
        d2_  = sum(Wy :^ 2)
    }

    // -- MCMC loop --------------------------------------------------------
    for (iter = 1; iter <= ndraw; iter++) {

        if (homo) {
            // ---- Closed-form β,σ² update via pre-built blocks ----------
            // XtX(λ)  = A0 - λ A1s + λ² A2          (filtered design Gram)
            // Xty(λ)  = c0 - λ c1  + λ² c2          (filtered design × y)
            // ys'ys(λ) = d0 - 2λ d1 + λ² d2         (filtered y norm squared)
            // No NT-sized matmuls needed in the inner loop.
            rho2     = rho :* rho
            XtX_rho  = A0_ :- rho :* A1s_ :+ rho2 :* A2_
            xty_rho  = c0_ :- rho :* c1_  :+ rho2 :* c2_
            ysys_rho = d0_ - 2 :* rho :* d1_ + rho2 :* d2_

            // Beta draw
            prec = XtX_rho :+ sige :* Q
            AI   = invsym(prec)
            b0   = AI * (xty_rho :+ sige :* Qpc)
            bhat = b0 + cholesky(sige :* AI) * rnormal(nvar, 1, 0, 1)

            // Sigma² draw via algebraic identity:
            //   sum(eps_filt²) = ys'ys - 2 Xty' β̂ + β̂' XtX β̂
            // (= residual sum-of-squares on filtered data)
            xty_b    = sum(xty_rho :* bhat)
            ssr_filt = ysys_rho - 2 * xty_b + sum(bhat :* (XtX_rho * bhat))
            nu1      = nobs + 2 :* nu
            d1       = 2 :* d0 + ssr_filt
            chi      = rchi2(1, 1, nu1)[1, 1]
            sige     = d1 / chi

            // Lambda draw — same constant blocks
            rho = _spmixw_griddy_rho_sem(detval_panel,
                                          A0_, A1s_, A2_,
                                          c0_, c1_, c2_,
                                          d0_, d1_, d2_,
                                          nobs, nvar, rmin, rmax)
        }
        else {
            // -- Heteroscedastic block -------------------------------------
            sqrtV = sqrt(V)

            // V-scale rows of demeaned (X, y), then recompute W on the
            // scaled data. R port does this every iteration because V moves.
            // xwith :* sqrtV broadcasts the column on rows.
            xs_v = xwith :* sqrtV
            ys_v = ywith :* sqrtV

            // W * (sqrtV :* X)
            Y_mat   = rowshape(ys_v, TT)'
            WY_mat  = W * Y_mat
            Wy_h    = vec(WY_mat)

            Wx_h = J(nobs, nvar, 0)
            for (j = 1; j <= nvar; j++) {
                Xj_mat   = rowshape(xs_v[., j], TT)'
                WXj_mat  = W * Xj_mat
                Wx_h[., j] = vec(WXj_mat)
            }

            xs_filt = xs_v :- rho :* Wx_h
            ys_filt = ys_v :- rho :* Wy_h

            XtX  = quadcross(xs_filt, xs_filt)
            prec = XtX :+ sige :* Q
            AI   = invsym(prec)
            b0   = AI * (quadcross(xs_filt, ys_filt) :+ sige :* Qpc)
            bhat = b0 + cholesky(sige :* AI) * rnormal(nvar, 1, 0, 1)

            // Sigma^2 draw — same W*e shortcut as in the homo branch but
            // using the V-scaled spatial lags (Wy_h, Wx_h) computed at the
            // top of this iteration.
            e        = ys_v - xs_v * bhat
            ev_filt  = Wy_h - Wx_h * bhat
            eps_filt = e :- rho :* ev_filt
            nu1      = nobs + 2 :* nu
            d1       = 2 :* d0 + sum(eps_filt :^ 2)
            chi      = rchi2(1, 1, nu1)[1, 1]
            sige     = d1 / chi

            // V_i draw
            ev   = ywith :* sqrtV - (xwith :* sqrtV) * bhat
            chiv = rchi2(nobs, 1, rval + 1)
            vi   = ((ev :^ 2) :/ sige :+ rval) :/ chiv
            V    = 1 :/ vi

            // Rebuild blocks for the new V — quadcross(A, V, B) computes
            // A' diag(V) B in one pass. The R port's "filter then scale"
            // griddy convention means we use (X, Wx, y, Wy) with diag(V),
            // not the V-scaled versions used by the beta draw above.
            A0_  = quadcross(xwith, V, xwith)
            A1s_ = quadcross(xwith, V, Wx) + quadcross(Wx, V, xwith)
            A2_  = quadcross(Wx, V, Wx)
            c0_  = quadcross(xwith, V, ywith)
            c1_  = quadcross(xwith, V, Wy) + quadcross(Wx, V, ywith)
            c2_  = quadcross(Wx, V, Wy)
            d0_  = sum(V :* (ywith :^ 2))
            d1_  = sum(V :* ywith :* Wy)
            d2_  = sum(V :* (Wy :^ 2))

            // Lambda draw
            rho = _spmixw_griddy_rho_sem(detval_panel,
                                          A0_, A1s_, A2_,
                                          c0_, c1_, c2_,
                                          d0_, d1_, d2_,
                                          nobs, nvar, rmin, rmax)
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

    beta_pm = mean(bdraw)'
    sige_pm = mean(sdraw)
    rho_pm  = mean(pdraw)
}

end
