*! version 0.1.0  2026-05-08  spmixw: Bayesian SAR panel with convex W
*!
*! Estimates a SAR panel where W = sum_m gamma_m W_m, gamma on the simplex.
*! Uses Metropolis-Hastings random walk for rho with adaptive step size,
*! and an MH "coin-flip" reversible-jump proposal for gamma. Log-det
*! approximated via Taylor series of order taylor_order (precomputed once).
*!
*! Reference: Debarsy & LeSage (2021), JBES 40(2), 547-558.
*! Mirrors R port spmixW::sar_conv_panel.
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_gamma_proposal_uniform()
//
// "Coin-flip" proposal used for the first `noo` MH iterations (matches
// LeSage's default warm-up step). For each j in 1..M-1: pick random in
// (0, gamma_j), keep gamma_j, or pick random in (gamma_j, 1). Set
// gnew[M] = 1 - sum(gnew[1..M-1]). If any negative, retry up to 50 times.
// Returns J(0, 1, 0) on failure.
// ----------------------------------------------------------------------
real colvector _spmixw_gamma_proposal_uniform(real colvector gamma, real scalar M)
{
    real scalar    j, attempts, coin
    real colvector gtst, gnew

    gtst = J(M - 1, 1, 0)
    attempts = 0
    while (attempts < 50) {
        for (j = 1; j <= M - 1; j++) {
            coin = runiform(1, 1)[1, 1]
            if (coin <= 1/3)      gtst[j] = runiform(1, 1)[1, 1] * gamma[j]
            else if (coin <= 2/3) gtst[j] = gamma[j]
            else                  gtst[j] = gamma[j] + runiform(1, 1)[1, 1] * (1 - gamma[j])
        }
        gnew = gtst \ (1 - sum(gtst))
        if (min(gnew) >= 0) return(gnew)
        attempts = attempts + 1
    }
    return(J(0, 1, 0))
}


// ----------------------------------------------------------------------
// _spmixw_gamma_proposal_adapted()
//
// After warm-up, propose with a window of width dd*gam_std around current
// gamma (matching the R port's adaptive proposal).
// ----------------------------------------------------------------------
real colvector _spmixw_gamma_proposal_adapted(real colvector gamma,
                                              real scalar    M,
                                              real colvector gam_std,
                                              real scalar    dd)
{
    real scalar    j, attempts, coin, lo, hi
    real colvector gtst, gnew

    gtst = J(M - 1, 1, 0)
    attempts = 0
    while (attempts < 50) {
        for (j = 1; j <= M - 1; j++) {
            coin = runiform(1, 1)[1, 1]
            if (coin <= 1/3) {
                lo = gamma[j] - dd * gam_std[j]
                gtst[j] = lo + runiform(1, 1)[1, 1] * (gamma[j] - lo)
            }
            else if (coin <= 2/3) {
                gtst[j] = gamma[j]
            }
            else {
                hi = gamma[j] + dd * gam_std[j]
                gtst[j] = gamma[j] + runiform(1, 1)[1, 1] * (hi - gamma[j])
            }
        }
        gnew = gtst \ (1 - sum(gtst))
        if (min(gnew) >= 0) return(gnew)
        attempts = attempts + 1
    }
    return(J(0, 1, 0))
}


// ----------------------------------------------------------------------
// _spmixw_eval_cond_sar_conv()
//
// Log-conditional kernel used in MH steps for both rho and gamma:
//   log p(rho, gamma | rest)  =  detm  -  detx  -  ((NT-k)/2) log(epe)
// where
//   omega = (1, rho*gamma_1, ..., rho*gamma_M)'
//   x_filt = xs * omega         (col j obtained by xs[j]_block * omega)
//   y_filt = ys * omega
//   epe    = || y_filt - x_filt * bhat ||^2
//   detx   = 0.5 * log(det(x_filt' x_filt))
//   detm   = T * Taylor-lndet (already T-scaled in `traces_panel`)
// ----------------------------------------------------------------------
real scalar _spmixw_eval_cond_sar_conv(real scalar    rho,
                                       real colvector gamma,
                                       real colvector bhat,
                                       real matrix    ys_filt_basis,
                                       real matrix    xs_filt_basis,
                                       real colvector traces_panel,
                                       real scalar    M,
                                       real scalar    max_order,
                                       real scalar    use_order,
                                       real scalar    nobs,
                                       real scalar    nvar)
{
    real scalar    nmk, detx, detm, epe
    real colvector omega, y_filt, e
    real matrix    x_filt, xpx
    real scalar    j, col_per_var

    nmk         = (nobs - nvar) / 2
    omega       = 1 \ (rho :* gamma)
    col_per_var = M + 1

    y_filt = ys_filt_basis * omega

    x_filt = J(nobs, nvar, 0)
    for (j = 1; j <= nvar; j++) {
        x_filt[., j] = xs_filt_basis[., ((j - 1) * col_per_var + 1) :: (j * col_per_var)] * omega
    }

    e   = y_filt - x_filt * bhat
    epe = sum(e :^ 2)

    xpx  = quadcross(x_filt, x_filt)
    detx = 0.5 * log(det(xpx))
    detm = _spmixw_eval_taylor_lndet(traces_panel, M, max_order,
                                     rho, gamma, use_order)

    return(detm - detx - nmk * log(epe))
}


// ----------------------------------------------------------------------
// _spmixw_sar_conv()
//
// Main MCMC loop. Mirrors sar_conv_panel.R.
//
// W matrices passed as horizontal stack `W_stack` of size N x (N*M).
//
// Outputs (passed by reference):
//   bdraw   n_save x k           retained beta draws
//   sdraw   n_save x 1           retained sigma^2 draws
//   pdraw   n_save x 1           retained rho draws
//   gdraw   n_save x M           retained gamma draws
//   beta_pm k x 1                posterior mean of beta
//   sige_pm scalar               posterior mean of sigma^2
//   rho_pm  scalar               posterior mean of rho
//   gam_pm  M x 1                posterior mean of gamma
//   rho_acc_rate, gam_acc_rate   MH acceptance rates
// ----------------------------------------------------------------------
void _spmixw_sar_conv(real colvector y,
                      real matrix    X,
                      real matrix    W_stack,
                      real scalar    M,
                      real scalar    N,
                      real scalar    TT,
                      real scalar    ndraw,
                      real scalar    nomit,
                      real scalar    thin,
                      real scalar    mod_,
                      real scalar    rmin,
                      real scalar    rmax,
                      real scalar    taylor_order,
                      real matrix    bdraw,
                      real colvector sdraw,
                      real colvector pdraw,
                      real matrix    gdraw,
                      real colvector beta_pm,
                      real scalar    sige_pm,
                      real scalar    rho_pm,
                      real colvector gam_pm,
                      real scalar    rho_acc_rate,
                      real scalar    gam_acc_rate,
                      real scalar    log_marg)
{
    real scalar    nobs, nvar, n_save, save_idx, iter, j, m
    real scalar    sige, rho, cc, dd, rho_acc, gam_acc, alpMH, p_for_save
    real scalar    rho_new, log_u, acc_rate_rho, acc_rate_gam, noo
    real scalar    use_order, col_per_var
    real colvector ywith, gam, gnew, omega, b0, bhat, e, gam_std, Wm_y
    real matrix    xwith, x_filt, AI, prec, xpx, ys_basis, xs_basis
    real matrix    Y_mat, Xj_mat, WY_mat, WXj_mat
    real matrix    meannx, meantx
    real colvector meanny, meanty
    real colvector traces_panel, traces_unscaled
    real colvector dummy_meanny, dummy_meanty
    real matrix    dummy_xwith, dummy_meannx, dummy_meantx, dummy_X
    real colvector y_filt
    real scalar    cond_old, cond_new

    nobs = N * TT
    nvar = cols(X)

    if (ndraw <= nomit) _error(3300, "_spmixw_sar_conv: ndraw <= nomit")
    if (rows(W_stack) != N | cols(W_stack) != N * M) {
        _error(3200, "_spmixw_sar_conv: W_stack must be N x (N*M)")
    }

    // -- Demean for FE -------------------------------------------------
    _spmixw_demean(y, X, N, TT, mod_,
                   ywith, xwith, meanny, meannx, meanty, meantx)

    // -- Pre-compute Taylor traces (T-scaled for panel) ----------------
    traces_unscaled = _spmixw_logdet_taylor(W_stack, M, taylor_order)
    traces_panel    = TT :* traces_unscaled
    use_order       = taylor_order

    // -- Build ys_basis = [ywith, -W_1 ywith, ..., -W_M ywith] ---------
    ys_basis = J(nobs, M + 1, 0)
    ys_basis[., 1] = ywith
    for (m = 1; m <= M; m++) {
        Y_mat  = rowshape(ywith, TT)'
        WY_mat = W_stack[., ((m - 1) * N + 1) :: (m * N)] * Y_mat
        Wm_y   = vec(WY_mat)
        ys_basis[., m + 1] = -Wm_y
    }

    // -- Build xs_basis: NT x ((M+1)*nvar) ----------------------------
    // Layout: cols ((j-1)*(M+1)+1)..(j*(M+1)) hold [x_j, -W_1 x_j, ..., -W_M x_j]
    col_per_var = M + 1
    xs_basis = J(nobs, col_per_var * nvar, 0)
    for (j = 1; j <= nvar; j++) {
        xs_basis[., (j - 1) * col_per_var + 1] = xwith[., j]
        for (m = 1; m <= M; m++) {
            Xj_mat  = rowshape(xwith[., j], TT)'
            WXj_mat = W_stack[., ((m - 1) * N + 1) :: (m * N)] * Xj_mat
            xs_basis[., (j - 1) * col_per_var + m + 1] = -vec(WXj_mat)
        }
    }

    // -- MCMC storage --------------------------------------------------
    n_save = floor((ndraw - nomit) / thin)
    bdraw  = J(n_save, nvar, 0)
    sdraw  = J(n_save, 1,    0)
    pdraw  = J(n_save, 1,    0)
    gdraw  = J(n_save, M,    0)

    rho     = 0.5
    sige    = 1
    gam     = J(M, 1, 1 / M)
    cc      = 0.1
    dd      = 3.0
    rho_acc = 0
    gam_acc = 0
    gam_std = J(M, 1, 0.1)
    noo     = 1000
    save_idx = 0

    // ----------------------- MCMC loop --------------------------------
    for (iter = 1; iter <= ndraw; iter++) {

        // ---- 1. Beta draw -------------------------------------------
        omega  = 1 \ (rho :* gam)
        x_filt = J(nobs, nvar, 0)
        for (j = 1; j <= nvar; j++) {
            x_filt[., j] = xs_basis[., ((j - 1) * col_per_var + 1) :: (j * col_per_var)] * omega
        }
        y_filt = ys_basis * omega

        AI   = invsym(quadcross(x_filt, x_filt))
        b0   = AI * quadcross(x_filt, y_filt)
        bhat = b0 + cholesky(sige :* AI) * rnormal(nvar, 1, 0, 1)

        // ---- 2. Sigma^2 draw ----------------------------------------
        e    = y_filt - x_filt * bhat
        sige = sum(e :^ 2) / rchi2(1, 1, nobs)[1, 1]

        // ---- 3. Rho MH random walk ----------------------------------
        rho_new = rho + cc * rnormal(1, 1, 0, 1)[1, 1]
        if (rho_new > rmin & rho_new < rmax) {
            cond_old = _spmixw_eval_cond_sar_conv(rho, gam, bhat,
                                                  ys_basis, xs_basis,
                                                  traces_panel,
                                                  M, taylor_order, use_order,
                                                  nobs, nvar)
            cond_new = _spmixw_eval_cond_sar_conv(rho_new, gam, bhat,
                                                  ys_basis, xs_basis,
                                                  traces_panel,
                                                  M, taylor_order, use_order,
                                                  nobs, nvar)
            alpMH = cond_new - cond_old
            log_u = log(runiform(1, 1)[1, 1])
            if (log_u < alpMH) {
                rho     = rho_new
                rho_acc = rho_acc + 1
            }
        }

        // Adapt cc
        acc_rate_rho = rho_acc / iter
        if (acc_rate_rho < 0.4) cc = cc / 1.1
        if (acc_rate_rho > 0.6) cc = cc * 1.1
        if (cc < 0.001 | cc > 1.0) cc = 0.1

        // ---- 4. Gamma MH (only if M >= 2) ---------------------------
        if (M >= 2) {
            if (iter <= noo) {
                gnew = _spmixw_gamma_proposal_uniform(gam, M)
            }
            else {
                gnew = _spmixw_gamma_proposal_adapted(gam, M, gam_std, dd)
            }

            if (rows(gnew) == M) {
                cond_old = _spmixw_eval_cond_sar_conv(rho, gam, bhat,
                                                      ys_basis, xs_basis,
                                                      traces_panel,
                                                      M, taylor_order, use_order,
                                                      nobs, nvar)
                cond_new = _spmixw_eval_cond_sar_conv(rho, gnew, bhat,
                                                      ys_basis, xs_basis,
                                                      traces_panel,
                                                      M, taylor_order, use_order,
                                                      nobs, nvar)
                alpMH = cond_new - cond_old
                log_u = log(runiform(1, 1)[1, 1])
                if (log_u < alpMH) {
                    gam     = gnew
                    gam_acc = gam_acc + 1
                }
            }

            // Adapt dd
            acc_rate_gam = gam_acc / iter
            if (acc_rate_gam < 0.10) dd = dd / 1.1
            if (acc_rate_gam > 0.40) dd = dd * 1.1
            if (dd > 3.0) dd = 3.0
            if (dd < 1.0) dd = 1.0

            // Refresh gam_std from recent draws
            if (iter >= noo & save_idx > 10 & mod(iter, noo) == 0) {
                p_for_save = save_idx
                if (p_for_save > 0) {
                    for (j = 1; j <= M; j++) {
                        cond_new = sqrt(variance(gdraw[(1::p_for_save), j]))
                        if (cond_new > 0.01) gam_std[j] = cond_new
                        else                 gam_std[j] = 0.01
                    }
                }
            }
        }

        // ---- Save draws ---------------------------------------------
        if (iter > nomit & mod(iter - nomit, thin) == 0) {
            save_idx = save_idx + 1
            bdraw[save_idx, .] = bhat'
            sdraw[save_idx]    = sige
            pdraw[save_idx]    = rho
            gdraw[save_idx, .] = gam'
        }
    }

    // -- Posterior summaries -------------------------------------------
    beta_pm = mean(bdraw)'
    sige_pm = mean(sdraw)
    rho_pm  = mean(pdraw)
    gam_pm  = mean(gdraw)'
    rho_acc_rate = rho_acc / ndraw
    gam_acc_rate = gam_acc / ndraw

    // -- Log-marginal likelihood (Chib-style; for BMA over W subsets) --
    //
    //   log p(y | M) ~= mean_i [log_cond(rho_i, gamma_i, beta_i)]  +  logC
    //
    // where logC = -log(D) + lgamma(dof) - dof*log(2*pi) - 0.5*log|X'X|_post
    //       D    = 1 - 1/rmin                  (uniform-on-(rmin,1) prior)
    //       dof  = (nobs - nvar) / 2
    //       X'X_post evaluated at omega_post = (1, rho_pm * gam_pm)'
    real scalar    drawpost_sum, dof, lndetx, logC, D
    real colvector omega_post, b_save_i, gam_save_i
    real matrix    x_filt_post, xpx_post
    real scalar    i_save

    drawpost_sum = 0
    for (i_save = 1; i_save <= n_save; i_save++) {
        b_save_i   = bdraw[i_save, .]'
        gam_save_i = gdraw[i_save, .]'
        drawpost_sum = drawpost_sum +
            _spmixw_eval_cond_sar_conv(pdraw[i_save], gam_save_i, b_save_i,
                                        ys_basis, xs_basis, traces_panel,
                                        M, taylor_order, use_order, nobs, nvar)
    }

    omega_post = 1 \ (rho_pm :* gam_pm)
    x_filt_post = J(nobs, nvar, 0)
    for (j = 1; j <= nvar; j++) {
        x_filt_post[., j] = xs_basis[., ((j - 1) * col_per_var + 1) :: (j * col_per_var)] * omega_post
    }
    xpx_post = quadcross(x_filt_post, x_filt_post)
    lndetx   = log(det(xpx_post))
    dof      = (nobs - nvar) / 2
    D        = 1 - 1 / rmin
    logC     = -log(D) + lngamma(dof) - dof * log(2 * pi()) - 0.5 * lndetx
    log_marg = drawpost_sum / n_save + logC
}

end
