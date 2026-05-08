*! version 0.1.0  2026-05-07  spmixw: Bayesian OLS panel MCMC
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_ols()
//
// Bayesian fixed-effects panel OLS via Gibbs sampling. Mirrors the R port
// spmixW::ols_panel. Heteroscedasticity follows Geweke (1993): each
// observation gets a chi-squared variance weight v_i with prior d.f. rval.
//
// Inputs:
//   y          NT x 1   dependent variable
//   X          NT x k   covariates
//   N          scalar   cross-sectional units
//   TT         scalar   time periods
//   ndraw      scalar   total MCMC draws (incl. burn-in)
//   nomit      scalar   burn-in
//   thin       scalar   thinning interval (>= 1)
//   mod_       scalar   FE: 0=pooled, 1=region, 2=time, 3=both
//   rval       scalar   chi-sq d.f. for heteroscedasticity (0 = homo)
//   nu         scalar   IG shape for sigma^2
//   d0         scalar   IG scale for sigma^2
//   c_beta     k x 1    prior mean of beta
//   C_beta     k x k    prior variance of beta
//
// Outputs (passed by reference):
//   bdraw      n_save x k    retained beta draws
//   sdraw      n_save x 1    retained sigma^2 draws
//   vmean      NT x 1        posterior mean of v_i (ones if homo)
//   beta_pm    k x 1         posterior mean of beta
//   sige_pm    scalar        posterior mean of sigma^2
//
// V convention (Geweke 1993):
//   vi = variance weight (large => more variance for that obs)
//   V  = 1/vi = precision weight
//   GLS pre-multiplication: row i of (y, X) scaled by sqrt(V_i)
// ----------------------------------------------------------------------
void _spmixw_ols(real colvector y,
                 real matrix    X,
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
                 real matrix    bdraw,
                 real colvector sdraw,
                 real colvector vmean,
                 real colvector beta_pm,
                 real scalar    sige_pm)
{
    real scalar    nobs, nvar, n_save, save_idx, iter, homo
    real scalar    sige, nu1, d1, chi
    real colvector ywith, meanny, meanty
    real matrix    xwith, meannx, meantx
    real matrix    Q, AI, prec, XtX
    real colvector Qpc, V, vi, sqrtV, ys, b0, bhat, e, ev, chiv

    nobs = N * TT
    nvar = cols(X)

    if (ndraw <= nomit) _error(3300, "_spmixw_ols: ndraw must exceed nomit")
    if (nomit < 0)      _error(3300, "_spmixw_ols: nomit must be >= 0")
    if (thin  < 1)      _error(3300, "_spmixw_ols: thin must be >= 1")

    // -- Demean for FE -----------------------------------------------------
    _spmixw_demean(y, X, N, TT, mod_,
                   ywith, xwith, meanny, meannx, meanty, meantx)

    // -- Prior precision ---------------------------------------------------
    Q   = invsym(C_beta)
    Qpc = Q * c_beta

    homo = (rval == 0)

    // -- MCMC storage ------------------------------------------------------
    n_save = floor((ndraw - nomit) / thin)
    bdraw  = J(n_save, nvar, 0)
    sdraw  = J(n_save, 1,    0)

    V    = J(nobs, 1, 1)
    vi   = J(nobs, 1, 1)
    sige = 1
    vmean = J(nobs, 1, 0)

    save_idx = 0

    // Pre-compute homoscedastic XtX once (it doesn't change across draws)
    if (homo) {
        XtX = quadcross(xwith, xwith)
    }

    // -- MCMC loop ---------------------------------------------------------
    for (iter = 1; iter <= ndraw; iter++) {

        // ---- Beta draw ---------------------------------------------------
        if (homo) {
            prec = XtX :+ sige :* Q
            AI   = invsym(prec)
            b0   = AI * (quadcross(xwith, ywith) :+ sige :* Qpc)
            bhat = b0 + cholesky(sige :* AI) * rnormal(nvar, 1, 0, 1)
        }
        else {
            sqrtV = sqrt(V)
            // Row-scale by sqrt(V_i): each row i of xwith multiplied by sqrtV[i]
            // Mata broadcasting: (NT x k) :* (NT x 1) does row-wise scaling
            // because element-wise ops auto-broadcast a column on rows.
            // We use explicit Kron-free scaling to be safe under matastrict:
            ys = ywith :* sqrtV
            // Scale X rows: outer product trick — multiply each col by sqrtV
            // Equivalent to diag(sqrtV) * xwith.
            prec = quadcross(xwith, sqrtV :^ 2, xwith) :+ sige :* Q
            AI   = invsym(prec)
            b0   = AI * (quadcross(xwith, sqrtV :^ 2, ywith) :+ sige :* Qpc)
            bhat = b0 + cholesky(sige :* AI) * rnormal(nvar, 1, 0, 1)
        }

        // ---- Sigma^2 draw ------------------------------------------------
        nu1 = nu + nobs
        if (homo) {
            e  = ywith - xwith * bhat
            d1 = d0 + sum(e :^ 2)
        }
        else {
            // Residuals on V-scaled data: e* = sqrtV :* (y - X*bhat)
            e  = ywith - xwith * bhat
            d1 = d0 + sum((sqrtV :* e) :^ 2)
        }
        chi  = rchi2(1, 1, nu1)[1, 1]
        sige = d1 / chi

        // ---- v_i draw (heteroscedastic only) -----------------------------
        if (!homo) {
            ev   = ywith - xwith * bhat
            chiv = rchi2(nobs, 1, rval + 1)
            vi   = ((ev :^ 2) :/ sige :+ rval) :/ chiv
            V    = 1 :/ vi
        }

        // ---- Save draws --------------------------------------------------
        if (iter > nomit & mod(iter - nomit, thin) == 0) {
            save_idx = save_idx + 1
            bdraw[save_idx, .] = bhat'
            sdraw[save_idx]    = sige
            vmean = vmean :+ vi
        }
    }

    vmean = vmean :/ n_save

    // -- Posterior means ---------------------------------------------------
    beta_pm = mean(bdraw)'
    sige_pm = mean(sdraw)
}

end
