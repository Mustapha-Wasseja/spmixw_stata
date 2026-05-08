*! spmixw_demo.do  --  end-to-end worked example
*!
*! Generates a synthetic spatial panel with realistic structure (50 regions
*! arranged on a 10x5 grid, 12 years), builds three meaningfully different
*! W matrices, and runs the package's main commands on it.
*!
*! Run after adding the package to your adopath:
*!     adopath ++ "/path/to/spmixw_stata/ado"
*!     do "/path/to/spmixw_stata/examples/spmixw_demo.do"
*!
*! Total runtime: roughly 4-6 minutes.

version 13.0
* Use `clear` (data only) rather than `clear all`, to preserve $SPMIXW_PKG
* and any prior `mata mlib index` so the package's .mlib stays resolvable.
clear
set more off
set seed 20260508

* If running from the dev tree, the package isn't installed -- wire adopath
* to the dev ado/ folder. Prefer the validation-style setup (sources every
* .mata file defensively) when SPMIXW_PKG is set; otherwise rely on the
* installed .mlib via the current adopath.
if "${SPMIXW_PKG}" != "" {
    capture confirm file "${SPMIXW_PKG}/tests/_setup.do"
    if !_rc {
        quietly do "${SPMIXW_PKG}/tests/_setup.do"
    }
    else {
        quietly adopath ++ "${SPMIXW_PKG}/ado"
        capture mata mlib index
    }
}
else {
    capture confirm file "ado/spmixw.ado"
    if !_rc {
        quietly adopath ++ "ado"
        capture mata mlib index
    }
    else {
        capture confirm file "../ado/spmixw.ado"
        if !_rc {
            quietly adopath ++ "../ado"
            capture mata mlib index
        }
    }
}

di as txt "{hline 70}"
di as txt "  spmixw demo -- synthetic spatial panel"
di as txt "{hline 70}"

// ----------------------------------------------------------------------
// 1. Generate region coordinates on a 10x5 grid (50 regions)
// ----------------------------------------------------------------------
local Nregions = 50
local Nrows    = 5
local Ncols    = 10
local TT       = 12

mata:
    Nregions = `Nregions'
    Nrows    = `Nrows'
    Ncols    = `Ncols'
    coords   = J(Nregions, 2, 0)
    for (i = 1; i <= Nregions; i++) {
        coords[i, 1] = ceil(i / Ncols)        // row (1..Nrows)
        coords[i, 2] = i - (coords[i, 1] - 1) * Ncols  // col (1..Ncols)
    }
    st_matrix("coords", coords)
end

// ----------------------------------------------------------------------
// 2. Build three W matrices
//    W_contig : rook contiguity (immediate N/S/E/W neighbours)
//    W_knn    : 4 nearest neighbours by Euclidean distance
//    W_dist   : exp(-d) inverse-distance weights (band of radius 3)
// ----------------------------------------------------------------------
mata:
    Nregions = `Nregions'
    coords   = st_matrix("coords")

    // -- Contiguity (rook)
    W_contig = J(Nregions, Nregions, 0)
    for (i = 1; i <= Nregions; i++) {
        for (j = 1; j <= Nregions; j++) {
            if (i != j) {
                dr = abs(coords[i, 1] - coords[j, 1])
                dc = abs(coords[i, 2] - coords[j, 2])
                if ((dr + dc) == 1) W_contig[i, j] = 1
            }
        }
    }

    // -- KNN (k = 4)
    k = 4
    W_knn = J(Nregions, Nregions, 0)
    for (i = 1; i <= Nregions; i++) {
        d = J(Nregions, 1, 0)
        for (j = 1; j <= Nregions; j++) {
            d[j] = sqrt((coords[i,1]-coords[j,1])^2 + (coords[i,2]-coords[j,2])^2)
        }
        d[i] = .   // exclude self
        ord = order(d, 1)
        for (kk = 1; kk <= k; kk++) {
            W_knn[i, ord[kk]] = 1
        }
    }

    // -- Distance-decay, banded (only neighbours within radius 3)
    W_dist = J(Nregions, Nregions, 0)
    for (i = 1; i <= Nregions; i++) {
        for (j = 1; j <= Nregions; j++) {
            if (i != j) {
                d = sqrt((coords[i,1]-coords[j,1])^2 + (coords[i,2]-coords[j,2])^2)
                if (d <= 3) W_dist[i, j] = exp(-d)
            }
        }
    }

    // Row-normalise each
    for (i = 1; i <= Nregions; i++) {
        rs = sum(W_contig[i, .]); if (rs > 0) W_contig[i, .] = W_contig[i, .] :/ rs
        rs = sum(W_knn[i, .]);    if (rs > 0) W_knn[i, .]    = W_knn[i, .]    :/ rs
        rs = sum(W_dist[i, .]);   if (rs > 0) W_dist[i, .]   = W_dist[i, .]   :/ rs
    }

    st_matrix("W_contig", W_contig)
    st_matrix("W_knn",    W_knn)
    st_matrix("W_dist",   W_dist)
end

// ----------------------------------------------------------------------
// 3. Generate the SAR-conv DGP:
//      y = (I - rho * W_c)^{-1} (X*beta + region_FE + year_FE + eps)
//      W_c = 0.5 * W_contig + 0.5 * W_knn  (W_dist has zero true weight)
//      rho = 0.5,  beta = (1, 0.5)
// ----------------------------------------------------------------------
mata:
    N    = `Nregions'
    TT   = `TT'
    nobs = N * TT

    rho_true   = 0.5
    beta_true  = (1 \ 0.5)
    gamma_true = (0.5 \ 0.5 \ 0)
    sigma2     = 0.5

    rseed(20260508)
    X = rnormal(nobs, 2, 0, 1)

    // Region FE: 0.1 * (region index / N).  Year FE: 0.05 * (year-1)
    region = J(TT, 1, 1) # (1::N)
    year   = (1::TT) # J(N, 1, 1)
    region_fe = 0.10 :* (region :/ N)
    year_fe   = 0.05 :* (year :- 1)

    eps = rnormal(nobs, 1, 0, sqrt(sigma2))
    mu  = X * beta_true + region_fe + year_fe + eps

    Wc = gamma_true[1] :* st_matrix("W_contig") ///
        + gamma_true[2] :* st_matrix("W_knn") ///
        + gamma_true[3] :* st_matrix("W_dist")
    A     = I(N) - rho_true :* Wc
    M_mat = rowshape(mu, TT)'
    Y_mat = lusolve(A, M_mat)
    y     = vec(Y_mat)

    M_dat = (region, year, y, X)

    st_addvar("int",    "region")
    st_addvar("int",    "year")
    st_addvar("double", "y")
    st_addvar("double", "x1")
    st_addvar("double", "x2")
    st_addobs(rows(M_dat))
    st_store(., ., M_dat)
end

xtset region year

di as txt _newline "{hline 70}"
di as txt "  Synthetic data:"
di as txt "    `Nregions' regions on a `Nrows' x `Ncols' grid"
di as txt "    `TT' years, balanced panel"
di as txt "    True DGP: SAR with W_c = 0.5*W_contig + 0.5*W_knn,  rho=0.5"
di as txt "    Three candidate W matrices: W_contig, W_knn, W_dist"
di as txt "{hline 70}"
summarize y x1 x2

// ======================================================================
// Demo 1: Misspecified single-W SAR (use W_dist, which has true weight 0)
// ----------------------------------------------------------------------
// Expected: rho biased (W misspecified). beta still close to truth thanks
// to fixed effects.
// ======================================================================
di as txt _newline "{hline 70}"
di as txt "  Demo 1: SAR with the wrong W (W_dist, true weight = 0)"
di as txt "{hline 70}"

spmixw y x1 x2, model(sar) w(W_dist) effects(twoway) ///
    ndraw(5000) nomit(1000) seed(20260508)

estat impact, seed(20260508)

// ======================================================================
// Demo 2: SAR convex with all three W's
// ----------------------------------------------------------------------
// Expected: gamma close to (0.5, 0.5, 0); rho close to 0.5
// ======================================================================
di as txt _newline "{hline 70}"
di as txt "  Demo 2: SAR convex with three candidate W matrices"
di as txt "{hline 70}"

spmixw y x1 x2, model(sar_conv) wmats(W_contig W_knn W_dist) ///
    effects(twoway) ndraw(8000) nomit(2000) seed(20260508) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

estat impact, seed(20260508)

// ======================================================================
// Demo 3: SDM convex
// ----------------------------------------------------------------------
// Adds W*X spillovers. Useful when X has direct (own-region) and indirect
// (neighbour) effects.
// ======================================================================
di as txt _newline "{hline 70}"
di as txt "  Demo 3: SDM convex (adds W*X spillovers)"
di as txt "{hline 70}"

spmixw y x1 x2, model(sdm_conv) wmats(W_contig W_knn W_dist) ///
    effects(twoway) ndraw(8000) nomit(2000) seed(20260508) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

estat impact, seed(20260508)

// ======================================================================
// Demo 4: BMA over all 7 non-empty W subsets
// ----------------------------------------------------------------------
// Expected: posterior probability concentrates on subsets containing both
// W_contig and W_knn (or close substitutes).
// ======================================================================
di as txt _newline "{hline 70}"
di as txt "  Demo 4: BMA over the 7 non-empty W subsets (sar_conv)"
di as txt "{hline 70}"

spmixw_bma y x1 x2, model(sar_conv) wmats(W_contig W_knn W_dist) ///
    effects(twoway) ndraw(6000) nomit(1500) seed(20260508) ///
    rmin(-0.99) rmax(0.99) taylororder(8)

di as txt _newline "{hline 70}"
di as txt "  Demo complete. Try modifying the DGP, swapping W matrices, or"
di as txt "  switching model(sar_conv) -> model(sdem_conv) to see how"
di as txt "  spillover decompositions change."
di as txt "{hline 70}"
