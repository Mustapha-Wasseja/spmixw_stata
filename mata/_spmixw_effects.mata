*! version 0.1.0  2026-05-08  spmixw: LeSage-Pace direct/indirect/total effects
version 13.0

mata:
mata set matastrict on

// ======================================================================
// LeSage-Pace scalar-summary spatial effects.
//
// For SAR:  y = ρ W y + X β + ε,  spatial multiplier (I - ρW)^{-1}
//   total_j(rho)    = β_j * Σ_{r=0..R} ρ^r          ≈ β_j / (1 - ρ)
//   direct_j(rho)   = β_j * Σ_{r=0..R} (tr(W^r)/N) * ρ^r
//   indirect_j(rho) = total_j - direct_j
//
// For SDM:  y = ρ W y + X β + W X θ + ε
//   total_j    = (β_j + θ_j) * Σ ρ^r
//   direct_j   = β_j * Σ (tr(W^r)/N) * ρ^r + θ_j * Σ (tr(W^{r+1})/N) * ρ^r
//   indirect_j = total_j - direct_j
//
// For SDEM, SLX (no spatial multiplier on X):
//   direct_j   = β_j     (exact, no MC error)
//   indirect_j = θ_j     (exact)
//   total_j    = β_j + θ_j
//
// Convention for `bdraw`:
//   SAR  : columns 1..p are the β coefficients (no θ).
//   SDM  : columns 1..p are β; columns (p+1)..(2p) are θ.
//   SDEM : same column layout as SDM.
//   SLX  : same column layout as SDM.
//
// Each effects function returns an n_save x (3p) matrix laid out as
//   [ direct(p cols) | indirect(p cols) | total(p cols) ].
// ======================================================================


// ----------------------------------------------------------------------
// SAR effects
// ----------------------------------------------------------------------
real matrix _spmixw_effects_sar(real matrix    bdraw,
                                real colvector pdraw,
                                real colvector tracew,
                                real scalar    p)
{
    real scalar    n_save, ntrs, i
    real colvector trs, ree, rmat
    real matrix    direct_d, indirect_d, total_d
    real rowvector beta_i, total_row, direct_row

    n_save = rows(bdraw)
    trs    = (1 \ tracew)               // length maxorder + 1
    ntrs   = rows(trs)
    ree    = (0::(ntrs - 1))            // 0, 1, ..., maxorder

    direct_d   = J(n_save, p, 0)
    indirect_d = J(n_save, p, 0)
    total_d    = J(n_save, p, 0)

    for (i = 1; i <= n_save; i++) {
        rmat       = pdraw[i] :^ ree
        beta_i     = bdraw[|i, 1 \ i, p|]
        total_row  = beta_i :* sum(rmat)
        direct_row = beta_i :* sum(trs :* rmat)

        total_d[i, .]    = total_row
        direct_d[i, .]   = direct_row
        indirect_d[i, .] = total_row :- direct_row
    }

    return((direct_d, indirect_d, total_d))
}


// ----------------------------------------------------------------------
// SDM effects
// ----------------------------------------------------------------------
real matrix _spmixw_effects_sdm(real matrix    bdraw,
                                real colvector pdraw,
                                real colvector tracew,
                                real scalar    p)
{
    real scalar    n_save, ntrs, i
    real colvector trs, trs_shift, ree, rmat
    real matrix    direct_d, indirect_d, total_d
    real rowvector beta_i, theta_i, total_row, direct_row

    n_save    = rows(bdraw)
    trs       = (1 \ tracew)            // [1, 0, tr(W^2)/N, tr(W^3)/N, ...]
    trs_shift = (tracew \ 0)            // [0, tr(W^2)/N, tr(W^3)/N, ..., 0]
    ntrs      = rows(trs)
    ree       = (0::(ntrs - 1))

    direct_d   = J(n_save, p, 0)
    indirect_d = J(n_save, p, 0)
    total_d    = J(n_save, p, 0)

    for (i = 1; i <= n_save; i++) {
        rmat       = pdraw[i] :^ ree
        beta_i     = bdraw[|i, 1     \ i, p     |]
        theta_i    = bdraw[|i, p + 1 \ i, 2 * p |]

        total_row  = (beta_i :+ theta_i) :* sum(rmat)
        direct_row = beta_i  :* sum(trs       :* rmat) ///
                   + theta_i :* sum(trs_shift :* rmat)

        total_d[i, .]    = total_row
        direct_d[i, .]   = direct_row
        indirect_d[i, .] = total_row :- direct_row
    }

    return((direct_d, indirect_d, total_d))
}


// ----------------------------------------------------------------------
// SDEM / SLX effects — no spatial multiplier on X.
//   direct  = β      (= bdraw[, 1..p])
//   indirect= θ      (= bdraw[, p+1..2p])
//   total   = β + θ
// ----------------------------------------------------------------------
real matrix _spmixw_effects_simple(real matrix bdraw, real scalar p)
{
    real matrix direct_d, indirect_d, total_d

    direct_d   = bdraw[., (1::p)]
    indirect_d = bdraw[., ((p+1)::(2*p))]
    total_d    = direct_d :+ indirect_d

    return((direct_d, indirect_d, total_d))
}

end
