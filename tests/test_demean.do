*! Test: _spmixw_demean reproduces R port spmixW::demean_panel for all 4 modes
version 13.0

if "${SPMIXW_PKG}" == "" global SPMIXW_PKG ".."
quietly do "${SPMIXW_PKG}/tests/_setup.do"

di as txt _newline "----- test_demean: _spmixw_demean parity checks -----"

mata:
    // Synthetic 4 x 3 panel (N=4, T=3) with known structure
    N    = 4
    TT   = 3
    nobs = N * TT
    rseed(20260507)
    y = rnormal(nobs, 1, 0, 1)
    X = rnormal(nobs, 2, 0, 1)

    // -- mode 0: pooled (no demeaning) -----------------------------------
    ywith = .; xwith = .
    meanny = .; meannx = .; meanty = .; meantx = .
    _spmixw_demean(y, X, N, TT, 0,
                   ywith, xwith, meanny, meannx, meanty, meantx)
    assert(all(ywith :== y))
    assert(all(rowsum(abs(xwith :- X)) :< 1e-12))

    // -- mode 1: region FE -----------------------------------------------
    _spmixw_demean(y, X, N, TT, 1,
                   ywith, xwith, meanny, meannx, meanty, meantx)
    // Region means of demeaned y must be zero (within numerical tol)
    ymat = rowshape(ywith, TT)'
    region_means = rowsum(ymat) :/ TT
    assert(all(abs(region_means) :< 1e-12))

    // -- mode 2: time FE -------------------------------------------------
    _spmixw_demean(y, X, N, TT, 2,
                   ywith, xwith, meanny, meannx, meanty, meantx)
    ymat = rowshape(ywith, TT)'
    time_means = (colsum(ymat) :/ N)'
    assert(all(abs(time_means) :< 1e-12))

    // -- mode 3: two-way FE ----------------------------------------------
    _spmixw_demean(y, X, N, TT, 3,
                   ywith, xwith, meanny, meannx, meanty, meantx)
    ymat = rowshape(ywith, TT)'
    region_means = rowsum(ymat) :/ TT
    time_means   = (colsum(ymat) :/ N)'
    assert(all(abs(region_means) :< 1e-12))
    assert(all(abs(time_means)   :< 1e-12))

    // -- meanny / meanty consistency with raw y --------------------------
    ywith = .; xwith = .
    _spmixw_demean(y, X, N, TT, 3,
                   ywith, xwith, meanny, meannx, meanty, meantx)
    raw_ymat = rowshape(y, TT)'
    raw_region = rowsum(raw_ymat) :/ TT
    raw_time   = (colsum(raw_ymat) :/ N)'
    assert(all(abs(meanny :- raw_region) :< 1e-12))
    assert(all(abs(meanty :- raw_time)   :< 1e-12))

    "demean tests passed"
end
