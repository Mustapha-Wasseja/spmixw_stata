# spmixw examples

Self-contained worked example demonstrating every major capability of the
package. Run it after adding `spmixw_stata/ado` to your Stata adopath:

```stata
adopath ++ "/path/to/spmixw_stata/ado"
do "/path/to/spmixw_stata/examples/spmixw_demo.do"
```

The script (no external data required):

1. **Generates** 50 synthetic regions on a 10×5 grid, 12 years.
2. **Builds three W matrices** from the region coordinates: rook contiguity
   (`W_contig`), 4-nearest-neighbours (`W_knn`), distance-decay (`W_dist`).
3. **Simulates a SAR-conv DGP** where the true convex weight on `W_dist` is
   zero, so `W_dist` is a methodological "control": correctly-specified
   models should down-weight it.
4. **Runs four demos:**
   - **Demo 1**: SAR with a misspecified W (using only `W_dist`). Shows
     biased ρ when the wrong W is chosen.
   - **Demo 2**: SAR convex with all three W's. Recovers γ ≈ (0.5, 0.5, 0).
   - **Demo 3**: SDM convex with W·X spillovers and `estat impact`.
   - **Demo 4**: BMA over all 7 non-empty subsets — posterior model
     probabilities concentrate on subsets containing both `W_contig` and
     `W_knn`.

Total runtime: roughly 4–6 minutes on a modern laptop (Stata 13+).

## What this is *not*

A real-world dataset. The synthetic DGP is deliberately structured so the
true answer is known and the package's inferences can be checked against it.
For applications, replace the synthetic block with `use yourdata.dta`,
construct `W_*` from your spatial structure (geography, trade flows, …), and
proceed with the same `spmixw` / `spmixw_bma` calls.
