*! version 0.1.0  2026-05-08  spmixw: posterior model probabilities (softmax)
*!
*! Converts a vector of log-marginal likelihoods to posterior model
*! probabilities under equal model priors:
*!     p(M_j | y) = exp(L_j - max(L)) / sum_i exp(L_i - max(L))
*! Subtracting the max before exponentiating prevents overflow.
version 13.0

mata:
mata set matastrict on

real colvector _spmixw_model_probs(real colvector lmarginal)
{
    real scalar    K
    real colvector x

    K = rows(lmarginal)
    if (K < 1) _error(3300, "_spmixw_model_probs: empty input")

    x = exp(lmarginal :- max(lmarginal))
    return(x :/ sum(x))
}

end
