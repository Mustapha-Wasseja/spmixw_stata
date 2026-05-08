{smcl}
{* version 0.1.0  2026-05-08}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "spmixw" "help spmixw"}{...}
{vieweralsosee "spmixw postestimation" "help spmixw_postestimation"}{...}
{vieweralsosee "[XT] xtset" "help xtset"}{...}
{viewerjumpto "Syntax" "spmixw_bma##syntax"}{...}
{viewerjumpto "Description" "spmixw_bma##description"}{...}
{viewerjumpto "Options" "spmixw_bma##options"}{...}
{viewerjumpto "Examples" "spmixw_bma##examples"}{...}
{viewerjumpto "Stored results" "spmixw_bma##results"}{...}
{viewerjumpto "References" "spmixw_bma##references"}{...}
{title:Title}

{phang}
{bf:spmixw_bma} {hline 2} Bayesian Model Averaging over W subsets for
convex spatial panel models

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:spmixw_bma}
{depvar} {indepvars}
{ifin}{cmd:,}
{cmdab:m:odel(}{it:name}{cmd:)}
{cmdab:wm:ats(}{it:matnames}{cmd:)}
[{it:options}]

{p 4 4 2}
The data must be {bf:xtset} as a balanced panel before invoking
{cmd:spmixw_bma}.

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:m:odel(}{it:name}{cmd:)}} one of {cmd:sar_conv},
{cmd:sem_conv}, {cmd:sdm_conv}, {cmd:sdem_conv}{p_end}
{synopt:{cmdab:wm:ats(}{it:matnames}{cmd:)}} space-separated list of
M >= 2 N x N weight matrices{p_end}

{syntab:Fixed effects}
{synopt:{cmdab:ef:fects(}{it:type}{cmd:)}} {opt none}, {opt region}, {opt time}, or {opt twoway}; default {cmd:none}{p_end}

{syntab:MCMC (per subset fit)}
{synopt:{cmdab:nd:raw(}{it:#}{cmd:)}} total MCMC draws (incl. burn-in); default 15000{p_end}
{synopt:{cmdab:no:mit(}{it:#}{cmd:)}} burn-in draws to discard; default 5000{p_end}
{synopt:{cmdab:t:hin(}{it:#}{cmd:)}} thinning interval; default 1{p_end}
{synopt:{cmdab:rmi:n(}{it:#}{cmd:)}} lower bound on rho/lambda; default -0.99{p_end}
{synopt:{cmdab:rma:x(}{it:#}{cmd:)}} upper bound on rho/lambda; default  0.99{p_end}
{synopt:{cmdab:taylor:order(}{it:#}{cmd:)}} Taylor truncation order; default 6{p_end}

{syntab:Other}
{synopt:{cmdab:se:ed(}{it:#}{cmd:)}} non-negative seed; if omitted, no {cmd:set seed} is issued{p_end}
{synopt:{cmdab:l:evel(}{it:#}{cmd:)}} credible-interval level; default 95{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:spmixw_bma} performs Bayesian Model Averaging over the {it:2^M - 1}
non-empty subsets of {it:M} candidate weight matrices. For each subset,
the matching convex-W spatial panel model is fit by MCMC and a Chib-style
log-marginal likelihood is computed. Posterior model probabilities are
obtained from softmax of those log-marginals, and BMA-weighted point
estimates of the structural parameters (rho/lambda, sigma^2, beta, gamma,
and W*X coefficients for SDM/SDEM convex) are reported.

{pstd}
The procedure follows Debarsy and LeSage (2021): for a given dependent
variable and covariates, you supply a candidate set of M weight matrices
({cmd:wmats(W1 W2 ... WM)}); the procedure enumerates every non-empty
subset, runs the convex-W estimator on each subset's W matrices, and
reports posterior probability over models. With M = 3 there are 7
non-empty subsets; with M = 4 there are 15. The procedure is
embarrassingly parallel across subsets, but the current implementation
runs them sequentially.

{pstd}
The four supported {cmd:model()} values dispatch to the matching
{it:_conv} estimator and accumulate model-specific results. The same
{cmd:wmats()} list is used to enumerate the subsets, with each element
either present (indicator = 1) or absent (indicator = 0) in a given
subset.

{marker options}{...}
{title:Options}

{phang}
{cmd:wmats(}{it:matnames}{cmd:)} expects M >= 2 space-separated names of
existing N x N Stata matrices, where N matches the panel's
cross-sectional dimension as set by {cmd:xtset}.

{phang}
{cmd:effects(}{it:type}{cmd:)}, {cmd:ndraw()}, {cmd:nomit()}, {cmd:thin()},
{cmd:rmin()}, {cmd:rmax()}, and {cmd:taylororder()} are passed verbatim to
the per-subset {cmd:spmixw} fit. See {help spmixw} for their semantics.

{phang}
{cmd:seed(}{it:#}{cmd:)} sets a single seed before the entire run; each
subset fit then uses Mata's RNG state inherited from Stata's. To replicate
a prior run, use the same seed.

{marker examples}{...}
{title:Examples}

{pstd}
SAR convex BMA with three weight matrices ({cmd:7} non-empty subsets):{p_end}
{phang2}{cmd:. xtset region year}{p_end}
{phang2}{cmd:. spmixw_bma y x1 x2, model(sar_conv) wmats(W1 W2 W3) ///}{p_end}
{phang2}{cmd:.     effects(twoway) ndraw(8000) nomit(2000) seed(42)}{p_end}

{pstd}
SDM convex BMA over two W matrices:{p_end}
{phang2}{cmd:. spmixw_bma y x1 x2, model(sdm_conv) wmats(W_contig W_dist) ///}{p_end}
{phang2}{cmd:.     effects(region) ndraw(10000) nomit(2500) seed(42)}{p_end}

{pstd}
After the run, posterior probabilities are in {cmd:e(probs)} and the BMA
point estimates are in {cmd:e(bma_beta)}, {cmd:e(bma_gamma)}, and (for
SDM/SDEM) {cmd:e(bma_theta)}.

{marker results}{...}
{title:Stored results}

{pstd}{cmd:spmixw_bma} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(M)}}number of W matrices supplied{p_end}
{synopt:{cmd:e(nmodels)}}{it:2^M - 1}{p_end}
{synopt:{cmd:e(bma_rho)}}BMA-weighted rho (SAR/SDM_conv) or lambda
(SEM/SDEM_conv){p_end}
{synopt:{cmd:e(bma_sige)}}BMA-weighted sigma^2{p_end}
{synopt:{cmd:e(time)}}elapsed seconds across all subset fits{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:spmixw_bma}{p_end}
{synopt:{cmd:e(model)}}model class{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(covars)}}covariates{p_end}
{synopt:{cmd:e(wmats)}}list of W matrix names{p_end}
{synopt:{cmd:e(effects)}}FE specification{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(probs)}}({it:nmodels} x 1) posterior model probabilities{p_end}
{synopt:{cmd:e(logm)}}({it:nmodels} x 1) Chib log-marginal likelihoods{p_end}
{synopt:{cmd:e(subsets)}}({it:nmodels} x M) 0/1 indicator matrix; row {it:i}
encodes which W matrices are in subset {it:i}{p_end}
{synopt:{cmd:e(bma_beta)}}(1 x k) BMA-weighted beta{p_end}
{synopt:{cmd:e(bma_gamma)}}(1 x M) BMA-weighted gamma. Posterior probability
on subsets that exclude W_m down-weights its component.{p_end}
{synopt:{cmd:e(bma_theta)}}(1 x M*k) BMA-weighted W*X coefficients
(SDM/SDEM convex only){p_end}
{synopt:{cmd:e(beta_all)}}({it:nmodels} x k) per-subset posterior-mean beta{p_end}
{synopt:{cmd:e(rho_all)}}({it:nmodels} x 1) per-subset rho/lambda{p_end}
{synopt:{cmd:e(sige_all)}}({it:nmodels} x 1) per-subset sigma^2{p_end}
{synopt:{cmd:e(gam_all)}}({it:nmodels} x M) per-subset gamma vectors,
zero-padded for absent W matrices{p_end}
{synopt:{cmd:e(theta_all)}}({it:nmodels} x M*k) per-subset W*X coefficients
(SDM/SDEM convex){p_end}

{marker references}{...}
{title:References}

{phang}
Debarsy, N. and LeSage, J. P. (2021). "Bayesian model averaging for spatial
autoregressive models based on convex combinations of different types of
connectivity matrices."
{it:Journal of Business & Economic Statistics}, 40(2), 547-558.

{title:Author}

{pstd}
Mustapha Wasseja Mohammed{break}
{browse "mailto:muswaseja@gmail.com":muswaseja@gmail.com}

{title:Also see}

{psee}
Online: {help spmixw}, {help spmixw_postestimation}
