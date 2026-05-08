{smcl}
{* version 0.1.0  2026-05-08}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[XT] xtset" "help xtset"}{...}
{vieweralsosee "spmixw_bma" "help spmixw_bma"}{...}
{vieweralsosee "spmixw postestimation" "help spmixw_postestimation"}{...}
{viewerjumpto "Syntax" "spmixw##syntax"}{...}
{viewerjumpto "Description" "spmixw##description"}{...}
{viewerjumpto "Options" "spmixw##options"}{...}
{viewerjumpto "Examples" "spmixw##examples"}{...}
{viewerjumpto "Stored results" "spmixw##results"}{...}
{viewerjumpto "References" "spmixw##references"}{...}
{title:Title}

{phang}
{bf:spmixw} {hline 2} Bayesian spatial panel data models with convex
combinations of weight matrices

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:spmixw}
{depvar} {indepvars}
{ifin}{cmd:,}
{cmdab:m:odel(}{it:name}{cmd:)}
[{it:W-options}]
[{it:options}]

{p 4 4 2}
The data must be {bf:xtset} as a balanced panel before invoking {cmd:spmixw}.

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:m:odel(}{it:name}{cmd:)}} model class — see below{p_end}

{syntab:Weight matrix (one of these is required, depending on the model)}
{synopt:{cmd:w(}{it:matname}{cmd:)}} N x N spatial weight matrix; required for
single-W models ({cmd:sar}, {cmd:sem}, {cmd:sdm}, {cmd:sdem}, {cmd:slx}){p_end}
{synopt:{cmdab:wm:ats(}{it:matnames}{cmd:)}} space-separated list of N x N
matrices; required for {it:_conv} models{p_end}

{syntab:Fixed effects}
{synopt:{cmdab:ef:fects(}{it:type}{cmd:)}} {opt none}, {opt region}, {opt time}, or {opt twoway}; default {cmd:none}{p_end}

{syntab:MCMC}
{synopt:{cmdab:nd:raw(}{it:#}{cmd:)}} total MCMC draws (incl. burn-in); default 5500{p_end}
{synopt:{cmdab:no:mit(}{it:#}{cmd:)}} burn-in draws to discard; default 1500{p_end}
{synopt:{cmdab:t:hin(}{it:#}{cmd:)}} thinning interval; default 1{p_end}
{synopt:{cmdab:r:val(}{it:#}{cmd:)}} chi-sq d.f. for heteroscedasticity; 0 = homoscedastic; default 4{p_end}
{synopt:{cmdab:nu(}{it:#}{cmd:)}} inverse-gamma shape for sigma^2 prior; default 0 (diffuse){p_end}
{synopt:{cmdab:d0(}{it:#}{cmd:)}} inverse-gamma scale for sigma^2 prior; default 0 (diffuse){p_end}
{synopt:{cmdab:bp:rior(}{it:matname}{cmd:)}} k x 1 prior mean for beta; default zero{p_end}
{synopt:{cmdab:bv:ar(}{it:matname}{cmd:)}} k x k prior variance for beta; default I*1e12{p_end}

{syntab:Spatial parameter (rho / lambda)}
{synopt:{cmdab:rmi:n(}{it:#}{cmd:)}} lower bound; default -1 (single-W) / -0.99 (conv){p_end}
{synopt:{cmdab:rma:x(}{it:#}{cmd:)}} upper bound; default  1 (single-W) /  0.99 (conv){p_end}
{synopt:{cmdab:grids:tep(}{it:#}{cmd:)}} grid step for griddy Gibbs sampler (single-W only); default 0.001{p_end}
{synopt:{cmdab:taylor:order(}{it:#}{cmd:)}} truncation order of the Taylor log-det approximation (conv only); default 6{p_end}

{syntab:Other}
{synopt:{cmdab:se:ed(}{it:#}{cmd:)}} non-negative integer seed; if omitted, no {cmd:set seed} is issued{p_end}
{synopt:{cmdab:sa:ving(}{it:filename}{cmd:)}} save raw MCMC draws to a {cmd:.dta}{p_end}
{synopt:{cmdab:l:evel(}{it:#}{cmd:)}} credible-interval level; default 95{p_end}

{synoptline}

{p 4 4 2}
{bf:Models supported:}
{p_end}
{p 8 8 2}
{cmd:ols} {hline 2} Bayesian fixed-effects OLS panel (no spatial component).{break}
{cmd:sar} {hline 2} Spatial autoregressive: y = rho W y + X beta + eps.{break}
{cmd:sem} {hline 2} Spatial error: y = X beta + u, u = lambda W u + eps.{break}
{cmd:sdm} {hline 2} Spatial Durbin: y = rho W y + X beta + W X theta + eps.{break}
{cmd:sdem} {hline 2} Spatial Durbin error: SDM with spatial errors.{break}
{cmd:slx} {hline 2} Spatial lag of X: y = X beta + W X theta + eps.{break}
{cmd:sar_conv}, {cmd:sem_conv}, {cmd:sdm_conv}, {cmd:sdem_conv} {hline 2}
counterparts using a convex combination W_c(gamma) = sum_m gamma_m W_m of M
weight matrices, with gamma drawn jointly with the spatial parameter
(Debarsy and LeSage 2021).
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:spmixw} fits Bayesian fixed-effects panel data models via Markov chain
Monte Carlo, following the methodology in LeSage and Pace (2009) and
Debarsy and LeSage (2021).

{pstd}
For single-W spatial models the spatial parameter rho (SAR / SDM) or lambda
(SEM / SDEM) is sampled by a griddy Gibbs scheme over a fine grid of the
exact log-determinant ln|I - rho W|. For convex-W ({it:_conv}) models the
log-determinant is approximated by a Taylor series in (rho, gamma), and rho
and gamma are jointly updated by Metropolis-Hastings.

{pstd}
For all spatial models {cmd:estat impact} reports LeSage-Pace direct,
indirect, and total effects. For convex models, {cmd:spmixw_bma} averages
results across all 2^M - 1 non-empty subsets of the supplied W matrices.

{pstd}
The coefficient table uses Bayesian column headers ({bf:Post. Mean},
{bf:Post. SD}, {bf:Cred. Interval}). Frequentist columns ({bf:z}, {bf:P>|z|})
are intentionally omitted because they imply a sampling distribution that
the Bayesian posterior does not have. Credible intervals are quantile-based,
read directly from the chain.

{pstd}
The spatial parameter is rendered in the same coefficient table as a
labelled row: {bf:Wy} for SAR / SDM (matching {cmd:xsmle}), {bf:We} for
SEM / SDEM. This makes the output {cmd:estout}- and
{cmd:coefplot}-friendly.

{marker options}{...}
{title:Options}

{phang}
{cmd:effects(}{it:type}{cmd:)} selects the within-transformation:
{cmd:none} = pooled, {cmd:region} = region (cross-sectional) FE,
{cmd:time} = time-period FE, {cmd:twoway} = both. Default {cmd:none}.

{phang}
{cmd:rval(}{it:#}{cmd:)} sets the chi-squared degrees of freedom controlling
heteroscedasticity in the Geweke (1993) student-t treatment. {cmd:rval(0)}
forces homoscedastic errors. Smaller values allow heavier tails / outlier
robustness; default 4. Convex models always run homoscedastic.

{phang}
{cmd:bprior()} and {cmd:bvar()} take Stata matrix names. Use this to impose
informative priors on beta. Diffuse defaults are zero mean and a diagonal
variance of 1e12.

{phang}
{cmd:rmin()}, {cmd:rmax()} bound the support of the spatial parameter. The
griddy Gibbs sampler used for single-W models requires a finite range; the
default {cmd:[-1, 1]} should be tightened (typically to {cmd:[1/lambda_min,
1/lambda_max]} of the W eigenvalues) when known. For row-normalised W with
positive weights the upper bound is exactly 1, so {cmd:rmax(0.99)} is a safe
proper-prior choice.

{phang}
{cmd:taylororder(}{it:#}{cmd:)}, used by convex models only, controls the
truncation order of the Taylor expansion of ln|I - rho W_c(gamma)| around
rho = 0. Higher orders are more accurate but slower; 6-8 is a sensible range
for moderately dense W. Set this lower (e.g. 4) for very sparse W and higher
(e.g. 10) when the spatial parameter is close to the unit boundary.

{marker examples}{...}
{title:Examples}

{pstd}
SAR with two-way fixed effects:{p_end}
{phang2}{cmd:. xtset region year}{p_end}
{phang2}{cmd:. spmixw y x1 x2, model(sar) w(W) effects(twoway) ndraw(5000) nomit(1000) seed(42)}{p_end}
{phang2}{cmd:. estat impact}{p_end}

{pstd}
SDM with informative prior:{p_end}
{phang2}{cmd:. matrix b0 = (0.5 \ 0.5)}{p_end}
{phang2}{cmd:. matrix B0 = I(2) * 0.001}{p_end}
{phang2}{cmd:. spmixw y x1 x2, model(sdm) w(W) effects(region) bprior(b0) bvar(B0)}{p_end}

{pstd}
Convex SAR with three weight matrices:{p_end}
{phang2}{cmd:. spmixw y x1 x2, model(sar_conv) wmats(W_contig W_knn W_dist) ///}{p_end}
{phang2}{cmd:.     effects(twoway) ndraw(15000) nomit(5000) seed(42)}{p_end}
{phang2}{cmd:. estat impact}{p_end}

{pstd}
Heteroscedasticity-robust SAR (rval = 5):{p_end}
{phang2}{cmd:. spmixw y x1 x2, model(sar) w(W) effects(region) rval(5)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:spmixw} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}total observations N*T{p_end}
{synopt:{cmd:e(N_g)}}cross-sectional units N{p_end}
{synopt:{cmd:e(T)}}time periods T{p_end}
{synopt:{cmd:e(k)}}number of covariates{p_end}
{synopt:{cmd:e(M)}}number of W matrices ({it:_conv} models only){p_end}
{synopt:{cmd:e(ndraw)}}total MCMC draws{p_end}
{synopt:{cmd:e(nomit)}}burn-in{p_end}
{synopt:{cmd:e(thin)}}thinning interval{p_end}
{synopt:{cmd:e(rval)}}chi-sq d.f. for heteroscedasticity{p_end}
{synopt:{cmd:e(rho)}}posterior mean of rho (SAR/SDM/SDM_conv) or lambda (SEM/SDEM/SEM_conv/SDEM_conv){p_end}
{synopt:{cmd:e(sige)}}posterior mean of sigma^2{p_end}
{synopt:{cmd:e(rho_acc_rate)}}MH acceptance rate for rho/lambda ({it:_conv} only){p_end}
{synopt:{cmd:e(gam_acc_rate)}}MH acceptance rate for gamma ({it:_conv} only){p_end}
{synopt:{cmd:e(logmarg)}}Chib-style log-marginal likelihood ({it:_conv} only){p_end}
{synopt:{cmd:e(time)}}elapsed seconds in MCMC{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:spmixw}{p_end}
{synopt:{cmd:e(model)}}model class ({cmd:sar}, {cmd:sdm_conv}, ...){p_end}
{synopt:{cmd:e(effects)}}FE specification{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(covars)}}covariates (full list, including W*X augmentations
for SDM / SDEM){p_end}
{synopt:{cmd:e(covars_orig)}}original covariates only ({cmd:sdm} / {cmd:sdem}
/ {cmd:slx} / {cmd:sdm_conv} / {cmd:sdem_conv}){p_end}
{synopt:{cmd:e(wmats)}}list of W matrix names ({it:_conv} models){p_end}
{synopt:{cmd:e(wmatrix)}}W matrix name (single-W models){p_end}
{synopt:{cmd:e(estat_cmd)}}{cmd:spmixw_estat}{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}posterior mean coefficient row vector (covariates plus
spatial parameter row){p_end}
{synopt:{cmd:e(V)}}posterior covariance{p_end}
{synopt:{cmd:e(draws_b)}}retained beta draws (n_save x k){p_end}
{synopt:{cmd:e(draws_s)}}retained sigma^2 draws{p_end}
{synopt:{cmd:e(draws_p)}}retained rho/lambda draws{p_end}
{synopt:{cmd:e(draws_g)}}retained gamma draws ({it:_conv}; n_save x M){p_end}
{synopt:{cmd:e(gam)}}posterior mean gamma row vector ({it:_conv}; 1 x M){p_end}

{p2col 5 22 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{marker references}{...}
{title:References}

{phang}
Debarsy, N. and LeSage, J. P. (2021). "Bayesian model averaging for spatial
autoregressive models based on convex combinations of different types of
connectivity matrices."
{it:Journal of Business & Economic Statistics}, 40(2), 547-558.

{phang}
Geweke, J. (1993). "Bayesian treatment of the independent Student-t linear
model." {it:Journal of Applied Econometrics}, 8(S1), S19-S40.

{phang}
LeSage, J. P. and Pace, R. K. (2009). {it:Introduction to Spatial Econometrics.}
Boca Raton, FL: Chapman & Hall/CRC.

{title:Related implementations}

{pstd}
The same methodology is also available as the R package {bf:spmixW} and
as J. P. LeSage's MATLAB Panel Data Toolbox ({bf:toolbox_panelg}).

{title:Author}

{pstd}
Mustapha Wasseja Mohammed{break}
{browse "mailto:muswaseja@gmail.com":muswaseja@gmail.com}

{title:Also see}

{psee}
Online: {help spmixw_bma}, {help spmixw_postestimation}
