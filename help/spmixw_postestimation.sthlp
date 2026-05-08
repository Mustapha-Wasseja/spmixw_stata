{smcl}
{* version 0.1.0  2026-05-08}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "spmixw" "help spmixw"}{...}
{vieweralsosee "spmixw_bma" "help spmixw_bma"}{...}
{viewerjumpto "Description" "spmixw_postestimation##description"}{...}
{viewerjumpto "Syntax for estat" "spmixw_postestimation##syntax"}{...}
{viewerjumpto "Method" "spmixw_postestimation##method"}{...}
{viewerjumpto "Examples" "spmixw_postestimation##examples"}{...}
{viewerjumpto "Stored results" "spmixw_postestimation##results"}{...}
{viewerjumpto "References" "spmixw_postestimation##references"}{...}
{title:Title}

{phang}
{bf:spmixw postestimation} {hline 2} Postestimation tools for {help spmixw}

{marker description}{...}
{title:Postestimation commands}

{pstd}
The following postestimation command is available after {cmd:spmixw}:

{synoptset 18 tabbed}{...}
{synopt:{helpb spmixw_postestimation##estatimpact:estat impact}}LeSage-Pace
direct, indirect, and total spatial effects{p_end}
{synoptline}

{marker syntax}{...}
{title:Syntax for estat}

{marker estatimpact}{...}
{p 4 2 2}
LeSage-Pace direct / indirect / total effects:

{p 8 17 2}
{cmd:estat impact}
[{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr:options}
{synoptline}
{synopt:{cmdab:l:evel(}{it:#}{cmd:)}}credible-interval level; default 95{p_end}
{synopt:{cmdab:iter(}{it:#}{cmd:)}}Hutchinson stochastic-trace samples; default 50{p_end}
{synopt:{cmdab:ord:er(}{it:#}{cmd:)}}truncation order of the spatial multiplier
power series; default 100{p_end}
{synopt:{cmdab:se:ed(}{it:#}{cmd:)}}seed for the stochastic-trace estimator;
if omitted, no {cmd:set seed} is issued{p_end}
{synoptline}

{p 4 4 2}
{cmd:estat impact} is supported after the following models:
{cmd:sar}, {cmd:sdm}, {cmd:sdem}, {cmd:slx}, {cmd:sar_conv}, {cmd:sdm_conv},
{cmd:sdem_conv}.
({cmd:sem} / {cmd:sem_conv} / {cmd:ols}: spatial spillovers either don't
exist or aren't separable into LeSage-Pace decomposition.)

{marker method}{...}
{title:Method}

{pstd}
For SAR-style models y = rho W y + X beta + ..., the marginal effect of
covariate j on y is the matrix (I - rho W)^{-1} (beta_j I + theta_j W).
The LeSage-Pace scalar summary averages diagonal and off-diagonal entries:

{phang2}direct_j   = average diagonal{p_end}
{phang2}indirect_j = average off-diagonal{p_end}
{phang2}total_j    = direct_j + indirect_j

{pstd}
For non-convex SAR / SDM these reduce to closed-form expressions involving
{it:tr(W^r)/N} for r = 0..order, which are estimated stochastically by
Hutchinson's method (Pace and Barry 1999, Barry and Pace 1999).

{pstd}
For SDEM / SLX (no spatial multiplier on X), direct and indirect effects
are exact: direct = beta, indirect = theta, total = beta + theta. No
stochastic trace estimation is needed.

{pstd}
For convex models, the formula above generalises by replacing {it:W}
with W_c(gamma) = sum_m gamma_m W_m. For SAR convex this requires only
{it:tr(W_c^r)/N}; for SDM convex it additionally requires the cross-trace
{it:tr(W_c^r W_m)/N} for each constituent W_m. Both are estimated by
Hutchinson's method, evaluated at the posterior-mean gamma; per-draw rho
and per-draw beta / theta are then plugged into the closed form.

{marker examples}{...}
{title:Examples}

{pstd}
SAR effects (single W):{p_end}
{phang2}{cmd:. xtset region year}{p_end}
{phang2}{cmd:. spmixw y x1 x2, model(sar) w(W) effects(twoway) seed(42)}{p_end}
{phang2}{cmd:. estat impact, seed(202605)}{p_end}

{pstd}
SDM convex effects:{p_end}
{phang2}{cmd:. spmixw y x1 x2, model(sdm_conv) wmats(W1 W2 W3) seed(42)}{p_end}
{phang2}{cmd:. estat impact, seed(202605) iter(80) order(150)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:estat impact} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:r(p)}}number of original covariates{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:r(model)}}model class on which the impacts were computed{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:r(direct)}}(p x 4) posterior mean / SD / lower / upper of direct
effects, one row per original covariate{p_end}
{synopt:{cmd:r(indirect)}}(p x 4) indirect effects{p_end}
{synopt:{cmd:r(total)}}(p x 4) total effects{p_end}
{synopt:{cmd:r(effects_draws)}}(n_save x 3p) raw effects draws, in column
order [direct (p) | indirect (p) | total (p)]; emitted only when
matsize permits{p_end}

{marker references}{...}
{title:References}

{phang}
Barry, R. P. and Pace, R. K. (1999). "Monte Carlo estimates of the log
determinant of large sparse matrices." {it:Linear Algebra and its
Applications} 289(1-3), 41-54.

{phang}
LeSage, J. P. and Pace, R. K. (2009). {it:Introduction to Spatial Econometrics.}
Boca Raton, FL: Chapman & Hall/CRC.

{phang}
Pace, R. K. and Barry, R. P. (1997). "Quick computation of regressions with
a spatially autoregressive dependent variable." {it:Geographical Analysis}
29(3), 232-247.

{title:Author}

{pstd}
Mustapha Wasseja Mohammed{break}
{browse "mailto:muswaseja@gmail.com":muswaseja@gmail.com}

{title:Also see}

{psee}
Online: {help spmixw}, {help spmixw_bma}
