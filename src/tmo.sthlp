{smcl}
{* *! version 1.0.0 2026-10-29}{...}
{title:Title}

{pstd}
{hi:tmo} {hline 2} Estimating standard errors via the Thresholding Multiple Outcomes (TMO) method.

{title:Syntax}

{p 8 16 2}
{cmd:tmo}, {cmd:cmd}({it:cmdline}) {cmd:x}({it:{help varname}}) {cmd:ylist}({it:{help varlist}}) {cmd:{opt i:dvar}}({it:{help varname}}) [{help tmo##options_table:options}]
{p_end}

{marker options_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt cmd(cmdline)}} {it:cmdline} is the command that produces the
regression of interest{p_end}
{synopt : } - {cmd:tmo} currently supports regressions using {help regress},
{help reghdfe}, {help ivreg2}, or {help ivreghdfe}{p_end}

{synopt :{opth x(varname)}} regressor of interest in {cmd:cmd()} for which to
estimate TMO standard errors{p_end}
{synopt : } - with {cmd:ivreg2} or {cmd:ivreghdfe}, it must be the only
endogenous regressor{p_end}

{synopt :{opth ylist(varlist)}} list of auxiliary outcomes to use in {cmd:tmo}
{p_end}
{synopt : } - must not include the outcome of {cmd:cmd()} or any variable used
in its regressors or absorbed fixed effects{p_end}

{synopt :{opth i:dvar(varname)}} numeric or string location identifier; must be unique
(within each {it:t} for panel case){p_end}

{synopt :{opt miss:limit(#)}} maximum share of locations at which an auxiliary
outcome may be missing; default is 0.1{p_end}
{synopt : } - outcomes missing more often are not used; in panels, a location
counts as missing only if the outcome is missing in every period{p_end}
{synopt : } - auxiliary outcomes with no variation left after the regression
are also dropped{p_end}

{dlgtab:Panel setting}

{synopt :{opth t:imevar(varname)}} time identifier variable; required for panel
data. {cmd:cmd()} must then cluster by location or by groups of locations{p_end}

{dlgtab:Distance-based settings}

{synopt :{opth lat:itude(varname)}} latitude variable in signed decimal degrees
{p_end}

{synopt :{opth lon:gitude(varname)}} longitude variable in signed decimal
degrees {p_end}

{synopt :{opt distthres:hold(#)}} pairs of locations within {it:#} kilometers
({it:#} miles with {cmd:miles}) are always kept {p_end}
{synopt : } - requires setting {it:latitude} and {it:longitude} {p_end}

{dlgtab:Saving figures and tables}

{synopt :{opt file:suffix(str)}} folder path and base filename for saving
figures and results {p_end}
{synopt : } - required for the plot and save options below{p_end}

{synopt :{opt savedyad}} save {it:filesuffix}_dyad.dta, with each location
pair's correlation and contribution to the variance {p_end}

{synopt :{opt plotq}} plot the target variance V0, the penalty xi, and V0-xi
against {it:eta}, marking the selected {it:eta}. Saves {it:filesuffix}_q.pdf {p_end}

{synopt :{opt plothist}} histogram of the eligible pairs' correlations, with a
Gaussian fit matched to their interquartile range. Saves _hist.pdf and
_hist.png {p_end}
{synopt :{opt plothistnbins(#)}} number of bins for histogram of correlations
(default 1000) {p_end}

{synopt :{opt plotse}} plot the TMO SE at each {it:eta} divided by the original
SE, marking the selected {it:eta}. Saves _se.pdf {p_end}

{synopt :{opt saveest}} save _est.dta, _path.dta, and _anchor_path.dta {p_end}

{dlgtab:Threshold selection}

{pstd}
The threshold is chosen by maximizing the outcome of interest's thresholded
variance minus a penalty built from the auxiliary outcomes' variances. The
candidate thresholds are indexed by {it:eta} in (0,1] and depend on the data
and a pilot variance.{p_end}

{synopt :{opt eta(numlist)}} list of eta to search over ending exactly at 1.
Default 0.1, 0.15, ..., 1 {p_end}

{synopt :{opt maxeligiblepct(#)}} consider only thresholds that keep at most
{it:#} percent of eligible pairs; default 5; {cmd:maxeligiblepct(100)} means no
cap {p_end}

{synopt :{opt pilot(#)}} the pilot variance. Default is the square of the
standard error the command in {cmd:cmd()} estimates for {cmd:x()} {p_end}
{synopt :{opt pcl:ust(varname)}} use a variance clustered on {it:varname} as
the pilot instead {p_end}

{synopt :{opt threads(#)}} number of threads (default 4){p_end}

{dlgtab:Custom threshold}

{synopt :{opt thres:hold(#)}} set a custom threshold instead of using a rule
{p_end}

{synopt :{opt thresholdoff}} turns off the {cmd:tmo} adjustment entirely {p_end}

{dlgtab:SCPC options}

{synopt :{opt scpc}} construct a compatible regression for
{help scpc} from {cmd:cmd()}. Supports OLS and 2SLS without weights; fixed
effects must be categorical. {p_end}
{synopt : } - requires setting {it:latitude} and {it:longitude} {p_end}

{synopt :{opt scpc_uncond}} with {cmd:scpc}, use unconditional SCPC inference
{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
This Stata package implements the Thresholding Multiple Outcomes (TMO) method
for estimating standard errors in {help tmo##paper:DellaVigna et al. (2026)}.
The TMO method accounts for spatial correlation between locations by using a
set of auxiliary outcomes to estimate the correlation of errors between
locations. TMO keeps a pair of locations, allowing their errors to be
correlated, when the absolute value of their estimated correlation is at least
the optimal threshold. Pairs in the same {cmd:cmd()} cluster or within
{cmd:distthreshold()} are always kept ({it:protected}); all other pairs are
{it:eligible}.

{marker examples}{...}
{title:Examples}

{pstd}Load US county dataset{p_end}
{phang2}{cmd:. use "https://raw.githubusercontent.com/wjnkim/tmo/main/example/county_differences.dta", clear}{p_end}

{pstd}Define list of auxiliary outcomes{p_end}
{phang2}{cmd:. qui ds fips stfips PIN_persincpc_d EDU_college_d, not}{p_end}
{phang2}{cmd:. local ylist `r(varlist)'}{p_end}

{pstd}Regression of interest: change in per capita income on change in college
educated with state fixed effects{p_end}
{phang2}{cmd:. reg PIN_persincpc_d EDU_college_d i.stfips, r}{p_end}

{pstd}Run TMO{p_end}
{phang2}{cmd:. tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips)}{p_end}

{pstd}Run TMO and save figures{p_end}
{phang2}{cmd:. tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips) file(./example) plotq plothist plotse}{p_end}
{hline}

{pstd}Panel example{p_end}
{phang2}{cmd:. use "https://raw.githubusercontent.com/wjnkim/tmo/main/example/county_panel.dta", clear}{p_end}
{phang2}{cmd:. qui ds fips stfips EMN_farm EDU_publicenroll year, not}{p_end}
{phang2}{cmd:. local ylist `r(varlist)'}{p_end}
{phang2}{cmd:. tmo, cmd(reghdfe EMN_farm EDU_publicenroll i.year, absorb(stfips) cluster(fips)) x(EDU_publicenroll) ylist(`ylist') i(fips) t(year) file(./example_panel) plotq plothist}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:tmo} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(beta)}} coefficient on {cmd:x()}{p_end}
{synopt:{cmd:r(tmo_se)}} TMO standard error estimate{p_end}
{synopt:{cmd:r(orig_se)}} original standard error from {cmd:cmd()}{p_end}
{synopt:{cmd:r(lb)}} lower bound of 95% confidence interval{p_end}
{synopt:{cmd:r(ub)}} upper bound of 95% confidence interval{p_end}
{synopt:{cmd:r(threshold)}} selected threshold ({it:delta} in the output){p_end}
{synopt:{cmd:r(pct_ge_thres)}} percentage of off-diagonal location pairs kept,
including protected pairs {p_end}
{synopt:{cmd:r(pct_ge_thres_nocl)}} percentage of eligible pairs kept {p_end}
{synopt:{cmd:r(penalty)}} selected xi {p_end}
{synopt:{cmd:r(objective)}} selected V0-xi {p_end}
{synopt:{cmd:r(max_eligible_pct)}} eligible-pair percentage cap {p_end}
{synopt: }{it:From the selection step:}{p_end}
{synopt:{cmd:r(eta)}} the {it:eta} at the selected threshold{p_end}
{synopt:{cmd:r(pilot)}} pilot variance used{p_end}
{synopt:{cmd:r(N_used)}} auxiliary outcomes used{p_end}
{synopt:{cmd:r(N_eligible_pairs)}} number of eligible off-diagonal location pairs{p_end}
{synopt: }{it:Sample:}{p_end}
{synopt:{cmd:r(N)}} number of observations{p_end}
{synopt:{cmd:r(N_loc)}} number of locations{p_end}
{synopt:{cmd:r(T)}} number of time periods{p_end}
{synopt:{cmd:r(N_clust)}} number of clusters in {cmd:cmd()}, if any{p_end}
{synopt:{cmd:r(N_outcomes)}} number of auxiliary outcomes in {cmd:ylist()}{p_end}
{synopt: }{it:With {cmd:scpc}:}{p_end}
{synopt:{cmd:r(scpc_se)}} SCPC standard error{p_end}
{synopt:{cmd:r(tmo_scpc_se)}} TMO-SCPC hybrid standard error{p_end}
{synopt:{cmd:r(scpc_cv)}} SCPC critical value{p_end}

{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(etagrid)}} candidate {it:eta} values{p_end}
{synopt:{cmd:r(dgrid)}} candidate thresholds{p_end}
{synopt:{cmd:r(vpath)}} target variance V0 at each candidate{p_end}
{synopt:{cmd:r(xipath)}} penalty xi at each candidate{p_end}
{synopt:{cmd:r(objectivepath)}} objective V0-xi at each candidate{p_end}
{synopt:{cmd:r(eligiblepctpath)}} percent of eligible pairs kept at each
candidate; returned when {cmd:maxeligiblepct()} is below 100, as by default{p_end}
{p2colreset}{...}


{title:Dependencies}

{pstd}
{cmd:gtools} is required:{p_end}
{phang2}{cmd:. ssc install gtools}{p_end}

{pstd}
{help reghdfe}, {help ivreghdfe}, and {help ivreg2} are needed only for those
estimators. {cmd:tmo} was tested with versions 6.13.1, 1.1.4, and 4.1.12;
earlier versions may fail.{p_end}

{pstd}
The {cmd:scpc} option requires a TMO-compatible copy of
{browse "https://github.com/ukmueller/SCPC":SCPC} (the original will not work).
Install it with:{p_end}
{phang2}{cmd:. net install scpc, from("https://raw.githubusercontent.com/wjnkim/tmo/main/scpc_tmo") replace}{p_end}

{marker support}{...}
{title:Support}

{pstd}
Please report any bugs to wjnkim@stanford.edu. {p_end}

{marker references}{...}
{title:References}

{marker paper}{...}
{phang}DellaVigna, Stefano, Guido Imbens, Woojin Kim, and David Ritzwoller.
(2026). "Using Multiple Outcomes to Adjust Standard Errors for Spatial
Correlation." {browse "https://arxiv.org/abs/2504.13295"}{p_end}

{marker scpcpaper}{...}
{phang}Müller, Ulrich K., and Mark W. Watson. (2022). "Spatial Correlation
Robust Inference." {it:Econometrica} 90: 2901–2935.
{browse "https://www.princeton.edu/~umueller/SHAR.pdf"}{p_end}
