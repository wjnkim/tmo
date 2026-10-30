*! version 1.0.0 2026-10-29
* Threshold selection maximizes the target variance minus the auxiliary-outcome penalty.

* Any helpers-only reload must skip this block: tmo, tmo_work and
* tmo_session_work are all on the call stack when one fires.
if !inlist("`1'","__tmo_helpers_only","__tmo_scpc_backend_helpers_only","__tmo_scpc_auto_helpers_only") {
capture program drop tmo
program define tmo, rclass
    version 13
    capture noisily tmo_session_work `0'
    local tmo_rc = _rc
    if `tmo_rc' {
        return add
        exit `tmo_rc'
    }
    return add
end

capture program drop tmo_work
program define tmo_work, rclass
    version 13

    local raw_options `"`0'"'
    syntax, ///
        cmd(str) x(varname) ylist(varlist) Idvar(varname) ///
        [Timevar(varname)] ///
        [LATitude(varname)] [LONgitude(varname)] [DISTTHREShold(real 0)] [miles] ///
        [THREShold(real -9)] [thresholdoff] ///
        [MISSlimit(real 0.1)] ///
        [FILEsuffix(str)] ///
        [savedyad] [load(str)] ///
        [plotq] ///
        [plothist] [plothistnbins(int 1000)] ///
        [plotse] [saveplotseest] ///
        [saveest] [DIAGnostics] ///
        [SCPC] [scpc_cmd(str)] [scpc_uncond] ///
        [ETA(numlist ascending >0 <=1 min=2)] [PILOT(real -9)] ///
        [MAXELIGIBLEPCT(real -9)] ///
        [PCLust(varname)] [NONORMalize] [THREADS(int 4)]

    local automatic = (`threshold'==-9 & "`thresholdoff'"=="")
    * -9 is the not-supplied sentinel: automatic selection defaults to 5%.
    if `maxeligiblepct'==-9 local maxeligiblepct = cond(`automatic',5,100)
    if missing(`maxeligiblepct') | `maxeligiblepct'<0 | `maxeligiblepct'>100 {
        di as error "tmo: maxeligiblepct() must be between 0 and 100, inclusive."
        exit 198
    }
    if !`automatic' & `maxeligiblepct'<100 {
        di as error "tmo: maxeligiblepct() applies only to automatic threshold selection."
        exit 198
    }
    scalar __tmo_outcome_normalized = ("`nonormalize'"=="")
    local scpc_auto = ("`scpc'"!="")
    if `scpc_auto' & "`scpc_cmd'"!="" {
        di as error "tmo: scpc and scpc_cmd() are alternatives."
        exit 198
    }
    local scpc_requested = (`scpc_auto' | "`scpc_cmd'"!="")
    if !`scpc_requested' & "`scpc_uncond'"!="" {
        di as error "tmo: scpc_uncond requires scpc or scpc_cmd()."
        exit 198
    }
    if `automatic' {
        * syntax's real option normalizes its text and can lose the final digit.
        * Preserve the validated original distance literal for exact <= ties.
        tmo_distance_literal `raw_options'
        if `"`r(distance)'"'!="" local distthreshold `"`r(distance)'"'
    }
    if `automatic' & "`load'"!="" {
        di as txt "tmo: automatic selection recomputes its inputs; load() is not used."
        local load
    }
    local outopts `diagnostics' `savedyad' `saveest' `plotq' `plotse' `plothist' `saveplotseest'
    capture which gtools
    if _rc {
        di as error "tmo requires gtools package -- please run: ssc install gtools"
        exit 499
    }
    scalar __tmo_s_tmo_scpc_se = .
    scalar __tmo_s_scpc_se = .
    scalar __tmo_version = "1.0.0"
    scalar __tmo_nnovar = 0
    scalar __tmo_ndepend = .
    if "`filesuffix'"!="" local savepath filesuffix(`"`filesuffix'"')

    * Resolve both source-checkout and net-installed Java kernels.
    quietly findfile tmo.ado
    local jarpath = subinstr(`"`r(fn)'"', char(92), "/", .)
    local jarpath = regexr(`"`jarpath'"', "/[^/]+$", "")
    capture confirm file `"`jarpath'/tmo.jar"'
    if _rc {
        capture findfile tmo.jar
        if !_rc {
            local jarpath = subinstr(`"`r(fn)'"', char(92), "/", .)
            local jarpath = regexr(`"`jarpath'"', "/[^/]+$", "")
        }
    }


    qui `cmd'
    local spec `e(cmd)'

    if !inlist("`spec'","regress","reghdfe","ivreghdfe","ivreg2") {
            di as error "`spec' not supported"
            exit 198
    }

    scalar __tmo_s_finite_k = .
    scalar __tmo_s_dof_adj = 1
    if `automatic' {
        tmo_finite_sample_adjustment, spec("`spec'")
        scalar __tmo_s_finite_k = r(k)
        scalar __tmo_s_dof_adj = r(factor)
    }

    local y `e(depvar)'
    tmo_outcome_controls, outcomes(`y' `ylist')
    scalar __tmo_s_beta = _b[`x']
    scalar __tmo_s_se = _se[`x']
    scalar __tmo_s_N_obs = e(N)
    scalar __tmo_s_N_clust = e(N_clust)
    scalar __tmo_s_df_r = e(df_r)
    tempvar tmo_sample
    gen byte `tmo_sample' = e(sample)
    quietly count if `tmo_sample'
    scalar __tmo_s_N_rows = r(N)

    * All numerical paths share one location code; the caller's ID is unchanged.
    local source_idvar `idvar'
    quietly count if `tmo_sample' & missing(`idvar')
    if r(N) {
        di as error "tmo: idvar() contains missing IDs in the main estimation sample."
        exit 459
    }
    local stringid
    capture confirm string variable `idvar'
    if !_rc {
        tempvar tmo_numeric_id tmo_string_id
        local group_source `idvar'
        if "`: type `idvar''"=="strL" {
            quietly clonevar `tmo_string_id' = `idvar'
            quietly compress `tmo_string_id', nocoalesce
            local group_source `tmo_string_id'
        }
        quietly gegen long `tmo_numeric_id' = group(`group_source') if `tmo_sample'
        assert !missing(`tmo_numeric_id') & `tmo_numeric_id'>=1 if `tmo_sample'
        local idvar `tmo_numeric_id'
        local stringid stringid
    }

    local cluster `e(clustvar)'

    local absorb_vars `e(absvars)'

    * Use the fitted design, including native omissions and expanded names.
    local native_batch 0
    local native_rhs
    local native_controls
    local native_nocons 1
    tempvar native_resid
    if inlist("`spec'","regress","reghdfe") {
        tempname native_b native_omit
        matrix `native_b' = e(b)
        _ms_omit_info `native_b'
            matrix `native_omit' = r(omit)
            local native_names : colnames `native_b'
            local native_j = 1
            local native_target 0
            foreach native_var of local native_names {
                if `native_omit'[1,`native_j']==0 {
                    if "`native_var'"=="_cons" local native_nocons 0
                    else {
                        local native_rhs `native_rhs' `native_var'
                        if "`native_var'"=="`x'" local native_target 1
                        else local native_controls `native_controls' `native_var'
                    }
                }
                local ++native_j
            }
            if !`native_target' {
                di as error "tmo: x() must identify a retained main-effect coefficient."
                exit 198
            }
            if "`spec'"=="regress" {
                quietly predict double `native_resid' if `tmo_sample', residuals
                local native_batch 1
            }
    }

    if "`e(wexp)'"!="" {
        local weight_expr = subinstr(`"`e(wexp)'"',"=","",1)
        local weight_type `e(wtype)'
        tempvar tmo_weight
        quietly generate double `tmo_weight' = `weight_expr' if `tmo_sample'
        local weightvar `tmo_weight'
        local weightexp [`weight_type'=`weightvar']
    }
    else {
        local weightvar
        local weightexp
    }





    if `scpc_requested' {
        global scpc_cv = .
        if "`longitude'"=="" | "`latitude'"=="" {
            di as error "longitude() and latitude() required for SCPC"
            exit 198
        }
        if "`weightvar'"!="" {
            di as error "weights not allowed for scpc"
            exit 198
        }
        capture which scpc
        if _rc {
            di as error "SCPC requires the compatible companion package; see help tmo for installation."
            exit 499
        }

        preserve
          if `scpc_auto' {
            local scpc_modeopt
            if "`scpc_uncond'"!="" local scpc_modeopt uncond
            quietly tmo_scpc_auto, cmd(`"`cmd'"') x(`x') sample(`tmo_sample') ///
                idvar(`idvar') latitude(`latitude') longitude(`longitude') `scpc_modeopt'
            local scpc_cluster `r(cluster)'
            scalar __tmo_s_scpc_se = r(scpc_se)
          }
          else {
            qui keep if !missing(`longitude') & !missing(`latitude') & `tmo_sample'

            rename `latitude' s_1
            rename `longitude' s_2
            order s_1 s_2

            qui hashsort `idvar'
            qui `scpc_cmd'
            qui keep if e(sample)
            local scpc_cluster `e(clustvar)'
            local scpc_target = colnumb(e(b),"`x'")
            if missing(`scpc_target') {
                di as error "tmo: scpc_cmd() must estimate the coefficient named in x()."
                exit 198
            }
            scalar __tmo_scpc_owned = 1
            tmo_scpc_backend
            if "`scpc_uncond'"=="" {
                qui scpc, k(`scpc_target') latlong
            }
            else {
                qui scpc, k(`scpc_target') latlong uncond
            }
            scalar __tmo_s_scpc_se=e(scpcstats)[`scpc_target',2]

            * Some SCPC distributions clear Mata when first loaded.
            capture mata: tmo_lowtri_gram(J(0,0,.))
            if _rc {
                quietly findfile tmo.ado
                quietly run `"`r(fn)'"' __tmo_helpers_only
            }

            tempvar scpc_group
            if "`scpc_cluster'"!="" {
                quietly gegen long `scpc_group' = group(`scpc_cluster'), missing
            }
            else quietly generate long `scpc_group' = _n
            mata: tmo_scpc_capture("`idvar'","`scpc_group'")
            tmo_scpc_release
          }
            if !`automatic' {
                clear
                mata: __tmo_legacy_Wfin_sum_vec = tmo_lowtri_gram(__tmo_scpc_basis) :/ cols(__tmo_scpc_basis)
                mata: __tmo_legacy_id_scpc_rowvec = __tmo_legacy_id_scpc_colvec = J(0,0,.)
                mata: lowtri_ids(__tmo_scpc_ids, __tmo_legacy_id_scpc_rowvec, __tmo_legacy_id_scpc_colvec)
                mata: st_local("scpc_obsN",strofreal(rows(__tmo_legacy_Wfin_sum_vec),"%50.0f"))
                gen double id1=.
                gen double id2=.
                gen double Wfin=.
                qui set obs `scpc_obsN'
                mata: st_store(.,.,(__tmo_legacy_id_scpc_rowvec,__tmo_legacy_id_scpc_colvec,__tmo_legacy_Wfin_sum_vec))
                qui compress
                tempfile Wfin
                qui save `Wfin'
            }
        restore
    }
    else {
        global scpc_cv = .
    }


    * Only exports and cache validation need the original location labels.
    * Capture after SCPC's first load, which can clear Mata session objects.
    if "`stringid'"!="" & "`savedyad'`load'"!="" {
        preserve
            quietly keep if `tmo_sample'
            keep `idvar' `source_idvar'
            quietly hashsort `idvar'
            quietly by `idvar': keep if _n==1
            assert `idvar'==_n
            mata: __tmo_idlabels = st_sdata(.,st_local("source_idvar"))
            mata: assert(all(__tmo_idlabels:!=""))
        restore
    }




    * The fitted coefficient identifies x; check only the dependent-variable token.
    local ycheck: word 2 of `cmd'
	if "`y'"!="`ycheck'" {
		di as error "cmd() must name `y' as its dependent variable"
		exit 198
	}

    local done ""
    local dups ""
    foreach var in `y' `ylist' {
        confirm var `var', exact
        if strpos(" `done' "," `var' ")>0 {
            local dups "`dups' `var'"
        }
        local done "`done' `var'"
    }
    if "`dups'"!="" {
        di as error "Duplicated variables in depvar/ylist: `dups'"
        exit 198
    }

    if `misslimit'<0 | `misslimit'>1 {
        di as error "misslimit() value must be between 0 and 1"
        exit 198
    }

    if `plothistnbins'<=0 {
        di as error "plothistnbins() value must be positive"
        exit 198
    }

    if "`plotse'"=="" & "`saveplotseest'"!="" {
        di as error "saveplotseest option requires plotse option"
        exit 198
    }

    if "`filesuffix'"=="" & ("`plotq'`plotse'`plothist'`savedyad'`saveplotseest'`saveest'"!="") {
        di as error "filesuffix() required for `plotse' `plothist' `savedyad' `saveplotseest' `saveest'"
        exit 198
    }

    if "`timevar'"!="" & "`cluster'"=="" {
        di as error "Clustering required for panel case (at least at `idvar' level)"
        exit 198
    }

    if inlist("`spec'","reghdfe","ivreghdfe") & "`absorb_vars'"=="" {
        di as error "absorb() required for `spec'"
        exit 198
    }

    confirm numeric variable `idvar'

	if "`idvar'"!="" & "`timevar'"=="" {
		gisid `idvar' if `tmo_sample'
	}
	if "`idvar'"!="" & "`timevar'"!="" {
		gisid `idvar' `timevar' if `tmo_sample'
	}

    if ("`longitude'"=="" | "`latitude'"=="") & `distthreshold'!=0 {
        di as error "longitude() and latitude() required for distthreshold() option"
        exit 198
    }

    if ("`longitude'"=="" | "`latitude'"=="") & `scpc_requested' {
        di as error "longitude() and latitude() required for SCPC"
        exit 198
    }

    if `scpc_requested' & "`weightvar'"!="" {
        di as error "weights not allowed for scpc"
        exit 198
    }

    * Zero disables distance protection; missing cannot be used as infinity.
    if missing(`distthreshold') | `distthreshold'<0 {
        di as error "distthreshold() must be finite and nonnegative"
        exit 198
    }

    if !`automatic' & "`longitude'"!="" & "`latitude'"!="" & (`distthreshold'>0 | "`savedyad'"!="") {
        cap which geodist
        if _rc {
            di as error "legacy distance dyads require geodist -- please run: ssc install geodist"
            exit 499
        }
    }

    * A missing cutoff is not a substitute for thresholdoff.
    if `threshold'!=-9 & (missing(`threshold') | `threshold'<0) {
        di as error "Custom threshold() must be finite and nonnegative"
        exit 198
    }

    qui gdistinct `idvar' if `tmo_sample'
    local N = r(ndistinct)
    scalar __tmo_s_N = `r(ndistinct)'
    mata: __tmo_legacy_N = `N'
    if `scpc_requested' {
        scalar __tmo_scpc_normN = scalar(__tmo_s_N)
        if "`scpc_cluster'"!="" & "`scpc_cluster'"!="`idvar'" {
            quietly gdistinct `scpc_cluster' if `tmo_sample'
            scalar __tmo_scpc_normN = r(ndistinct)
        }
        if missing(scalar(__tmo_scpc_normN)) | scalar(__tmo_scpc_normN)<scalar(__tmo_scpc_neff) {
            di as error "tmo: SCPC normalization count is incompatible with its estimation sample."
            exit 459
        }
    }
    if "`timevar'"!="" {
        qui gdistinct `timevar' if `tmo_sample'
        local T = r(ndistinct)
        scalar __tmo_s_T = r(ndistinct)
        mata: __tmo_legacy_T = `T'
    }
    else {
        local T = 1
        scalar __tmo_s_T = 1
        mata: __tmo_legacy_T = 1
    }





    if inlist("`spec'","reghdfe","ivreghdfe") {
        mata: st_local("comma_start", strofreal(tmo_optcomma(st_local("cmd"))))
        local before_comma = substr("`cmd'", 1, `comma_start'-1)
        local after_comma = substr("`cmd'", `comma_start'+1, .)
        tmo_fe_options, text(`"`after_comma'"')
        local after_comma `"`r(options)'"'

        local cmd `before_comma', absorb(`absorb_vars') `after_comma'
    }





    * Parse token boundaries and the option comma outside parentheses/quotes;
    * legal if expressions and similarly named variables must remain intact.
    mata: st_local("ccomma", strofreal(tmo_optcomma(st_local("cmd"))))

    if `ccomma'>0 {
        local cmdhead = substr(`"`cmd'"', 1, `ccomma'-1)
        local cmdopts = substr(`"`cmd'"', `ccomma'+1, .)
    }
    else {
        local cmdhead `"`cmd'"'
        local cmdopts
    }

    gettoken cmdname cmdrest : cmdhead
    gettoken cmddep  cmdrhs  : cmdrest

    local cmdvars
    local cmdqual
    local inqual 0
    foreach tk of local cmdrhs {
        if !`inqual' & (inlist("`tk'","if","in") | substr("`tk'",1,1)=="[") local inqual 1
        if `inqual' local cmdqual `cmdqual' `tk'
        else        local cmdvars `cmdvars' `tk'
    }

    * nocons means the partialled-out regressor is not mean zero, so the
    * verification regression further down has to suppress the constant too
    local orignocons
    if regexm("`cmdopts'", "(^|[ ,])noc(o(n(s(t(a(n(t)?)?)?)?)?)?)?([ ,]|$)") ///
        local orignocons , nocons


    * Compute the clustered pilot on the original estimation sample.
    tempname clpilot
    scalar `clpilot' = .
    if "`pclust'"!="" & (`threshold'!=-9 | "`thresholdoff'"!="") {
        di as txt "tmo: pclust() is ignored because threshold()/thresholdoff bypass threshold selection."
    }
    if "`pclust'"!="" & `threshold'==-9 & "`thresholdoff'"=="" {
        if `pilot'!=-9 {
            di as error "tmo: pilot() and pclust() are alternatives"
            exit 198
        }
        qui count if `tmo_sample' & missing(`pclust')
        if r(N)>0 {
            di as error "tmo: pclust() must be nonmissing on the estimation sample."
            exit 459
        }
        qui gdistinct `pclust' if `tmo_sample'
        if missing(r(ndistinct)) | r(ndistinct)<2 {
            di as error "tmo: pclust() requires at least two clusters on the estimation sample."
            exit 459
        }
        tmo_cluster_options, `cmdopts'
        local pilotopts `"`r(options)'"'
        if inlist("`spec'","regress","reghdfe") local pilotvce vce(cluster `pclust')
        else local pilotvce cluster(`pclust')
        qui `cmdname' `cmddep' `cmdvars' if `tmo_sample' `weightexp', `pilotopts' `pilotvce'
        capture assert e(sample)==`tmo_sample'
        if _rc | reldif(_b[`x'],scalar(__tmo_s_beta))>1e-8 {
            di as error "tmo: pclust() changed the estimation sample or coefficient; supply pilot() instead."
            exit 459
        }
        scalar `clpilot' = _se[`x']^2
        if missing(scalar(`clpilot')) | scalar(`clpilot')<=0 {
            di as error "tmo: pclust() produced a missing or nonpositive pilot variance; supply pilot() instead."
            exit 459
        }
    }

    * The original if/in restrictions are already recorded in the sample marker.
    local cmdqual `weightexp'
    local cmd `cmdname' `cmddep' `cmdvars' `cmdqual'
    if `ccomma'>0 local cmd `cmd', `cmdopts'



    local retention diagonal_only
    if `automatic' {
        mata: __tmo_Bcl = __tmo_Bcoords = J(0,0,.)
        mata: __tmo_Bclids = __tmo_Bcoordids = J(0,1,.)
        mata: __tmo_Bpairs = J(0,2,.)
    }
    if "`cluster'"!="" & "`cluster'"!="`idvar'" {
        preserve
            qui keep if `tmo_sample'
            keep `idvar' `cluster'

            qui gduplicates drop
            cap gisid `idvar'

            if _rc {
                di as error "Only clustering by groups of locations is supported (`cluster' must be constant within `idvar')"
                exit 459
            }

            * Cluster codes preserve equality within each dimension, including
            * missing levels, which compare equal for pair protection.
            qui hashsort `idvar'
            local clvars
            foreach var in `cluster' {
                tempvar cluster_group
                qui gegen double `cluster_group' = group(`var')
                local clvars `clvars' `cluster_group'
            }
            mata: __tmo_clmat = st_data(., "`clvars'")
            if `automatic' {
                mata: __tmo_Bcl = __tmo_clmat
                mata: __tmo_Bclids = st_data(., "`idvar'")
                local retention cluster
            }
        restore
    }





    * Automatic selection needs only one coordinate row per observed location.
    * It does not build a distance dyad dataset merely to construct protection.
    if `automatic' & `distthreshold'>0 {
        if "`retention'"=="cluster" local retention cluster_distance
        else local retention distance
        preserve
            quietly keep if `tmo_sample'
            keep `idvar' `latitude' `longitude'
            quietly gduplicates drop
            capture gisid `idvar'
            if _rc {
                di as error "tmo: coordinates must be constant within each location."
                exit 459
            }
            capture assert inrange(`latitude',-90,90) | missing(`latitude')
            local badcoords = _rc
            capture assert inrange(`longitude',-180,180) | missing(`longitude')
            if `badcoords' | _rc {
                di as error "tmo: latitude must be in [-90,90] and longitude in [-180,180]."
                exit 459
            }
            quietly count if missing(`latitude',`longitude')
            if r(N)>0 di as txt "tmo: missing coordinates do not establish distance protection; cluster protection still applies."
            quietly hashsort `idvar'
            assert `idvar'>`idvar'[_n-1] if _n>1
            mata: __tmo_Bcoords = st_data(.,"`latitude' `longitude'")
            mata: __tmo_Bcoordids = st_data(.,"`idvar'")
        restore
    }

    * Only override selection and dyad exports need distance dyads;
    * SCPC uses coordinates directly.
    if "`longitude'"!="" & "`latitude'"!="" & !`automatic' & (`distthreshold'>0 | "`savedyad'"!="") {
        preserve
            qui keep if `tmo_sample'
            keep `idvar' `longitude' `latitude'

            qui gduplicates drop
            cap gisid `idvar'

            if _rc {
                di as error "`longitude' and/or `latitude' not constant within some `idvar'"
                exit 459
            }

            * Construct the lower triangle without an N x N intermediate.
            qui hashsort `idvar'
            mata: __dst_id = st_data(.,"`idvar'")
            mata: __dst_la = st_data(.,"`latitude'")
            mata: __dst_lo = st_data(.,"`longitude'")

            clear
            mata: __dst_i1 = __dst_i2 = J(0,0,.)
            mata: lowtri_ids(__dst_id, __dst_i1, __dst_i2)
            mata: __dst_la1 = __dst_la2 = J(0,0,.)
            mata: lowtri_ids(__dst_la, __dst_la1, __dst_la2)
            mata: __dst_lo1 = __dst_lo2 = J(0,0,.)
            mata: lowtri_ids(__dst_lo, __dst_lo1, __dst_lo2)
            mata: st_local("npair",strofreal(rows(__dst_i1),"%50.0f"))

            gen double id1=.
            gen double id2=.
            gen double lat1=.
            gen double lon1=.
            gen double lat2=.
            gen double lon2=.

            qui set obs `npair'
            mata: st_store(.,.,(__dst_i1,__dst_i2,__dst_la1,__dst_lo1,__dst_la2,__dst_lo2))
            mata: mata drop __dst_id __dst_la __dst_lo __dst_i1 __dst_i2 __dst_la1 __dst_la2 __dst_lo1 __dst_lo2

            qui geodist lat1 lon1 lat2 lon2, gen(dist) `miles' sphere
            qui replace dist=0 if id1==id2

            keep id1 id2 dist

            qui compress
            tempfile dist
            qui save `dist'
        restore
    }





    tempvar tmo_xtilde tmo_xhat tmo_xb
    if  inlist("`spec'","regress","reghdfe") {
        * Materialize the retained native design, not the raw command tokens.
        local xtildevars
        if `"`native_controls'"'!="" {
            local native_expand
            foreach term of local native_controls {
                local expanded
                local pieces : subinstr local term "#" " ", all
                foreach piece of local pieces {
                    * A base level is written 0b. or 0bn., so the marker has to
                    * be optional here. Without it the base term keeps its own
                    * base declaration while the other levels are rewritten to
                    * bn., and fvrevar rejects the pair as a base conflict.
                    if regexm("`piece'","^([0-9]+)b?n?\.(.+)$") {
                        local piece `=regexs(1)'bn.`=regexs(2)'
                    }
                    if "`expanded'"=="" local expanded `piece'
                    else local expanded `expanded'#`piece'
                }
                local native_expand `native_expand' `expanded'
            }
            quietly fvrevar `native_expand'
            local xtildevars `r(varlist)'
        }
        local xtildecmd `cmdname' `x' `xtildevars' `cmdqual'
        if `ccomma'>0 {
            local xtildecmd `xtildecmd', `cmdopts'
        }
        if "`spec'"=="reghdfe" {
            local xtildecmd `xtildecmd' residuals(`tmo_xtilde')
        }
    }
    if inlist("`spec'","ivreghdfe","ivreg2") {
        * Bind parenthesized tokens so factor groups and IV clauses stay intact.
        local iv_rest `"`cmdvars'"'
        local iv_block
        local controls
        while strtrim(`"`iv_rest'"')!="" {
            gettoken iv_token iv_rest : iv_rest, bind
            if substr(strtrim(`"`iv_token'"'),1,1)=="(" & strpos(`"`iv_token'"',"=")>0 {
                if `"`iv_block'"'!="" {
                    di as error "Multiple endogenous-regressor clauses are not supported"
                    exit 198
                }
                local iv_token = strtrim(`"`iv_token'"')
                local iv_block = substr(`"`iv_token'"',2,strlen(`"`iv_token'"')-2)
            }
            else local controls `controls' `iv_token'
        }
        local instr_start = strpos(`"`iv_block'"',"=")
        if `instr_start'==0 {
            di as error "cmd() must contain an (endogenous = instruments) clause"
            exit 198
        }
        local endog = strtrim(substr(`"`iv_block'"',1,`instr_start'-1))
        if wordcount(`"`endog'"')>1 {
            di as error "Multiple endogenous regressors not supported"
            exit 198
        }
        if "`endog'"!="`x'" {
            di as error "Endogenous regressor `endog' is not `x'"
            exit 198
        }
        local instr = substr(`"`iv_block'"',`instr_start'+1,.)
        if "`spec'"=="ivreghdfe" {
            local xtildecmd1 reghdfe `endog' `instr' `controls' `cmdqual', `cmdopts'
            local xtildecmd2 reghdfe `tmo_xhat' `controls' `cmdqual', `cmdopts' residuals(`tmo_xtilde')
        }
        else {
            local xtildecmd1 regress `endog' `instr' `controls' `cmdqual' `orignocons'
            local xtildecmd2 regress `tmo_xhat' `controls' `cmdqual' `orignocons'
        }
    }




    preserve
        * Cached dyads do not contain the per-anchor selection inputs.
        if "`load'"=="" | (`threshold'==-9 & "`thresholdoff'"=="") {
            qui keep if `tmo_sample'

            if  inlist("`spec'","regress","reghdfe") {
                qui `xtildecmd'
                if "`spec'"=="regress" qui predict double `tmo_xtilde', resid
            }
            if inlist("`spec'","ivreghdfe","ivreg2") {
                qui `xtildecmd1'
                qui predict double `tmo_xhat'
                qui `xtildecmd2'
                if "`spec'"=="ivreg2" qui predict double `tmo_xtilde', resid
                drop `tmo_xhat'
            }

            qui reg `y' `tmo_xtilde' `weightexp' `orignocons'
            if abs(scalar(__tmo_s_beta)-_b[`tmo_xtilde'])>1e-5 {
                di as error "The partialled-out regressor does not reproduce the original coefficient"
                exit 459
            }

            * Allocate residual columns in outcome order, with the main one first.
            local allresids
            foreach dv in `y' `ylist' {
                tempvar outcome_resid
                local allresids `allresids' `outcome_resid'
            }
            gettoken mainresid residlist : allresids
            * Each outcome's coefficient on x, for the dependence check below.
            tempname auxb
            matrix `auxb' = J(1,`: word count `y' `ylist'',.)

            if  inlist("`spec'","regress","ivreg2") {
                * Variance-only options do not affect residuals; nocons does.
                * Replace the outcome token without changing similarly named variables.
                local batchvars
                local plainols = ("`spec'"=="regress" & `native_batch')
                foreach rhs of local cmdvars {
                    capture confirm numeric variable `rhs'
                    if _rc local plainols 0
                }
                if `plainols' {
                    foreach dv in `y' `ylist' {
                        quietly count if missing(`dv')
                        if r(N)==0 local batchvars `batchvars' `dv'
                    }
                }
                local nbatch : word count `batchvars'
                if `nbatch'>1 {
                    local batchres
                    local bj = 1
                    foreach dv in `y' `ylist' {
                        if strpos(" `batchvars' "," `dv' ") {
                            local current_resid : word `bj' of `allresids'
                            quietly gen double `current_resid' = .
                            local batchres `batchres' `current_resid'
                        }
                        local ++bj
                    }
                    local batch_ok 0
                    mata: tmo_ols_batch("`native_rhs'", "`batchvars'", "`batchres'", "`weightvar'", `native_nocons', "`native_resid'", "`x'")
                    if !`batch_ok' {
                        drop `batchres'
                        local nbatch 0
                    }
                    else {
                        local bk = 1
                        local bj = 1
                        foreach dv in `y' `ylist' {
                            if strpos(" `batchvars' "," `dv' ") {
                                matrix `auxb'[1,`bj'] = __tmo_batch_bx[1,`bk']
                                local ++bk
                            }
                            local ++bj
                        }
                        capture matrix drop __tmo_batch_bx
                    }
                }
                local ynum=1
                foreach dv in `y' `ylist' {
                    local current_resid : word `ynum' of `allresids'
                    if `nbatch'<=1 | !strpos(" `batchvars' "," `dv' ") {
                        capture quietly `cmdname' `dv' `cmdvars' `cmdqual' `orignocons'
                        local residual_rc = _rc
                        if inlist(`residual_rc',2000,2001) & "`dv'"!="`y'" {
                            quietly count if !missing(`dv')
                            if r(N)>0 exit `residual_rc'
                            quietly generate double `current_resid' = .
                        }
                        else {
                            if `residual_rc' exit `residual_rc'
                            qui predict double `current_resid', resid
                            capture matrix `auxb'[1,`ynum'] = _b[`x']
                        }
                    }
                    local ++ynum
                }
            }

            if "`spec'"=="reghdfe" {
                tmo_cluster_options, `cmdopts'
                local loopopts `"`r(options)'"'

                local ynum=1
                foreach dv in `y' `ylist' {
                    local current_resid : word `ynum' of `allresids'
                    local cmd_inloop `cmdname' `dv' `cmdvars' `cmdqual', `loopopts'
                    capture quietly `cmd_inloop' residuals(`current_resid')
                    local residual_rc = _rc
                    if inlist(`residual_rc',2000,2001) & "`dv'"!="`y'" {
                        quietly count if !missing(`dv')
                        if r(N)>0 exit `residual_rc'
                        quietly generate double `current_resid' = .
                    }
                    else {
                        if `residual_rc' exit `residual_rc'
                        capture matrix `auxb'[1,`ynum'] = _b[`x']
                    }
                    local ++ynum
                }
            }

            if "`spec'"=="ivreghdfe" {
                tmo_cluster_options, `cmdopts'
                local loopopts `"`r(options)'"'

                local absorb_vars_savefe
                local absorb_vars_fesum
                local fevars
                foreach fe_var in `absorb_vars' {
                    tempvar saved_fe
                    local absorb_vars_savefe `absorb_vars_savefe' `saved_fe'=`fe_var'
                    local absorb_vars_fesum `absorb_vars_fesum' + `saved_fe'
                    local fevars `fevars' `saved_fe'
                }
                tmo_fe_options, text(`"`loopopts'"')
                local loopopts `"`r(options)' absorb(`absorb_vars_savefe')"'

                local ynum=1
                foreach dv in `y' `ylist' {
                    local current_resid : word `ynum' of `allresids'
                    local cmd_inloop `cmdname' `dv' `cmdvars' `cmdqual', `loopopts'
                    capture drop `fevars'
                    capture quietly `cmd_inloop'
                    local residual_rc = _rc
                    if inlist(`residual_rc',2000,2001) & "`dv'"!="`y'" {
                        quietly count if !missing(`dv')
                        if r(N)>0 exit `residual_rc'
                        quietly generate double `current_resid' = .
                    }
                    else {
                        if `residual_rc' exit `residual_rc'
                        capture matrix `auxb'[1,`ynum'] = _b[`x']
                        capture drop `tmo_xb'
                        qui predict double `tmo_xb'
                        qui gen double `current_resid' = `dv' - (`tmo_xb'`absorb_vars_fesum')
                    }
                    local ++ynum
                }
            }
            * -2 because of the main outcome of interest and ++ynum at end of loop
            scalar __tmo_s_D = `ynum' - 2

            * Drop auxiliaries that the regression leaves as rounding noise, and
            * check the rest for exact linear dependence.
            tempname novar
            local novar_resids
            mata: tmo_aux_checks("`y' `ylist'", "`allresids'", "`weightvar'", "`idvar'", "`tmo_xtilde'", `misslimit', "`auxb'", "`novar'")
            local residlist : list residlist - novar_resids
            if "`residlist'"=="" {
                di as error "tmo: no auxiliary outcome has variation after the regression."
                exit 459
            }

            keep `idvar' `timevar' `weightvar' `pclust' `tmo_xtilde' `allresids'
            qui hashsort `idvar' `timevar'
            quietly gisid `idvar' `timevar'
            assert `idvar'>=`idvar'[_n-1] if _n>1

            * Capture only observed estimation rows, before panel filling.
            assert !missing(`mainresid', `tmo_xtilde')
            if "`weightvar'"!="" assert !missing(`weightvar') & `weightvar'>0
            mata: __tmo_input_R = st_data(.,"`allresids'")
            mata: __tmo_input_x = st_data(.,"`tmo_xtilde'")
            mata: __tmo_input_a = J(rows(__tmo_input_x),1,1)
            if "`weightvar'"!="" {
                mata: __tmo_input_a = st_data(.,"`weightvar'")
            }
            mata: __tmo_input_w = __tmo_input_x :* __tmo_input_a
            mata: __tmo_input_id = st_data(.,"`idvar'")
        }

        if `automatic' {
            mata: tmo_build_protection(`distthreshold', "`miles'"!="")
            scalar __tmo_need_counts = ///
                ("`diagnostics'`saveest'`savedyad'`plothist'`saveplotseest'"!="" | ///
                `maxeligiblepct'<100)
            if "`load'"!="" di as txt "tmo: automatic selection recomputes inputs; load() is not used."
            tmo_estimate, eta(`eta') pilot(`pilot') pclust(`pclust') ///
                clpilot(`clpilot') `nonormalize' misslimit(`misslimit') ///
                threads(`threads') libpath(`"`jarpath'"') ///
                dofadj(`=scalar(__tmo_s_dof_adj)') ///
                maxeligiblepct(`maxeligiblepct') protected novar(`novar')
            tmo_selection_capture
            * The SCPC hybrid screens with the selection's own pair scores.
            if `scpc_requested' {
                mata: tmo_scpc_stream(tmo_selection_corr("`selkernel'", st_local("jarpath")), ///
                    __tmo_input_id, __tmo_input_x, __tmo_input_R[,1], "`scpc_uncond'"!="")
            }
            tmo_selection_output, `outopts' `savepath' `stringid' nbins(`plothistnbins') ///
                xname(`x') kernel(`selkernel') pilotsource(`selpilot') libpath(`"`jarpath'"') retention(`retention')
            return add
            restore
            tmo_cleanup
            exit
        }

        if "`load'"=="" {
            * A custom threshold screens with automatic selection's statistic;
            * thresholdoff screens no pairs, so it builds no screen.
            local screen_nrm = cond("`nonormalize'"=="","normalize","")
            mata: __tmo_legacy_id = st_data(.,"`idvar'")
            mata: __tmo_legacy_xtilde = st_data(.,"`tmo_xtilde'")
            if "`weightvar'"!="" {
                mata: __tmo_legacy_wgt = st_data(.,"`weightvar'")
            }
            else {
                mata: __tmo_legacy_wgt = J(rows(__tmo_legacy_xtilde),1,1)
            }
            mata: __tmo_legacy_xtilde_wgt = __tmo_legacy_xtilde:*__tmo_legacy_wgt
            mata: __tmo_legacy_Res1_x = st_data(.,"`mainresid'")

            mata: __tmo_legacy_id_widerowvec = __tmo_legacy_id_widecolvec = J(0,0,.)
            mata: lowtri_ids(uniqrows(__tmo_legacy_id), __tmo_legacy_id_widerowvec, __tmo_legacy_id_widecolvec)
            if "`thresholdoff'"=="" {
                mata: tmo_screen_inputs("`screen_nrm'", `misslimit', "`novar'")
                mata: __tmo_legacy_CovEpsVec = lowtri_vec(tmo_target_screen(__tmo_Xs, st_numscalar("__tmo_hasmiss")))
                scalar __tmo_s_D_cols = scalar(__tmo_nout)
            }
            else {
                mata: __tmo_legacy_CovEpsVec = J(rows(__tmo_legacy_id_widerowvec),1,.)
                scalar __tmo_s_D_cols = .
            }

            mata: __tmo_legacy_ResXtildeVec = J(0,0,.)
            if "`timevar'"=="" {
                mata: sandwich_crosssec(__tmo_legacy_id_widerowvec, __tmo_legacy_id_widecolvec, __tmo_legacy_xtilde_wgt, __tmo_legacy_Res1_x, __tmo_legacy_ResXtildeVec)
            }
            else {
                mata: sandwich_panel_collapse(__tmo_legacy_N, __tmo_legacy_id, __tmo_legacy_xtilde_wgt, __tmo_legacy_Res1_x, __tmo_legacy_id_widerowvec, __tmo_legacy_id_widecolvec, __tmo_legacy_ResXtildeVec)
            }

            mata: resxtildenorm(__tmo_legacy_ResXtildeVec, __tmo_legacy_xtilde, __tmo_legacy_xtilde_wgt)

            if `scpc_requested' {
                mata: __tmo_legacy_y_scpc = (sqrt(st_numscalar("__tmo_scpc_normN"))/(__tmo_legacy_xtilde'*__tmo_legacy_xtilde)):*__tmo_legacy_Res1_x:*__tmo_legacy_xtilde
                mata: st_local("scpc_obsN",strofreal(rows(__tmo_legacy_id),"%50.0f"))

                clear
                gen double id1=.
                gen double y_scpc1=.

                qui set obs `scpc_obsN'

                mata: st_store(.,.,(__tmo_legacy_id,__tmo_legacy_y_scpc))

                cap gisid id1
                if _rc gcollapse (sum) y_scpc1, by(id1)

                qui compress
                tempfile scpc1
                qui save `scpc1'

                rename id1 id2
                rename y_scpc1 y_scpc2

                tempfile scpc2
                qui save `scpc2'
            }

            if "`cluster'"!="" & "`cluster'"!="`idvar'" {
                local clp samecl
            }
            else {
                local clp
            }
            if "`filesuffix'"!="" {
                local savepath filesuffix(`filesuffix')
            }
            else {
                local savepath
            }
            if "`dist'"!="" {
                local distp distpath(`dist')
            }
            else {
                local distp
            }
            if `scpc_requested' {
                local scpcw scpcwfin(`Wfin')
                local scpcy1 scpcy1path(`scpc1')
                local scpcy2 scpcy2path(`scpc2')
            }
            else {
                local scpcw
                local scpcy1
                local scpcy2
            }
            local legacy_savedyad `savedyad'
            if `automatic' local legacy_savedyad
            local cachemode fixed_threshold
            if "`thresholdoff'"!="" local cachemode thresholdoff
            load_data, nloc(`N') `clp' `distp' `legacy_savedyad' `savepath' `scpcw' `scpcy1' `scpcy2' `stringid' cachemode(`cachemode')
        }
        else {
            if "`filesuffix'"!="" {
                local savepath filesuffix(`filesuffix')
            }
            else {
                local savepath
            }

            use "`load'_est.dta", clear
            tmo_cache_normalization
            capture assert diag_schema==1 & diag_method=="legacy_override"
            if _rc {
                di as error "tmo: incompatible dyad cache; regenerate using the same legacy override method."
                exit 459
            }
            scalar __tmo_s_D = N_outcomes[1]
            scalar __tmo_nnovar = .
            capture scalar __tmo_nnovar = N_novar[1]

            use "`load'_dyad.dta", clear
            tmo_cache_normalization
            capture assert diag_schema==1 & diag_method=="legacy_override"
            if _rc {
                di as error "tmo: incompatible dyad cache schema."
                exit 459
            }
            if "`thresholdoff'"=="" {
                local cached_screen : char _dta[tmo_screen]
                if "`cached_screen'"!="selection" {
                    di as error "tmo: this dyad cache has no selection screen; regenerate it with threshold()."
                    exit 459
                }
            }
            tmo_dyad_cache_ids, nloc(`N') `stringid'
        }

        mkkeep, nloc(`N') dist_cutoff(`distthreshold')

        * Back out the finite-sample adjustment from the original SE. mkkeep
        * summed the original specification's own pairs.
        if scalar(__tmo_s_v_orig)<=0 {
            di as error "the original specification's variance estimate is not positive"
            exit 459
        }
        scalar __tmo_s_dof_adj = scalar(__tmo_s_se)^2/scalar(__tmo_s_v_orig)

        scalar __tmo_s_df_iqr = .
        scalar __tmo_s_sd_iqr = .
        local selection_mode fixed_threshold
        if "`thresholdoff'"!="" local selection_mode thresholdoff

        scalar __tmo_s_sel_eta   = .
        scalar __tmo_s_sel_gmax  = .
        scalar __tmo_s_sel_v0    = .
        scalar __tmo_s_sel_pilot = .
        scalar __tmo_s_sel_logf  = .
        local selkernel "none"
        local selpilot  "none"

        foreach m in __tmo_vpath_r __tmo_xipath_r __tmo_dgrid_r __tmo_etagrid_r {
            capture matrix drop `m'
        }
        est_tmo_se, custom_thres(`threshold') `thresholdoff' `scpc_uncond'

        if "`plotq'`plotse'`plothist'"!="" {
            di as txt "tmo: automatic-selection plots are skipped for legacy threshold overrides."
        }
        if `scpc_requested' {
            scalar __tmo_s_tmo_scpc_se = scalar(__tmo_s_tmo_se)
            tempvar legacy_Wfin
            rename Wfin `legacy_Wfin'
            est_tmo_se, custom_thres(`threshold') `thresholdoff'
            di as txt "tmo: legacy TMO and TMO-SCPC hybrid are returned separately."
        }




        clear
        return clear
        ereturn clear

        tempname tmo_results tmo_details
        mat `tmo_results' = J(1,6,.)
        mat rownames `tmo_results' = tmo
        mat colnames `tmo_results' = "Coef" "TMO SE" "t" "P>|t|" "95% Conf" "Interval"
        if missing(scalar(__tmo_s_df_r)) mat colnames `tmo_results' = "Coef" "TMO SE" "z" "P>|z|" "95% Conf" "Interval"
        mat `tmo_results'[1,1] = scalar(__tmo_s_beta)
        mat `tmo_results'[1,2] = scalar(__tmo_s_tmo_se)
        mat `tmo_results'[1,3] = scalar(__tmo_s_beta)/scalar(__tmo_s_tmo_se)
        tempname critical_value
        if missing(scalar(__tmo_s_df_r)) {
            scalar `critical_value' = invnormal(.975)
            mat `tmo_results'[1,4] = 2*normal(-abs(scalar(__tmo_s_beta)/scalar(__tmo_s_tmo_se)))
        }
        else {
            scalar `critical_value' = invttail(scalar(__tmo_s_df_r),.025)
            mat `tmo_results'[1,4] = 2*ttail(scalar(__tmo_s_df_r),abs(scalar(__tmo_s_beta)/scalar(__tmo_s_tmo_se)))
        }
        mat `tmo_results'[1,5] = scalar(__tmo_s_beta) - scalar(`critical_value')*scalar(__tmo_s_tmo_se)
        mat `tmo_results'[1,6] = scalar(__tmo_s_beta) + scalar(`critical_value')*scalar(__tmo_s_tmo_se)
        scalar __tmo_s_lb = `tmo_results'[1,5]
        scalar __tmo_s_ub = `tmo_results'[1,6]

        di as txt _n "TMO standard errors for `x'"
        matlist `tmo_results', border(all) cspec(o2& %20s | %9.3f o2 & %9.3f o2 & %6.2f o2 & %4.3f o2 & %9.3f o2 & %9.3f o2 &) rspec(&-&)

        mat `tmo_details' = J(8,1,.)
        mat rowname `tmo_details' = "Threshold" "eta at the threshold" "Largest candidate threshold" "% of off-diag in SE est." "% >= threshold (excl. clusters/Conley)" "# outcomes" "Degrees of freedom" "Pilot variance"
        mat colnames `tmo_details' = tmo
        mat `tmo_details'[1,1] = scalar(__tmo_s_thres)
        mat `tmo_details'[2,1] = scalar(__tmo_s_sel_eta)
        mat `tmo_details'[3,1] = scalar(__tmo_s_sel_gmax)
        mat `tmo_details'[4,1] = scalar(__tmo_s_offdP)*100
        mat `tmo_details'[5,1] = scalar(__tmo_s_offdPnocl)*100
        * cond() evaluates both arms, and __tmo_nout exists only after a
        * selection run, so branch explicitly.
        if scalar(__tmo_s_sel_v0)<. mat `tmo_details'[6,1] = scalar(__tmo_nout)
        else        mat `tmo_details'[6,1] = scalar(__tmo_s_D) - cond(scalar(__tmo_nnovar)<., scalar(__tmo_nnovar), 0)
        if scalar(__tmo_s_sel_v0)<. mat `tmo_details'[7,1] = 1/scalar(__tmo_s_sel_v0)
        else        mat `tmo_details'[7,1] = scalar(__tmo_s_df_iqr)
        mat `tmo_details'[8,1] = scalar(__tmo_s_sel_pilot)
        matlist `tmo_details', cspec(& %38s | %9.4f &) rspec(& & & & & & & - & &) coleqonly
        if "`selkernel'"!="none" {
            di as txt "  threshold selected on a `=colsof(__tmo_etagrid_r)'-point eta grid" ///
                ", pilot from the `selpilot' variance, `selkernel' kernel"
        }

        return scalar beta = scalar(__tmo_s_beta)
        return scalar orig_se = scalar(__tmo_s_se)
        return scalar tmo_se = scalar(__tmo_s_tmo_se)
        return scalar variance = scalar(__tmo_s_tmo_variance)
        return scalar lb = scalar(__tmo_s_lb)
        return scalar ub = scalar(__tmo_s_ub)
        return scalar threshold = scalar(__tmo_s_thres)
        return scalar pct_ge_thres = scalar(__tmo_s_offdP)*100
        return scalar pct_ge_thres_nocl = scalar(__tmo_s_offdPnocl)*100

        return scalar eta          = scalar(__tmo_s_sel_eta)
        return scalar grid_max     = scalar(__tmo_s_sel_gmax)
        return scalar v0hat        = scalar(__tmo_s_sel_v0)
        return scalar pilot        = scalar(__tmo_s_sel_pilot)
        return scalar logfactor    = scalar(__tmo_s_sel_logf)
        if "`selkernel'"!="none" {
            return scalar N_used     = scalar(__tmo_nout)
            return scalar N_dropped  = scalar(__tmo_ndropped)
            return scalar incomplete = scalar(__tmo_hasmiss)
            return scalar max_eligible_pct = ///
                scalar(__tmo_s_sel_max_eligible_pct)
            return scalar admissible = 1
            return matrix vpath   = __tmo_vpath_r
            return matrix xipath  = __tmo_xipath_r
            return matrix dgrid   = __tmo_dgrid_r
            return matrix etagrid = __tmo_etagrid_r
            if scalar(__tmo_need_counts) {
                return matrix screened_countpath = __tmo_screened_counts_r
                return matrix eligiblepctpath = __tmo_eligible_pctpath_r
                return matrix admissiblepath = __tmo_admissiblepath_r
            }
        }
        return local  pilot_source "`selpilot'"
        return local  kernel       "`selkernel'"

        return scalar T = scalar(__tmo_s_T)
        return scalar N_loc = scalar(__tmo_s_N)
        return scalar N_clust = scalar(__tmo_s_N_clust)
        return scalar N_outcomes = scalar(__tmo_s_D)
        return scalar N_novar = scalar(__tmo_nnovar)
        return scalar N_depend = scalar(__tmo_ndepend)
        * auxiliary outcomes entering the screen; missing under thresholdoff
        cap return scalar N_cols_used = scalar(__tmo_s_D_cols)
        return scalar N = scalar(__tmo_s_N_obs)
        return scalar dof = scalar(__tmo_s_df_iqr)
        return scalar finite_sample_dof = scalar(__tmo_s_dof_adj)
        return scalar df_r = scalar(__tmo_s_df_r)
        return scalar scpc_cv = ${scpc_cv}
        return scalar diag_schema = 1
        return local diag_method "legacy_override"
        return local selection_mode "`selection_mode'"
        return scalar scpc_se = scalar(__tmo_s_scpc_se)
        return scalar tmo_scpc_se = scalar(__tmo_s_tmo_scpc_se)
    restore

    if "`saveest'"!="" {
        tmo_save, `savepath' mode(`selection_mode')
    }

    tmo_cleanup
end

* Inside the guard: this program is on the call stack when tmo_work re-runs
* this file as __tmo_helpers_only, so it must not be redefined there.
capture program drop tmo_session_work
program define tmo_session_work, rclass
    version 13
    local invocation : copy local 0
    tempname incoming_results fit_results
    _return hold `incoming_results'

    * Only the SCPC path can overwrite the external backend's matrix names.
    capture syntax, [SCPC scpc_cmd(string asis) *]
    local parse_rc = _rc
    local scpc_requested 0
    if !`parse_rc' local scpc_requested = ("`scpc'"!="" | `"`scpc_cmd'"'!="")

    local owned_matrices
    if `scpc_requested' local owned_matrices ///
        b coef scpc_Bread scpc_Bread0 scpc_cvs scpc_cvs_c scpcstats scpctab

    local caller_globals : all globals
    local cv_name scpc_cv
    local had_cv : list cv_name in caller_globals
    if `had_cv' local caller_cv : copy global scpc_cv

    * Stata tempnames survive an external package's Mata clear.
    local scalar_names : all scalars
    local saved_scalars
    local scalar_copies
    foreach object of local scalar_names {
        if substr("`object'",1,6)=="__tmo_" {
            tempname copy
            capture scalar `copy' = scalar(`object')
            local copy_rc = _rc
            if `copy_rc' {
                _return restore `incoming_results'
                return add
                exit `copy_rc'
            }
            local saved_scalars `saved_scalars' `object'
            local scalar_copies `scalar_copies' `copy'
        }
    }
    local matrix_names : all matrices
    local saved_matrices
    local matrix_copies
    foreach object of local matrix_names {
        local exact_owned : list object in owned_matrices
        if `exact_owned' | substr("`object'",1,6)=="__tmo_" {
            tempname copy
            capture matrix `copy' = `object'
            local copy_rc = _rc
            if `copy_rc' {
                _return restore `incoming_results'
                return add
                exit `copy_rc'
            }
            local saved_matrices `saved_matrices' `object'
            local matrix_copies `matrix_copies' `copy'
        }
    }

    _return restore `incoming_results', hold
    scalar __tmo_scpc_owned = 0
    capture noisily tmo_work `invocation'
    local fit_rc = _rc
    if !`fit_rc' _return hold `fit_results'
    else capture tmo_cleanup

    * Enumerate exact names: never use an abbreviated scalar reference to drop.
    local scalar_names : all scalars
    foreach object of local scalar_names {
        if substr("`object'",1,6)=="__tmo_" {
            scalar drop `object'
        }
    }
    local matrix_names : all matrices
    foreach object of local matrix_names {
        local exact_owned : list object in owned_matrices
        if `exact_owned' | substr("`object'",1,6)=="__tmo_" {
            matrix drop `object'
        }
    }
    foreach object of local saved_scalars {
        gettoken copy scalar_copies : scalar_copies
        scalar `object' = scalar(`copy')
    }
    foreach object of local saved_matrices {
        gettoken copy matrix_copies : matrix_copies
        matrix `object' = `copy'
    }
    if `had_cv' global scpc_cv : copy local caller_cv
    else capture macro drop scpc_cv

    if `fit_rc' {
        _return restore `incoming_results'
        return add
        exit `fit_rc'
    }
    _return restore `fit_results'
    return add
end

}

capture program drop tmo_outcome_controls
program define tmo_outcome_controls
    version 13
    syntax, OUTcomes(varlist)
    tempname design_b fe_dof
    matrix `design_b' = e(b)
    local terms : colnames `design_b'
    local partial `e(inexog0)' `e(partial1)' `e(partial)'
    local absorbed
    capture matrix `fe_dof' = e(dof_table)
    if !_rc local absorbed : rownames `fe_dof'
    else local absorbed `e(extended_absvars)'
    if "`absorbed'"=="" local absorbed `e(absvars)'
    local terms `terms' `partial' `absorbed'
    local terms : list uniq terms
    local used
    foreach term of local terms {
        if "`term'"!="_cons" {
            capture _ms_parse_parts `term'
            if !_rc {
                local kind `r(type)'
                if inlist("`kind'","interaction","product") {
                    local count = r(k_names)
                    forvalues j=1/`count' {
                        local used `used' `r(name`j')'
                    }
                }
                else local used `used' `r(name)'
            }
            else {
                * Compact native FV lists can lack coefficient-stripe levels.
                quietly fvrevar `term', list
                local used `used' `r(varlist)'
            }
        }
    }
    local used : list uniq used
    local overlap : list outcomes & used
    if "`overlap'"!="" {
        di as error "tmo: outcome variable(s) also appear in the regression controls: `overlap'"
        di as error "Remove these variables from the controls or from ylist()."
        exit 198
    }
end

capture program drop tmo_finite_sample_adjustment
program define tmo_finite_sample_adjustment, rclass
    version 13
    syntax, SPEC(string)
    tempname b omitted k n dfm dfa nested rankxx dofminus sdofminus
    scalar `n' = e(N)
    scalar `k' = .

    if "`spec'"=="regress" {
        matrix `b' = e(b)
        _ms_omit_info `b'
        matrix `omitted' = r(omit)
        mata: st_numscalar("`k'", cols(st_matrix("`b'"))-sum(st_matrix("`omitted'")))
    }
    else if "`spec'"=="reghdfe" {
        scalar `dfm' = e(df_m)
        scalar `dfa' = e(df_a)
        capture scalar `nested' = e(df_a_nested)
        if _rc | missing(`nested') scalar `nested' = 0
        scalar `k' = `dfm' + `dfa' + (`nested'>0)
    }
    else {
        capture scalar `rankxx' = e(rankxx)
        if _rc | missing(`rankxx') scalar `rankxx' = e(rank)
        capture scalar `dofminus' = e(dofminus)
        if _rc | missing(`dofminus') scalar `dofminus' = 0
        capture scalar `sdofminus' = e(sdofminus)
        if _rc | missing(`sdofminus') scalar `sdofminus' = 0
        scalar `k' = `rankxx' + `dofminus' + `sdofminus'
    }

    if missing(`n') | missing(`k') | `n'<=`k' | `k'<=0 {
        di as error "tmo: cannot form the finite-sample factor n/(n-k)."
        di as error "     The fitted command reports n=" %12.0g `n' " and k=" %12.0g `k' "."
        exit 459
    }
    return scalar k = `k'
    return scalar factor = `n'/(`n'-`k')
end


if "`1'"!="__tmo_scpc_backend_helpers_only" {
capture program drop tmo_scpc_backend
program define tmo_scpc_backend
    version 13
    capture quietly program list scpc
    local needload = _rc
    capture mata: assert(tmo_scpc_backend_ready()==1)
    if _rc local needload 1
    if !`needload' {
        capture mata: assert(getavc(1,J(2,1,1))>0)
        local needload = _rc
    }
    if `needload' {
        quietly findfile scpc.ado
        local backend `"`r(fn)'"'
        quietly checksum `"`backend'"'
        * Compatible SCPC sources differ only in their critical-value export comment.
        local compatible = (r(checksum)==2275767291 & r(filelen)==27017) | ///
            (r(checksum)==989185685 & r(filelen)==27016)
        if !`compatible' {
            di as error "tmo: this SCPC source is not supported by the runtime adapter."
            di as error "Install the compatible SCPC companion listed in help tmo (Dependencies)."
            exit 198
        }
        quietly run `"`backend'"'
        quietly findfile tmo.ado
        quietly run `"`r(fn)'"' __tmo_scpc_backend_helpers_only
    }
    capture mata: tmo_lowtri_gram(J(0,0,.))
    if _rc {
        quietly findfile tmo.ado
        quietly run `"`r(fn)'"' __tmo_helpers_only
    }
end
}

if "`1'"=="__tmo_scpc_backend_helpers_only" {
capture mata: mata drop tmo_scpc_backend_ready()
capture mata: mata drop set_Wx_cluster()
mata:
real scalar tmo_scpc_backend_ready()
{
    external real scalar __tmo_scpc_backend_patch
    return(__tmo_scpc_backend_patch==1)
}

void set_Wx_cluster(string scalar Vname, string scalar sel)
{
    external real scalar condflag
    external real matrix Wfin, Wx
    external real vector permfin
    real vector r
    real scalar i
    string scalar cluster

    cluster=st_global("e(clustvar)")
    if(cluster=="" | cols(tokens(cluster))!=1) {
        errprintf("tmo: conditional SCPC requires the current single cluster variable.\n")
        _error(459)
    }
    condflag=1
    r=J(rows(Wfin),1,0)
    Wx=Wfin
    for(i=1;i<=cols(Wfin);i++) {
        r[permfin]=Wfin[.,i]
        stata("capture drop scpc_r")
        stata("quietly generate scpc_r=.")
        stata("quietly replace scpc_r=0 if e(sample)")
        st_store(.,"scpc_r",sel,r)
        stata("capture drop scpc_rx")
        stata("quietly gegen scpc_rx=total(scpc_r) if e(sample), by("+cluster+")")
        stata("capture drop scpc_r")
        if(i==1) {
            stata("quietly generate scpc_r=scpc_rx*scpc_xs if e(sample)")
        }
        else {
            stata("quietly replace scpc_rx=scpc_rx*scpc_xs")
            stata("_estimates hold oreg, restore")
            stata(Vname)
            stata("quietly predict scpc_r, residuals")
            stata("_estimates unhold oreg")
        }
        stata("capture drop scpc_rx")
        stata("quietly hashsort "+cluster)
        stata("quietly gegen scpc_rx=total(scpc_r*scpc_x) if e(sample), by("+cluster+")")
        Wx[.,i]=st_data(.,"scpc_rx",sel)[permfin]
    }
    set_ccspc()
}
__tmo_scpc_backend_patch=1
end
}


if "`1'"!="__tmo_scpc_auto_helpers_only" {
capture program drop tmo_scpc_auto
program define tmo_scpc_auto, rclass
    version 13
    syntax, CMD(string) X(varname) SAMPLE(varname) IDVAR(varname) ///
        LATitude(varname) LONgitude(varname) [UNCOND]
    capture mata: tmo_scpc_auto_ready()
    if _rc {
        quietly findfile tmo.ado
        quietly run `"`r(fn)'"' __tmo_scpc_auto_helpers_only
    }

    local source `e(cmd)'
    local y `e(depvar)'
    local cluster `e(clustvar)'
    local absorb `e(absvars)'
    local isiv = inlist("`source'","ivreg2","ivreghdfe")
    if !inlist("`source'","regress","reghdfe","ivreg2","ivreghdfe") {
        di as error "tmo: scpc does not support estimator `source'."
        exit 198
    }
    if "`e(wtype)'`e(wexp)'"!="" {
        di as error "tmo: weights are not supported by scpc."
        exit 198
    }
    if `: word count `cluster''>1 {
        di as error "tmo: scpc requires a single spatial clustering variable."
        exit 198
    }
    mata: st_local("cmd_comma",strofreal(tmo_optcomma(st_local("cmd"))))
    if inlist("`source'","reghdfe","ivreghdfe") & `cmd_comma'>0 {
        local option_tail = substr(`"`cmd'"',`cmd_comma'+1,.)
        while strtrim(`"`option_tail'"')!="" {
            gettoken option_part option_tail : option_tail, bind
            local option_open = strpos(`"`option_part'"',"(")
            local option_name = lower(strtrim(`"`option_part'"'))
            if `option_open'>0 local option_name = lower(strtrim(substr(`"`option_part'"',1,`option_open'-1)))
            local team_option = (strlen("`option_name'")>=1 & ///
                "`option_name'"==substr("group_id",1,strlen("`option_name'"))) | ///
                (strlen("`option_name'")>=1 & ///
                "`option_name'"==substr("individual_id",1,strlen("`option_name'")))
            if `team_option' {
                di as error "tmo: scpc does not support group()/individual() fixed effects."
                exit 198
            }
        }
    }
    if `isiv' {
        if "`e(model)'"!="iv" | "`e(partial)'`e(partial1)'"!="" {
            di as error "tmo: scpc supports ordinary 2SLS without partial(), not LIML or GMM."
            exit 198
        }
        if "`e(instd1)'"!="`x'" {
            di as error "tmo: scpc requires x() to be the sole retained endogenous variable."
            exit 198
        }
    }
    local vce `e(vce)'
    if !inlist("`vce'","ols","unadjusted","robust","cluster","robust cluster","") {
        di as error "tmo: scpc does not support this covariance estimator (`vce')."
        exit 198
    }
    if "`e(kernel)'"!="" | (!missing(e(bw)) & e(bw)>0) {
        di as error "tmo: scpc does not translate HAC covariance specifications."
        exit 198
    }

    tempname native_b native_V native_omit native_beta native_var native_rss ///
        native_n native_df native_rank native_df_a native_df_m native_small ///
        native_dofminus native_sdofminus native_nested native_g
    matrix `native_b' = e(b)
    matrix `native_V' = e(V)
    scalar `native_beta' = _b[`x']
    scalar `native_var' = _se[`x']^2
    scalar `native_rss' = e(rss)
    scalar `native_n' = e(N)
    scalar `native_df' = e(df_r)
    scalar `native_rank' = e(rank)
    if `isiv' scalar `native_rank' = e(rankxx)
    scalar `native_df_a' = e(df_a)
    scalar `native_df_m' = e(df_m)
    scalar `native_dofminus' = e(dofminus)
    scalar `native_sdofminus' = e(sdofminus)
    scalar `native_nested' = e(df_a_nested)
    scalar `native_g' = e(N_clust)
    local native_small_flag `e(small)'
    scalar `native_small' = ("`native_small_flag'"!="")
    local exinst `e(exexog1)'
    _ms_omit_info `native_b'
    matrix `native_omit' = r(omit)
    local native_names : colnames `native_b'
    local native_cons 0
    local retained
    local native_j 0
    foreach term of local native_names {
        local ++native_j
        if `native_omit'[1,`native_j']==0 {
            if "`term'"=="_cons" local native_cons 1
            else local retained `retained' `term'
        }
    }
    if !`: list x in retained' {
        di as error "tmo: scpc requires a retained main-effect coefficient in x()."
        exit 198
    }
    if missing(`native_var') | `native_var'<=0 {
        di as error "tmo: scpc requires a positive finite conventional target variance."
        exit 459
    }
    quietly count if `sample'!=e(sample)
    if r(N) {
        di as error "tmo: scpc source sample does not match the main estimation sample."
        exit 459
    }
    quietly count if `sample'
    if r(N)!=`native_n' {
        di as error "tmo: scpc source observation count does not match e(N)."
        exit 459
    }
    quietly count if `sample' & (missing(`latitude',`longitude',`idvar') | ///
        abs(`latitude')>90 | abs(`longitude')>180)
    if r(N) {
        di as error "tmo: scpc requires valid coordinates and location IDs on every main-sample row."
        exit 459
    }
    if "`cluster'"!="" {
        quietly count if `sample' & missing(`cluster')
        if r(N) {
            di as error "tmo: scpc requires nonmissing main-sample cluster IDs."
            exit 459
        }
    }

    * Obtain structural residuals from the original estimator, not fitted treatment.
    tempvar native_resid native_es
    quietly generate byte `native_es' = e(sample)
    capture quietly predict double `native_resid' if `sample', residuals
    if _rc {
        if !inlist("`source'","reghdfe","ivreghdfe") {
            di as error "tmo: scpc could not obtain original structural residuals."
            exit 459
        }
        local delim = cond(`cmd_comma'>0,"",",")
        quietly `cmd' `delim' residuals(`native_resid')
        quietly count if `native_es'!=e(sample)
        if r(N) | reldif(_b[`x'],`native_beta')>1e-9 {
            di as error "tmo: obtaining HDFE residuals changed the original sample or coefficient."
            exit 459
        }
        tempname refit_V
        matrix `refit_V' = e(V)
        mata: assert(mreldif(st_matrix("`native_V'"),st_matrix("`refit_V'"))<1e-9)
    }
    quietly count if `sample' & missing(`native_resid')
    if r(N) {
        di as error "tmo: original structural residuals are missing in the main sample."
        exit 459
    }

    * Materialize native retained columns; the external score parser needs variables.
    tempvar target
    quietly generate double `target' = `x' if `sample'
    local controls
    foreach term of local retained {
        if "`term'"!="`x'" {
            local fvterm
            local pieces : subinstr local term "#" " ", all
            foreach piece of local pieces {
                if regexm("`piece'","^([0-9]+)\.(.+)$") {
                    local piece `=regexs(1)'bn.`=regexs(2)'
                }
                if "`fvterm'"=="" local fvterm `piece'
                else local fvterm `fvterm'#`piece'
            }
            quietly fvrevar `fvterm'
            local expanded `r(varlist)'
            if `: word count `expanded''!=1 {
                di as error "tmo: scpc cannot unambiguously materialize retained term `term'."
                exit 198
            }
            tempvar col
            quietly generate double `col' = `expanded' if `sample'
            local controls `controls' `col'
        }
    }
    local instvars
    foreach term of local exinst {
        local fvterm
        local pieces : subinstr local term "#" " ", all
        foreach piece of local pieces {
            if regexm("`piece'","^([0-9]+)\.(.+)$") {
                local piece `=regexs(1)'bn.`=regexs(2)'
            }
            if "`fvterm'"=="" local fvterm `piece'
            else local fvterm `fvterm'#`piece'
        }
        quietly fvrevar `fvterm'
        foreach z of varlist `r(varlist)' {
            tempvar col
            quietly generate double `col' = `z' if `sample'
            local instvars `instvars' `col'
        }
    }
    local translated_cons `native_cons'
    foreach fe of local absorb {
        if "`fe'"=="_cons" {
            local translated_cons 1
            continue
        }
        if strpos("`fe'",".") | strpos("`fe'","=") | strpos("`fe'","(") | strpos("`fe'",")") {
            di as error "tmo: scpc currently supports categorical absorbed effects only."
            exit 198
        }
        local fevars = subinstr("`fe'","#"," ",.)
        capture confirm variable `fevars'
        if _rc {
            di as error "tmo: scpc cannot resolve absorbed effect `fe'."
            exit 198
        }
        tempvar fegroup
        quietly gegen long `fegroup' = group(`fevars') if `sample', missing
        quietly gstats sum `fegroup' if `sample', meanonly
        local levels = r(max)
        if `levels'+c(k)>c(maxvar)-32 {
            di as error "tmo: scpc fixed-effect expansion exceeds the available variable limit."
            exit 198
        }
        if `levels'>1 {
            forvalues k=2/`levels' {
                tempvar col
                quietly generate byte `col' = (`fegroup'==`k') if `sample'
                local controls `controls' `col'
            }
        }
        local translated_cons 1
    }
    if !`translated_cons' & "`controls'"!="" {
        mata: st_local("translated_cons",strofreal(tmo_scpc_auto_constant("`controls'","`sample'")))
    }
    if !`translated_cons' & ("`uncond'"=="" | `isiv') {
        di as error "tmo: conditional SCPC and IV SCPC require a constant in the model space."
        exit 198
    }
    local consopt
    if !`translated_cons' local consopt noconstant

    tempvar ycopy idcopy latcopy loncopy group native_score
    quietly generate double `ycopy' = `y'
    quietly generate double `idcopy' = `idvar'
    quietly generate double `latcopy' = `latitude'
    quietly generate double `loncopy' = `longitude'
    if "`cluster'"!="" quietly gegen long `group' = group(`cluster') if `sample'
    else quietly generate long `group' = _n if `sample'
    quietly keep if `sample'
    tempvar coordtag
    quietly gegen byte `coordtag' = tag(`latcopy' `loncopy')
    quietly count if `coordtag'
    if r(N)<2 {
        di as error "tmo: scpc requires at least two distinct coordinate pairs."
        exit 459
    }
    quietly hashsort `group'
    capture quietly by `group': assert `latcopy'==`latcopy'[1] & `loncopy'==`loncopy'[1]
    if _rc {
        di as error "tmo: scpc requires constant coordinates within each cluster."
        exit 459
    }
    quietly hashsort `idcopy' `group'
    capture quietly by `idcopy': assert `group'==`group'[1]
    if _rc {
        di as error "tmo: scpc requires repeated locations to belong to a single cluster."
        exit 459
    }
    quietly drop `native_es'
    quietly gegen byte `native_es' = tag(`group')
    quietly count if `native_es'
    local ngroups = r(N)
    if `ngroups'<5 {
        di as error "tmo: scpc requires at least five distinct spatial units."
        exit 198
    }
    if "`cluster'"!="" & `ngroups'!=`native_g' {
        di as error "tmo: scpc cluster count differs from the main estimator."
        exit 459
    }
    quietly hashsort `group' `idcopy'
    quietly generate double `native_score' = .
    tempname diag
    mata: tmo_scpc_auto_native("`target'","`controls'","`instvars'", ///
        "`native_resid'","`group'","`native_score'",`translated_cons',`isiv', ///
        "`source'","`vce'","`cluster'",st_numscalar("`native_beta'"), ///
        st_numscalar("`native_var'"),st_numscalar("`native_rss'"), ///
        st_numscalar("`native_df'"),st_numscalar("`native_rank'"), ///
        st_numscalar("`native_df_a'"),st_numscalar("`native_df_m'"), ///
        st_numscalar("`native_small'"),st_numscalar("`native_dofminus'"), ///
        st_numscalar("`native_sdofminus'"),st_numscalar("`native_nested'"),"`diag'")

    local vcetrans vce(robust)
    if "`cluster'"!="" local vcetrans vce(cluster `group')
    if `isiv' local adapter ivregress 2sls `ycopy' `controls' (`target' = `instvars'), `vcetrans' `consopt'
    else local adapter regress `ycopy' `target' `controls', `vcetrans' `consopt'
    quietly `adapter'
    quietly count if !e(sample)
    if r(N) | reldif(_b[`target'],`native_beta')>1e-7 {
        di as error "tmo: translated SCPC sample or target coefficient differs from the main fit."
        exit 459
    }
    tempvar translated_resid
    quietly predict double `translated_resid', residuals
    mata: tmo_scpc_auto_translated("`target'","`controls'","`instvars'", ///
        "`native_resid'","`translated_resid'","`native_score'",`translated_cons',`isiv',"`diag'")

    * Drop collinear translated columns before SCPC's literal score parser runs.
    tempname tb omit
    matrix `tb' = e(b)
    _ms_omit_info `tb'
    matrix `omit' = r(omit)
    local tnames : colnames `tb'
    local cleancontrols
    local j 0
    foreach term of local tnames {
        local ++j
        if `omit'[1,`j']==0 & !inlist("`term'","`target'","_cons") local cleancontrols `cleancontrols' `term'
    }
    if "`cleancontrols'"!="`controls'" {
        if `isiv' local adapter ivregress 2sls `ycopy' `cleancontrols' (`target' = `instvars'), `vcetrans' `consopt'
        else local adapter regress `ycopy' `target' `cleancontrols', `vcetrans' `consopt'
        quietly `adapter'
        if reldif(_b[`target'],`native_beta')>1e-7 {
            di as error "tmo: removing omitted SCPC columns changed the target coefficient."
            exit 459
        }
        quietly drop `translated_resid'
        quietly predict double `translated_resid', residuals
        mata: tmo_scpc_auto_translated("`target'","`cleancontrols'","`instvars'", ///
            "`native_resid'","`translated_resid'","`native_score'",`translated_cons',`isiv',"`diag'")
    }
    local target_index = colnumb(e(b),"`target'")
    if `target_index'!=1 {
        di as error "tmo: translated SCPC target is not the first fitted coefficient."
        exit 459
    }
    * Work in an isolated variable namespace; SCPC owns s_* and scpc_* names.
    local keepvars `ycopy' `target' `cleancontrols' `instvars' `group' `idcopy' `latcopy' `loncopy'
    local keepvars : list uniq keepvars
    quietly keep `keepvars'
    capture confirm variable s_1
    if !_rc {
        di as error "tmo: scpc cannot use s_1 as a model variable."
        exit 198
    }
    capture confirm variable s_2
    if !_rc {
        di as error "tmo: scpc cannot use s_2 as a model variable."
        exit 198
    }
    rename `latcopy' s_1
    rename `loncopy' s_2
    order s_1 s_2
    scalar __tmo_scpc_owned = 1
    tmo_scpc_backend
    quietly scpc, k(1) latlong `uncond'
    tempname result_se
    scalar `result_se' = e(scpcstats)[1,2]
    if missing(`result_se') | `result_se'<=0 {
        di as error "tmo: scpc returned a nonpositive or missing standard error."
        exit 459
    }
    capture mata: tmo_lowtri_gram(J(0,0,.))
    if _rc {
        quietly findfile tmo.ado
        quietly run `"`r(fn)'"' __tmo_helpers_only
    }
    if "`cluster'"=="" {
        quietly hashsort `group'
        quietly replace `group' = _n
    }
    mata: tmo_scpc_capture("`idcopy'","`group'")
    tmo_scpc_release
    return scalar scpc_se = `result_se'
    return scalar sample_n = `native_n'
    return scalar spatial_n = `ngroups'
    return scalar finite_correction = `diag'[1,1]
    return scalar variance_relative_error = `diag'[1,2]
    return scalar residual_relative_error = `diag'[1,3]
    return scalar influence_relative_error = `diag'[1,4]
    return local cluster `cluster'
    return local source_cmd `source'
    return local adapter_cmd `"`adapter'"'
end
}

capture mata: mata drop tmo_scpc_auto_ready()
capture mata: mata drop tmo_scpc_auto_constant()
capture mata: mata drop tmo_scpc_auto_native()
capture mata: mata drop tmo_scpc_auto_translated()
mata:
real scalar tmo_scpc_auto_ready()
{
    return(1)
}

real scalar tmo_scpc_auto_constant(string scalar cvars, string scalar samplevar)
{
    real matrix C
    real colvector one,r
    C=st_data(.,tokens(cvars),samplevar)
    if(cols(C)==0 | rows(C)==0 | hasmissing(C)) return(0)
    one=J(rows(C),1,1)
    r=one-C*(pinv(C)*one)
    return(max(abs(r))<1e-10)
}

void tmo_scpc_auto_native(string scalar xvar, string scalar cvars, string scalar zvars,
    string scalar uvar, string scalar gvar, string scalar scorevar,
    real scalar cons, real scalar isiv, string scalar source, string scalar vce,
    string scalar cluster, real scalar beta, real scalar variance, real scalar rss,
    real scalar dfr, real scalar native_rank, real scalar dfa, real scalar dfm,
    real scalar small, real scalar dofminus, real scalar sdofminus, real scalar nested,
    string scalar out)
{
    real matrix C, X, Z, H, bread, info
    real colvector x,u,q,s,g,total
    real scalar n,G,f,base,err
    x=st_data(.,xvar); u=st_data(.,uvar); g=st_data(.,gvar); n=rows(x)
    C=cvars=="" ? J(n,0,.) : st_data(.,tokens(cvars))
    if(cons) C=(C,J(n,1,1))
    X=(x,C)
    if(isiv) {
        Z=(st_data(.,tokens(zvars)),C)
        H=Z*(pinv(Z)*X)
    }
    else H=X
    bread=pinv(cross(X,H)); q=H*bread[,1]
    if(hasmissing(q) | max(abs(q))==0 | abs(cross(q,x)-1)>1e-7) {
        errprintf("tmo: SCPC translated design cannot identify the original target.\n"); _error(459)
    }
    if(cols(C)>0 & max(abs(cross(C,q)))>1e-7*max((1,max(abs(C))))) {
        errprintf("tmo: SCPC target influence fails the native control-space assertion.\n"); _error(459)
    }
    if(abs(cross(u,u)-rss)>1e-7*max((1,rss))) {
        errprintf("tmo: SCPC original structural residuals do not reproduce native RSS.\n"); _error(459)
    }
    s=q:*u; info=panelsetup(g,1); G=rows(info)
    if(missing(dofminus)) dofminus=0
    if(missing(sdofminus)) sdofminus=0
    if(cluster!="") {
        total=panelsum(s,info); base=cross(total,total)
        if(isiv) f=small ? ((n-1)/(n-native_rank-sdofminus))*G/(G-1) : 1
        else if(source=="reghdfe") {
            if(missing(nested)) nested=0
            f=(n-1)/(n-dfm-dfa-(nested>0))*G/(G-1)
        }
        else f=(n-1)/(n-native_rank)*G/(G-1)
    }
    else if(vce=="robust") {
        base=cross(s,s)
        if(isiv) f=small ? (n-dofminus)/(n-native_rank-dofminus-sdofminus) : 1
        else f=n/dfr
    }
    else {
        base=cross(q,q)*rss
        if(isiv) f=small ? 1/(n-native_rank-dofminus-sdofminus) : 1/(n-dofminus)
        else f=1/dfr
    }
    err=abs(f*base-variance)/max((1e-30,abs(variance)))
    if(missing(f) | f<=0 | err>2e-6) {
        errprintf("tmo: SCPC source variance/finite-sample assertion failed (relative error %g).\n",err)
        _error(459)
    }
    st_store(.,scorevar,s)
    st_matrix(out,(f,err,.,.))
}

void tmo_scpc_auto_translated(string scalar xvar, string scalar cvars, string scalar zvars,
    string scalar native_uvar, string scalar uvar, string scalar native_scorevar,
    real scalar cons, real scalar isiv, string scalar out)
{
    real matrix C,X,Z,H,V,D
    real colvector x,u,u0,q,s0
    real scalar n,re,se
    x=st_data(.,xvar); u=st_data(.,uvar); u0=st_data(.,native_uvar)
    s0=st_data(.,native_scorevar); n=rows(x)
    C=cvars=="" ? J(n,0,.) : st_data(.,tokens(cvars))
    if(cons) C=(C,J(n,1,1))
    X=(x,C)
    if(isiv) {
        Z=(st_data(.,tokens(zvars)),C)
        H=Z*(pinv(Z)*X)
    }
    else H=X
    V=st_matrix("e(V_modelbased)")
    if(cols(V)!=cols(H)) {
        errprintf("tmo: SCPC translated bread and retained design dimensions disagree.\n")
        _error(459)
    }
    q=H*V[,1]
    re=max(abs(u-u0))/max((1,max(abs(u0))))
    se=max(abs(q:*u-s0))/max((1e-30,max(abs(s0))))
    if(re>1e-7 | se>2e-6 | missing(re) | missing(se)) {
        errprintf("tmo: SCPC structural residual/influence assertion failed (%g, %g).\n",re,se)
        _error(459)
    }
    D=st_matrix(out); D[1,3]=re; D[1,4]=se; st_matrix(out,D)
}
end





capture program drop tmo_distance_literal
program define tmo_distance_literal, rclass
    syntax, [DISTTHREShold(string) *]
    return local distance `"`distthreshold'"'
end

capture program drop tmo_selection_capture
program define tmo_selection_capture
    scalar __tmo_s_sel_thres = r(threshold)
    scalar __tmo_s_sel_se = r(tmo_se)
    scalar __tmo_s_sel_var = r(tmo_var)
    scalar __tmo_s_sel_eta = r(eta)
    scalar __tmo_s_sel_index = r(selected_index)
    scalar __tmo_s_sel_gmax = r(grid_max)
    scalar __tmo_s_sel_v0 = r(v0hat)
    scalar __tmo_s_sel_pilot = r(pilot)
    scalar __tmo_s_sel_gamma0 = r(gamma0)
    scalar __tmo_s_sel_logf = r(logfactor)
    scalar __tmo_s_sel_q25 = r(q25)
    scalar __tmo_s_sel_q75 = r(q75)
    scalar __tmo_s_sel_max_eligible_pct = r(max_eligible_pct)
    scalar __tmo_s_sel_interior = r(interior)
    scalar __tmo_s_sel_nprotected = r(N_protected_pairs)
    scalar __tmo_s_sel_neligible = r(N_eligible_pairs)
    scalar __tmo_s_sel_calweight = r(calibration_weight)
    scalar __tmo_s_sel_unsupported = r(unsupported_pairs)
    scalar __tmo_s_sel_screened = .
    if scalar(__tmo_need_counts) {
        matrix __tmo_screened_counts_r = r(screened_counts)
        matrix __tmo_eligible_pctpath_r = r(eligiblepctpath)
        matrix __tmo_admissiblepath_r = r(admissiblepath)
        scalar __tmo_s_sel_screened = __tmo_screened_counts_r[1,scalar(__tmo_s_sel_index)]
    }
    matrix __tmo_vpath_r = r(vpath)
    matrix __tmo_xipath_r = r(xipath)
    matrix __tmo_dgrid_r = r(dgrid)
    matrix __tmo_etagrid_r = r(etagrid)
    matrix __tmo_anchorpath_r = r(anchorpath)
    matrix __tmo_anchorstats_r = r(anchorstats)
    matrix __tmo_objectivepath_r = r(objectivepath)
    scalar __tmo_s_sel_penalty = __tmo_xipath_r[1,scalar(__tmo_s_sel_index)]
    scalar __tmo_s_sel_objective = __tmo_objectivepath_r[1,scalar(__tmo_s_sel_index)]
    c_local selkernel `r(kernel)'
    c_local selpilot `r(pilot_source)'
end

capture program drop tmo_selection_output
program define tmo_selection_output, rclass
    syntax, Xname(string) KERNEL(string) PILOTSOURCE(string) LIBpath(string) RETENTION(string) ///
        [DIAGnostics SAVEdyad SAVEest PLOTq PLOTse PLOThist SAVEPLOTSEest ///
        FILEsuffix(string) NBINS(integer 10000) STRINGID]
    scalar __tmo_s_tmo_se = scalar(__tmo_s_sel_se)
    scalar __tmo_s_thres = scalar(__tmo_s_sel_thres)
    scalar __tmo_s_offdP = .
    scalar __tmo_s_offdPnocl = .
    scalar __tmo_retained = .
    scalar __tmo_screened = .
    scalar __tmo_unsupported = .
    scalar __tmo_histweight = scalar(__tmo_s_sel_calweight)
    scalar __tmo_diagnostic_preparations = 0
    tempname critical_value
    if missing(scalar(__tmo_s_df_r)) scalar `critical_value' = invnormal(.975)
    else scalar `critical_value' = invttail(scalar(__tmo_s_df_r),.025)
    scalar __tmo_s_lb = scalar(__tmo_s_beta)-scalar(`critical_value')*scalar(__tmo_s_tmo_se)
    scalar __tmo_s_ub = scalar(__tmo_s_beta)+scalar(`critical_value')*scalar(__tmo_s_tmo_se)
    local pairwork = ///
        ("`diagnostics'`savedyad'`saveest'`plothist'`saveplotseest'"!="" | ///
        scalar(__tmo_s_sel_max_eligible_pct)<100)
    local pairrecompute = ("`savedyad'`plothist'"!="")
    if `pairwork' {
        scalar __tmo_screened = scalar(__tmo_s_sel_screened)
        scalar __tmo_retained = scalar(__tmo_s_sel_screened)+scalar(__tmo_s_sel_nprotected)
        scalar __tmo_unsupported = scalar(__tmo_s_sel_unsupported)
        assert scalar(__tmo_screened)>=0 & scalar(__tmo_screened)<=scalar(__tmo_s_sel_neligible)
        assert scalar(__tmo_unsupported)>=0 & scalar(__tmo_unsupported)<=scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2
        scalar __tmo_s_offdP = scalar(__tmo_retained)/(scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2)
        scalar __tmo_s_offdPnocl = scalar(__tmo_screened)/scalar(__tmo_s_sel_neligible)
    }
    if `pairrecompute' {
        clear
        if "`savedyad'"!="" {
            foreach v in id1 id2 corr diagonal contribution retained selected_contribution missing_support shared_support histogram_weight protected eligible screened_included {
                quietly gen double `v' = .
            }
        }
        mata: tmo_selection_pairs("`savedyad'"!="", "`plothist'"!="", `nbins', "`kernel'", st_local("libpath"))
        assert scalar(__tmo_screened)==scalar(__tmo_s_sel_screened)
        assert scalar(__tmo_unsupported)==scalar(__tmo_s_sel_unsupported)
        if "`savedyad'"!="" {
            assert diagonal==(id1==id2)
            assert retained==1 if diagonal
            assert protected==1 if diagonal
            assert eligible==1-protected
            assert retained==1 if protected
            assert screened_included==eligible*(abs(corr)>=scalar(__tmo_s_sel_thres))
            assert retained==protected+screened_included
            assert histogram_weight==2*eligible
            assert selected_contribution==contribution*retained
            quietly gen byte diag_schema = 2
            quietly gen str16 diag_method = "tmo"
            quietly gen str16 selection_mode = "automatic"
            quietly gen str16 retention_rule = "`retention'"
            quietly gen str24 calibration_scope = "outside_protected"
            quietly gen double selected_index = scalar(__tmo_s_sel_index)
            quietly gen double threshold = scalar(__tmo_s_sel_thres)
            quietly gen double gamma0 = scalar(__tmo_s_sel_gamma0)
            label var contribution "Signed variance contribution; off-diagonals doubled"
            label var missing_support "Undefined pair normalization, scored zero off-diagonal"
            label var shared_support "Number of jointly observed projected auxiliary residuals"
            label var histogram_weight "Calibration multiplicity (eligible 2, protected 0)"
            label var protected "Always retained: diagonal, main cluster, or distance window"
            label var eligible "Outside the protected set; included in IQR calibration"
            label var screened_included "Eligible pair retained by the correlation threshold"
            quietly compress
            tmo_dyad_save, file(`"`filesuffix'_dyad.dta"') `stringid'
        }
    }
    if "`plothist'"!="" {
        * The eligible pairs' correlations and the mean-zero normal with TMO's
        * null SD sqrt(v0hat), drawn as in the paper.
        clear
        quietly set obs `nbins'
        quietly gen double corr = .
        quietly gen double frequency = .
        mata: st_store(.,("corr","frequency"),__tmo_hist)
        local binwidth = 2/`nbins'
        quietly gen double density = frequency/(scalar(__tmo_s_sel_calweight)*`binwidth')
        local nullsd : display %21.17g sqrt(scalar(__tmo_s_sel_v0))
        quietly gen double reference = normalden(corr/`nullsd')/`nullsd'
        local ytop = 0
        foreach v in density reference {
            quietly summarize `v' if abs(corr)<=.8, meanonly
            local ytop = max(`ytop',r(max))
        }
        local ystep = cond(`ytop'<=4,.5,cond(`ytop'<=8,1,2))
        local ytop = `ystep'*ceil(`ytop'/`ystep')
        twoway (bar density corr if abs(corr)<=.8, base(0) barwidth(`binwidth') color(midgreen%30)) ///
            (line reference corr if abs(corr)<=.8, sort lcolor(blue) lwidth(medthin)), ///
            xtitle("Correlation") ytitle("Density") ///
            xlabel(-.8(.2).8, format(%02.1f)) ylabel(0(`ystep')`ytop', format(%3.1f) angle(horizontal)) ///
            xscale(range(-.8 .8)) yscale(range(0 `ytop')) ///
            legend(off) scheme(s2color) graphregion(color(white) lcolor(white)) bgcolor(white) ///
            xsize(16) ysize(6.5)
        graph export "`filesuffix'_hist.pdf", as(pdf) replace
        graph export "`filesuffix'_hist.png", as(png) width(2400) replace
    }
    if "`plotq'`plotse'`saveest'`saveplotseest'"!="" {
        clear
        local ng = colsof(__tmo_etagrid_r)
        quietly set obs `ng'
        foreach v in index eta delta variance penalty objective {
            quietly gen double `v' = .
        }
        mata: st_store(.,("index","eta","delta","variance","penalty","objective"), ///
            ((1::cols(st_matrix("__tmo_etagrid_r"))),st_matrix("__tmo_etagrid_r")', ///
            st_matrix("__tmo_dgrid_r")',st_matrix("__tmo_vpath_r")', ///
            st_matrix("__tmo_xipath_r")',st_matrix("__tmo_objectivepath_r")'))
        if scalar(__tmo_need_counts) {
            quietly gen double screened_pairs = .
            quietly gen double eligible_pct = .
            quietly gen byte admissible = .
            mata: st_store(.,"screened_pairs",st_matrix("__tmo_screened_counts_r")')
            mata: st_store(.,"eligible_pct",st_matrix("__tmo_eligible_pctpath_r")')
            mata: st_store(.,"admissible",st_matrix("__tmo_admissiblepath_r")')
            assert screened_pairs>=0 & ///
                screened_pairs<=scalar(__tmo_s_sel_neligible)
            assert abs(eligible_pct-100*screened_pairs/ ///
                scalar(__tmo_s_sel_neligible))<1e-12
        }
        quietly gen double se = sqrt(variance) if variance>0
        quietly gen double se_ratio = se/scalar(__tmo_s_se)
        quietly gen byte selected = index==scalar(__tmo_s_sel_index)
        quietly gen byte nonpositive_variance = variance<=0
        assert eta>=eta[_n-1] if _n>1
        assert delta>=delta[_n-1] if _n>1
        assert objective==variance-penalty
        quietly count if selected
        assert r(N)==1
        assert variance==scalar(__tmo_s_sel_var) if selected
        if scalar(__tmo_need_counts) assert admissible==1 if selected
        quietly gen byte diag_schema = 2
        quietly gen str16 diag_method = "tmo"
        quietly gen str16 tmo_version = scalar(__tmo_version)
        quietly gen str16 selection_mode = "automatic"
        quietly gen str16 retention_rule = "`retention'"
        quietly gen str24 calibration_scope = "outside_protected"
        local selected_eta : display %21.17g scalar(__tmo_s_sel_eta)
        local selected_delta : display %9.5g scalar(__tmo_s_sel_thres)
        if "`plotq'"!="" {
            twoway (connected variance eta, lcolor(navy) mcolor(navy)) ///
                (connected penalty eta, lcolor(maroon) mcolor(maroon)) ///
                (connected objective eta, lcolor(forest_green) mcolor(forest_green)) ///
                (scatter objective eta if selected, msymbol(Dh) msize(large) mcolor(black)), ///
                xline(`selected_eta', lpattern(dash)) xtitle("eta") ytitle("Variance units") ///
                legend(order(1 "V0" 2 "xi (penalty)" 3 "V0 - xi" 4 "Returned selected index") position(6) rows(2) size(small)) ///
                note("Actual candidate grid; selected delta = `selected_delta'.")
            graph export "`filesuffix'_q.pdf", as(pdf) replace
        }
        if "`plotse'"!="" {
            quietly count if nonpositive_variance
            local badse = r(N)
            twoway (connected se_ratio eta, cmissing(n) lcolor(navy)) ///
                (scatter se_ratio eta if selected, msymbol(Dh) msize(large) mcolor(maroon)), ///
                xline(`selected_eta', lpattern(dash)) yline(1, lpattern(dot)) ///
                xtitle("eta") ytitle("Unpenalized target SE / original SE") legend(off) ///
                note("Actual candidate grid; selected delta = `selected_delta'. `badse' nonpositive variances have undefined SE.")
            graph export "`filesuffix'_se.pdf", as(pdf) replace
        }
        if "`saveest'`saveplotseest'"!="" {
            char _dta[tmo_outcome_normalized] "`=scalar(__tmo_outcome_normalized)'"
            quietly save "`filesuffix'_path.dta", replace
            if "`saveplotseest'"!="" quietly save "`filesuffix'_se_over_thres.dta", replace
            clear
            quietly mata: tmo_selection_anchor_export()
            quietly gen byte diag_schema = 2
            quietly gen str16 diag_method = "tmo"
            quietly gen str16 tmo_version = scalar(__tmo_version)
            quietly gen str16 selection_mode = "automatic"
            quietly gen str16 retention_rule = "`retention'"
            quietly gen str24 calibration_scope = "outside_protected"
            char _dta[tmo_outcome_normalized] "`=scalar(__tmo_outcome_normalized)'"
            quietly save "`filesuffix'_anchor_path.dta", replace
        }
    }
    if "`saveest'"!="" {
        clear
        quietly set obs 1
        foreach v in beta tmo_se lb ub df_r scpc_se tmo_scpc_se {
            quietly gen double `v' = scalar(__tmo_s_`v')
        }
        quietly gen double orig_se = scalar(__tmo_s_se)
        quietly gen double variance = scalar(__tmo_s_sel_var)
        quietly gen double selected_index = scalar(__tmo_s_sel_index)
        quietly gen double eta = scalar(__tmo_s_sel_eta)
        quietly gen double threshold = scalar(__tmo_s_sel_thres)
        quietly gen double grid_max = scalar(__tmo_s_sel_gmax)
        quietly gen double penalty = scalar(__tmo_s_sel_penalty)
        quietly gen double objective = scalar(__tmo_s_sel_objective)
        quietly gen double pilot = scalar(__tmo_s_sel_pilot)
        quietly gen double gamma0 = scalar(__tmo_s_sel_gamma0)
        quietly gen double logfactor = scalar(__tmo_s_sel_logf)
        quietly gen double q25 = scalar(__tmo_s_sel_q25)
        quietly gen double q75 = scalar(__tmo_s_sel_q75)
        quietly gen double max_eligible_pct = ///
            scalar(__tmo_s_sel_max_eligible_pct)
        quietly gen double v0hat = scalar(__tmo_s_sel_v0)
        quietly gen double effective_null_dimension = 1/scalar(__tmo_s_sel_v0)
        quietly gen double N_loc = scalar(__tmo_s_N)
        quietly gen double N = scalar(__tmo_s_N_obs)
        quietly gen double N_observed_rows = scalar(__tmo_s_N_rows)
        quietly gen double T = scalar(__tmo_s_T)
        quietly gen double N_clust = scalar(__tmo_s_N_clust)
        quietly gen double finite_sample_k = scalar(__tmo_s_finite_k)
        quietly gen double finite_sample_dof = scalar(__tmo_s_dof_adj)
        quietly gen double N_outcomes = scalar(__tmo_s_D)
        quietly gen double N_used = scalar(__tmo_nout)
        quietly gen double N_dropped = scalar(__tmo_ndropped)
        quietly gen double N_novar = scalar(__tmo_nnovar)
        quietly gen double incomplete = scalar(__tmo_hasmiss)
        quietly gen double pct_ge_thres = 100*scalar(__tmo_s_offdP)
        quietly gen double pct_ge_thres_nocl = 100*scalar(__tmo_s_offdPnocl)
        quietly gen double pct_protected = 100*scalar(__tmo_s_sel_nprotected)/(scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2)
        quietly gen double pct_screened = 100*scalar(__tmo_screened)/(scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2)
        quietly gen double N_protected_pairs = scalar(__tmo_s_sel_nprotected)
        quietly gen double N_eligible_pairs = scalar(__tmo_s_sel_neligible)
        quietly gen double retained_pairs = scalar(__tmo_retained)
        quietly gen double screened_pairs = scalar(__tmo_screened)
        quietly gen byte admissible = ///
            100*screened_pairs<=max_eligible_pct*N_eligible_pairs
        assert admissible==1
        quietly gen double missing_support_pairs = scalar(__tmo_unsupported)
        quietly gen double calibration_weight = scalar(__tmo_histweight)
        quietly gen double scpc_cv = ${scpc_cv}
        quietly gen byte diag_schema = 2
        quietly gen str16 diag_method = "tmo"
        quietly gen str16 tmo_version = scalar(__tmo_version)
        quietly gen str16 selection_mode = "automatic"
        quietly gen str16 retention_rule = "`retention'"
        quietly gen str24 calibration_scope = "outside_protected"
        quietly gen str24 kernel = "`kernel'"
        quietly gen str24 pilot_source = "`pilotsource'"
        quietly compress
        char _dta[tmo_outcome_normalized] "`=scalar(__tmo_outcome_normalized)'"
        quietly save "`filesuffix'_est.dta", replace
    }
    clear
    ereturn clear
    di as txt _n "TMO standard error for `xname': " as res %12.7g scalar(__tmo_s_tmo_se)
    di as txt "  eta = " as res %8.5g scalar(__tmo_s_sel_eta) as txt ", delta = " as res %9.5g scalar(__tmo_s_sel_thres) ///
        as txt ", selected index = " as res %4.0f scalar(__tmo_s_sel_index)
    di as txt "  effective null dimension = " as res %10.5g 1/scalar(__tmo_s_sel_v0) ///
        as txt "; regression df_r = " as res %10.5g scalar(__tmo_s_df_r)
    di as txt "  finite-sample factor n/(n-k) = " as res %10.7g scalar(__tmo_s_dof_adj) ///
        as txt "; k = " as res %8.0f scalar(__tmo_s_finite_k)
    di as txt "  `kernel' kernel; pilot from `pilotsource'; retained outcomes = " as res scalar(__tmo_nout)
    di as txt "  protected off-diagonal pairs (cluster/distance) = " as res %8.3f ///
        100*scalar(__tmo_s_sel_nprotected)/(scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2) as txt "% of all off-diagonal pairs"
    if `pairwork' {
        di as txt "  additionally kept by threshold = " as res %8.3f ///
            100*scalar(__tmo_screened)/(scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2) as txt "% of all off-diagonal pairs"
        di as txt "  total kept = " as res %8.3f 100*scalar(__tmo_s_offdP) ///
            as txt "%; threshold pass rate among eligible pairs = " as res %8.3f 100*scalar(__tmo_s_offdPnocl) as txt "%"
    }
    if scalar(__tmo_s_tmo_scpc_se)<. {
        di as txt "  Separate comparisons: SCPC SE = " as res %12.7g scalar(__tmo_s_scpc_se) ///
            as txt "; TMO-SCPC hybrid SE = " as res %12.7g scalar(__tmo_s_tmo_scpc_se)
        di as txt "  Hybrid uses the TMO-selected pairs; not optimized by TMO."
    }
    return scalar beta = scalar(__tmo_s_beta)
    return scalar orig_se = scalar(__tmo_s_se)
    return scalar tmo_se = scalar(__tmo_s_tmo_se)
    return scalar variance = scalar(__tmo_s_sel_var)
    return scalar lb = scalar(__tmo_s_lb)
    return scalar ub = scalar(__tmo_s_ub)
    return scalar eta = scalar(__tmo_s_sel_eta)
    return scalar threshold = scalar(__tmo_s_sel_thres)
    return scalar selected_index = scalar(__tmo_s_sel_index)
    return scalar grid_max = scalar(__tmo_s_sel_gmax)
    return scalar penalty = scalar(__tmo_s_sel_penalty)
    return scalar objective = scalar(__tmo_s_sel_objective)
    return scalar v0hat = scalar(__tmo_s_sel_v0)
    return scalar effective_null_dimension = 1/scalar(__tmo_s_sel_v0)
    return scalar pilot = scalar(__tmo_s_sel_pilot)
    return scalar gamma0 = scalar(__tmo_s_sel_gamma0)
    return scalar logfactor = scalar(__tmo_s_sel_logf)
    return scalar q25 = scalar(__tmo_s_sel_q25)
    return scalar q75 = scalar(__tmo_s_sel_q75)
    return scalar max_eligible_pct = scalar(__tmo_s_sel_max_eligible_pct)
    return scalar admissible = 1
    return scalar interior = scalar(__tmo_s_sel_interior)
    return scalar N_used = scalar(__tmo_nout)
    return scalar N_dropped = scalar(__tmo_ndropped)
    return scalar N_novar = scalar(__tmo_nnovar)
    return scalar N_depend = scalar(__tmo_ndepend)
    return scalar incomplete = scalar(__tmo_hasmiss)
    return scalar N = scalar(__tmo_s_N_obs)
    return scalar N_observed_rows = scalar(__tmo_s_N_rows)
    return scalar N_loc = scalar(__tmo_s_N)
    return scalar N_clust = scalar(__tmo_s_N_clust)
    return scalar N_outcomes = scalar(__tmo_s_D)
    return scalar T = scalar(__tmo_s_T)
    return scalar df_r = scalar(__tmo_s_df_r)
    return scalar finite_sample_k = scalar(__tmo_s_finite_k)
    return scalar finite_sample_dof = scalar(__tmo_s_dof_adj)
    return scalar pct_ge_thres = 100*scalar(__tmo_s_offdP)
    return scalar pct_ge_thres_nocl = 100*scalar(__tmo_s_offdPnocl)
    return scalar pct_protected = 100*scalar(__tmo_s_sel_nprotected)/(scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2)
    return scalar pct_screened = 100*scalar(__tmo_screened)/(scalar(__tmo_s_N)*(scalar(__tmo_s_N)-1)/2)
    return scalar N_protected_pairs = scalar(__tmo_s_sel_nprotected)
    return scalar N_eligible_pairs = scalar(__tmo_s_sel_neligible)
    return scalar retained_pairs = scalar(__tmo_retained)
    return scalar screened_pairs = scalar(__tmo_screened)
    return scalar missing_support_pairs = scalar(__tmo_unsupported)
    return scalar calibration_weight = scalar(__tmo_histweight)
    return scalar scpc_se = scalar(__tmo_s_scpc_se)
    return scalar tmo_scpc_se = scalar(__tmo_s_tmo_scpc_se)
    return scalar scpc_cv = ${scpc_cv}
    return scalar diag_schema = 2
    return scalar diagnostics_computed = `pairwork'
    return scalar pair_scores_recomputed = `pairrecompute'
    return scalar diagnostic_preparations = scalar(__tmo_diagnostic_preparations)
    return local diag_method "tmo"
    return local selection_mode "automatic"
    return local retention_rule "`retention'"
    return local calibration_scope "outside_protected"
    return local kernel "`kernel'"
    return local pilot_source "`pilotsource'"
    return matrix vpath = __tmo_vpath_r
    return matrix xipath = __tmo_xipath_r
    return matrix dgrid = __tmo_dgrid_r
    return matrix etagrid = __tmo_etagrid_r
    return matrix objectivepath = __tmo_objectivepath_r
    return matrix anchorpath = __tmo_anchorpath_r
    return matrix anchorstats = __tmo_anchorstats_r
    if scalar(__tmo_need_counts) {
        return matrix screened_countpath = __tmo_screened_counts_r
        return matrix eligiblepctpath = __tmo_eligible_pctpath_r
        return matrix admissiblepath = __tmo_admissiblepath_r
    }
end

capture program drop tmo_scpc_release
program define tmo_scpc_release
    capture confirm scalar __tmo_scpc_owned
    if !_rc {
        if scalar(__tmo_scpc_owned) {
            foreach m in Wfin distmat Wx Omsfin Omsxfin s permfin {
                capture mata: mata drop `m'
            }
            scalar __tmo_scpc_owned = 0
        }
    }
end

capture program drop tmo_cleanup
program define tmo_cleanup
    capture confirm scalar __tmo_diagnostic_ready
    if !_rc {
        if scalar(__tmo_diagnostic_ready)==1 {
            capture quietly javacall tmo.tmo releaseTargetScores, jars(tmo.jar)
        }
    }
    foreach m in __tmo_input_R __tmo_input_x __tmo_input_a __tmo_input_w __tmo_input_id __tmo_idlabels __tmo_input_names ///
        __tmo_R __tmo_Xs __tmo_W __tmo_keepcols __tmo_locid __tmo_paths __tmo_hist ///
        __tmo_etagrid __tmo_dgrid ///
        __tmo_pair_i __tmo_pair_j __tmo_pair_scores __tmo_pair_supported __tmo_pair_shared ///
        __tmo_Bcl __tmo_Bcoords __tmo_Bclids __tmo_Bcoordids __tmo_Bpairs __tmo_Browinfo __tmo_Bgroup ///
        __tmo_clmat __tmo_legacy_id __tmo_legacy_t __tmo_legacy_N __tmo_legacy_T ///
        __tmo_legacy_idfull __tmo_legacy_id_wide __tmo_legacy_xtilde __tmo_legacy_wgt ///
        __tmo_legacy_xtilde_wgt __tmo_legacy_Res __tmo_legacy_Res1 __tmo_legacy_Res1_x __tmo_legacy_Dn ///
        __tmo_legacy_CovEpsVec __tmo_legacy_ResXtildeVec __tmo_legacy_id_widerowvec __tmo_legacy_id_widecolvec ///
        __tmo_legacy_y_scpc __tmo_legacy_id_scpc __tmo_legacy_id_scpc_uniq ///
        __tmo_legacy_id_scpc_rowvec __tmo_legacy_id_scpc_colvec __tmo_legacy_Wfin_sum_vec ///
        __tmo_scpc_ids __tmo_scpc_basis ///
        __dst_id __dst_la __dst_lo __dst_i1 __dst_i2 __dst_la1 __dst_la2 __dst_lo1 __dst_lo2 {
        capture mata: mata drop `m'
    }
    foreach m in __tmo_screened_counts __tmo_screened_counts_r ///
        __tmo_eligible_pctpath __tmo_eligible_pctpath_r ///
        __tmo_admissiblepath __tmo_admissiblepath_r __tmo_paths __tmo_dgrid_out ///
        __tmo_vpath __tmo_xipath __tmo_dgrid_m __tmo_etagrid_m __tmo_anchorpath __tmo_anchorstats __tmo_objectivepath ///
        __tmo_vpath_r __tmo_xipath_r __tmo_dgrid_r __tmo_etagrid_r __tmo_anchorpath_r __tmo_anchorstats_r __tmo_objectivepath_r {
        capture matrix drop `m'
    }
    tmo_scpc_release
    capture scalar drop __tmo_scpc_owned
    capture scalar drop __tmo_scpc_neff
    capture scalar drop __tmo_scpc_normN
    capture scalar drop __tmo_need_counts
end

capture mata: mata drop tmo_ols_batch()
capture mata: mata drop tmo_aux_checks()
capture mata: mata drop tmo_lowtri_gram()
capture mata: mata drop tmo_selection_pairs()
capture mata: mata drop tmo_selection_corr()
capture mata: mata drop tmo_selection_anchor_export()
capture mata: mata drop tmo_build_protection()
capture mata: mata drop tmo_mask_guard()
mata:
void tmo_mask_guard(real scalar capacity, real scalar n)
{
    real scalar budget, need
    budget=st_numscalar("__tmo_host_free_bytes")
    need=32*capacity+128*n
    if (missing(need) | capacity>2^52 | (budget>0 & budget<. & need>.8*budget)) {
        errprintf("tmo: the protected-pair index would require approximately %g MiB during allocation.\n",ceil(need/2^20))
        errprintf("     Reduce the location count or protection window, or free memory before retrying.\n")
        _error(909)
    }
}
void tmo_build_protection(real scalar cutoff, real scalar miles)
{
    external real matrix __tmo_Bcl, __tmo_Bcoords, __tmo_Bpairs
    external real colvector __tmo_input_id, __tmo_Bclids, __tmo_Bcoordids, __tmo_Bgroup
    real matrix info
    real colvector ids,keep,jj,d,h,valid,sel
    real scalar n,i,c,q,k,capacity,nextcap,radius,rad
    ids=uniqrows(__tmo_input_id); n=rows(ids)
    assert(n==st_numscalar("__tmo_s_N"))
    __tmo_Bpairs=J(0,2,.)
    __tmo_Bgroup=J(0,1,.)
    if (cols(__tmo_Bcl)) {
        assert(rows(__tmo_Bcl)==n & rows(__tmo_Bclids)==n)
        assert(all(__tmo_Bclids:==ids))
        for (c=1;c<=cols(__tmo_Bcl);c++) {
            if (all(__tmo_Bcl[,c]:==__tmo_Bcl[1,c])) {
                errprintf("tmo: no unprotected off-diagonal pairs remain for calibration.\n")
                errprintf("     A main-specification cluster contains every location.\n")
                exit(459)
            }
        }
    }
    if (cutoff>0) {
        assert(rows(__tmo_Bcoords)==n & cols(__tmo_Bcoords)==2)
        assert(rows(__tmo_Bcoordids)==n & all(__tmo_Bcoordids:==ids))
    }
    if (!cols(__tmo_Bcl) & cutoff==0) return
    if (cols(__tmo_Bcl)==1 & cutoff==0) {
        jj=order(__tmo_Bcl,1); info=panelsetup(__tmo_Bcl[jj,],1)
        __tmo_Bgroup=J(n,1,.)
        for (c=1;c<=rows(info);c++)
            __tmo_Bgroup[jj[|info[c,1]\info[c,2]|]]=J(info[c,2]-info[c,1]+1,1,c)
        assert(!hasmissing(__tmo_Bgroup))
        return
    }
    st_numscalar("__tmo_host_free_bytes",.)
    if (st_local("jarpath")!="")
        stata("capture quietly javacall tmo.tmo hostMemory, jars("+char(34)+st_local("jarpath")+"/tmo.jar"+char(34)+")")
    // Spherical geodist convention. Pairwise temporary storage is only O(n);
    // the retained sparse index grows with protected pairs, not all dyads.
    radius=6371/(miles ? 1.609344 : 1); rad=pi()/180
    capacity=min((max((n,1024)),n*(n-1)/2)); q=0
    tmo_mask_guard(capacity,n)
    __tmo_Bpairs=J(capacity,2,.)
    for (i=2;i<=n;i++) {
        keep=J(i-1,1,0)
        for (c=1;c<=cols(__tmo_Bcl);c++)
            keep=keep:|(__tmo_Bcl[|1,c\i-1,c|]:==__tmo_Bcl[i,c])
        if (cutoff>0) {
            if (!hasmissing(__tmo_Bcoords[i,.])) {
                valid=rowsum(__tmo_Bcoords[|1,1\i-1,2|]:<.):==2
                jj=selectindex(valid)
                if (length(jj)) {
                    h=sin((__tmo_Bcoords[jj,1]:-__tmo_Bcoords[i,1]):*(rad/2)):^2 :+
                        cos(__tmo_Bcoords[jj,1]:*rad):*cos(__tmo_Bcoords[i,1]*rad):*
                        sin((__tmo_Bcoords[jj,2]:-__tmo_Bcoords[i,2]):*(rad/2)):^2
                    h=rowmin((J(rows(h),1,1),sqrt(h)))
                    d=2:*asin(h):*radius
                    keep[jj]=keep[jj]:|(d:<=cutoff)
                }
            }
        }
        sel=selectindex(keep); k=length(sel)
        if (!k) continue
        if (q+k>capacity) {
            nextcap=min((n*(n-1)/2,max((q+k,2*capacity))))
            tmo_mask_guard(nextcap,n)
            __tmo_Bpairs=__tmo_Bpairs\J(nextcap-capacity,2,.)
            capacity=nextcap
        }
        __tmo_Bpairs[|q+1,1\q+k,2|]=(J(k,1,i),sel)
        q=q+k
    }
    __tmo_Bpairs=q ? __tmo_Bpairs[|1,1\q,2|] : J(0,2,.)
}
void tmo_ols_batch(string scalar xv, string scalar yv, string scalar rv,
    string scalar wv, real scalar nocons, | string scalar checkrv, string scalar target)
{
    real matrix X,Y,B,R
    real colvector w,reference
    real scalar rank,scale,pos
    st_local("batch_ok","0")
    if (args()<6) checkrv=""
    if (args()<7) target=""
    X=st_data(.,tokens(xv)); Y=st_data(.,tokens(yv))
    if (!nocons) X=X,J(rows(X),1,1)
    w=wv=="" ? J(rows(X),1,1) : sqrt(st_data(.,wv))
    assert(!hasmissing(X) & !hasmissing(Y) & !hasmissing(w))
    if (rows(X)<cols(X)) return
    B=qrsolve(X:*w,Y:*w,rank)
    if (rank!=cols(X) | hasmissing(B)) return
    R=Y-X*B
    if (checkrv!="") {
        reference=st_data(.,checkrv)
        if (hasmissing(reference)) return
        scale=max((epsilon(1),max(abs(reference)),max(abs(Y[,1]))))
        if (max(abs(R[,1]-reference))>1e-10*scale) return
        // Preserve the original fit's target residual exactly, including rounding.
        R[,1]=reference
    }
    if (target!="") {
        pos=selectindex(tokens(xv):==target)
        if (length(pos)!=1) return
        st_matrix("__tmo_batch_bx",B[pos,])
    }
    st_store(.,tokens(rv),R)
    st_local("batch_ok","1")
}

// Auxiliaries whose spread after the regression is below 1e-8 of their spread
// before it (or that are constant before it) carry only rounding noise and are
// dropped. The rest are checked for exact linear dependence on complete rows:
// auxiliaries that rebuild the outcome of interest or the regressor stop the
// command; identities among the auxiliaries alone are reported.
void tmo_aux_checks(string scalar dvs, string scalar rvs, string scalar wv,
    string scalar idv, string scalar xv, real scalar misslim,
    string scalar bname, string scalar flagname)
{
    external string rowvector __tmo_input_names
    real matrix Y, E, Z, V, info, Eo
    real colvector a, id, o, ok, t, bb, r, c, xs, top
    real rowvector flags, miss, aux, sdz, b, L, keepz
    real scalar D, j, k, sw, my, sy, me, se, d, nc, gap, beta, sx, nnull
    string rowvector nm, rn
    string scalar dropped, list

    nm = tokens(dvs); rn = tokens(rvs); D = cols(nm)
    __tmo_input_names = nm
    st_numscalar("__tmo_ndepend", .)
    Y = st_data(., nm); E = st_data(., rn)
    a = wv=="" ? J(rows(Y),1,1) : st_data(., wv)
    flags = J(1,D,0); dropped = ""
    for (j=1; j<=D; j++) {
        ok = selectindex((E[,j]:<.) :& (Y[,j]:<.))
        if (!rows(ok)) continue
        sw = quadsum(a[ok])
        my = quadsum(a[ok]:*Y[ok,j])/sw
        sy = sqrt(quadsum(a[ok]:*(Y[ok,j]:-my):^2)/sw)
        me = quadsum(a[ok]:*E[ok,j])/sw
        se = sqrt(quadsum(a[ok]:*(E[ok,j]:-me):^2)/sw)
        if (sy<=1e-12*max((max(abs(Y[ok,j])),1e-300)) | se<=1e-8*sy) {
            flags[j] = 1
            if (j>1) dropped = dropped + " " + rn[j]
            if (j>1) printf("{txt}tmo: %s has no variation after the regression; dropped.\n", nm[j])
        }
        else if (j>1 & se<1e-5*sy) {
            printf("{txt}tmo: note: %s keeps only %9.2e of its spread after the regression.\n", nm[j], se/sy)
        }
    }
    if (flags[1]) {
        errprintf("tmo: the outcome of interest has no variation after the regression.\n")
        exit(459)
    }
    st_matrix(flagname, flags)
    st_numscalar("__tmo_nnovar", sum(flags))
    st_local("novar_resids", strtrim(dropped))

    // Auxiliaries that misslimit() will drop are left out of the check.
    id = st_data(., idv); o = order(id,1); info = panelsetup(id[o],1)
    Eo = E[o,]
    miss = colsum(panelsum(Eo:<., info):==0):/rows(info)
    aux = selectindex((!flags) :& (miss:<=misslim))
    if (sum(aux:>1)<2) {
        st_numscalar("__tmo_ndepend", 0)
        return
    }
    aux = select(aux, aux:>1)
    d = cols(aux)
    if (d<2) {
        st_numscalar("__tmo_ndepend", 0)
        return
    }
    ok = selectindex(rowmissing(E[,(1,aux)]):==0)
    nc = rows(ok)
    if (nc<d+2) {
        printf("{txt}tmo: too few complete rows to check the auxiliary outcomes for exact linear dependence.\n")
        return
    }
    Z = E[ok,aux]
    Z = Z :- mean(Z)
    sdz = sqrt(colsum(Z:^2):/nc)
    keepz = selectindex(sdz:>0)
    aux = aux[keepz]; Z = Z[,keepz]; sdz = sdz[keepz]; d = cols(aux)
    if (d<2) return
    Z = Z :/ sdz
    t = E[ok,1]; t = t :- mean(t)
    if (!(quadcross(t,t)>0)) return

    // Auxiliaries that add up to the outcome of interest.
    bb = qrsolve(Z, t)
    r = t - Z*bb
    gap = quadcross(r,r)/quadcross(t,t)
    if (gap<1e-6) {
        top = order(-abs(bb),1)
        list = invtokens(rowshape(nm[aux[rowshape(top[|1\min((8,d))|],1)]],1))
        if (gap<1e-10) {
            errprintf("tmo: auxiliary outcomes rebuild the outcome of interest exactly (1-R2 = %9.2e),\n", gap)
            errprintf("     led by: %s\n", list)
            errprintf("     Drop at least one of them from ylist().\n")
            exit(459)
        }
        printf("{txt}tmo: auxiliary outcomes nearly rebuild the outcome of interest (1-R2 = %9.2e), led by: %s\n", gap, list)
    }

    // Identities among the auxiliaries. One that also rebuilds x stops the command.
    symeigensystem(quadcross(Z,Z):/nc, V, L)
    b = st_matrix(bname)
    xs = st_data(., xv)
    sx = sqrt(quadvariance(xs[ok]))
    nnull = 0
    for (k=1; k<=d; k++) {
        if (L[k]>=1e-6) continue
        c = V[,k]
        top = order(-abs(c),1)
        top = top[|1\max((1,min((8,rows(select(top, abs(c[top]):>=0.05))))))|]
        list = invtokens(rowshape(nm[aux[rowshape(top,1)]],1))
        beta = .
        if (cols(b)==D) {
            if (!hasmissing(b[aux])) beta = quadsum(c :* (b[aux]' :/ sdz'))
        }
        if (L[k]<1e-10) {
            if (beta<. & abs(beta)*sx>1e-6) {
                errprintf("tmo: auxiliary outcomes rebuild the regressor exactly: %s\n", list)
                errprintf("     Drop at least one of them from ylist().\n")
                exit(459)
            }
            nnull++
            printf("{txt}tmo: note: exact linear dependence among auxiliary outcomes: %s\n", list)
            if (beta>=.) printf("{txt}tmo: note: x() coefficients unavailable; not checked against the regressor.\n")
        }
        else printf("{txt}tmo: note: near linear dependence among auxiliary outcomes: %s\n", list)
    }
    st_numscalar("__tmo_ndepend", nnull)
}
real colvector tmo_lowtri_gram(real matrix X)
{
    real scalar i,k,n
    real colvector out
    n=rows(X); out=J(n*(n+1)/2,1,.); k=1
    for (i=1;i<=n;i++) {
        out[|k\k+n-i|]=X[|i,1\n,cols(X)|]*X[i,.]'
        k=k+n-i+1
    }
    return(out)
}
// Target-screen scores in lowtri_vec order, from the kernel that selected.
real colvector tmo_selection_corr(string scalar kernel, string scalar libpath)
{
    external real matrix __tmo_Xs
    external real colvector __tmo_pair_i, __tmo_pair_j, __tmo_pair_scores
    real colvector corr, j
    real scalar n, i, k, first, last, nb
    n=rows(__tmo_Xs)
    if (kernel!="java")
        return(lowtri_vec(tmo_target_screen(__tmo_Xs,st_numscalar("__tmo_hasmiss"))))
    corr=J(n*(n+1)/2,1,.)
    stata("quietly javacall tmo.tmo prepareTargetScores, jars("+char(34)+libpath+"/tmo.jar"+char(34)+")")
    assert(st_numscalar("__tmo_diagnostic_api")==2 & st_numscalar("__tmo_diagnostic_ready")==1)
    for (first=1;first<=n;first=last+1) {
        last=first
        while (last<n & (last+1)*(last+2)/2-first*(first-1)/2<=65536) last++
        nb=last*(last+1)/2-first*(first-1)/2
        __tmo_pair_i=__tmo_pair_j=J(nb,1,.)
        k=1
        for (i=first;i<=last;i++) {
            __tmo_pair_i[|k\k+i-1|]=J(i,1,i)
            __tmo_pair_j[|k\k+i-1|]=(1::i)
            k=k+i
        }
        stata("quietly javacall tmo.tmo targetScores, jars("+char(34)+libpath+"/tmo.jar"+char(34)+")")
        assert(rows(__tmo_pair_scores)==nb)
        k=1
        for (i=first;i<=last;i++) {
            j=1::i
            corr[(j:-1):*(2*n:-j:+2):/2:+i:-j:+1]=__tmo_pair_scores[|k\k+i-1|]
            k=k+i
        }
    }
    stata("quietly javacall tmo.tmo releaseTargetScores, jars("+char(34)+libpath+"/tmo.jar"+char(34)+")")
    assert(st_numscalar("__tmo_diagnostic_ready")==0 & !hasmissing(corr))
    return(corr)
}
void tmo_selection_pairs(real scalar savepairs, real scalar histogram, real scalar nbins,
    string scalar kernel, string scalar libpath)
{
    external real matrix __tmo_Xs, __tmo_W, __tmo_input_id, __tmo_hist
    external real colvector __tmo_pair_i, __tmo_pair_j, __tmo_pair_scores
    external real colvector __tmo_pair_supported, __tmo_pair_shared
    real matrix M,En,A,out,C
    real colvector ids,ci,den,ss,xx,flags,unsupported,weight,bprotect,eligible,screened,hweight
    real scalar n,p,i,j,k,pos,nkeep,nscreen,nsupport,hw,total,xxden,delta,b,h,first,last,nb,abssum
    n=rows(__tmo_Xs); p=cols(__tmo_Xs)
    ids=uniqrows(__tmo_input_id)
    assert(rows(ids)==n & n==st_numscalar("__tmo_s_N"))
    if (kernel!="java") {
        tmo_project(__tmo_Xs,1,En,M)
        C=tmo_target_screen(__tmo_Xs,st_numscalar("__tmo_hasmiss"))
        st_numscalar("__tmo_diagnostic_preparations",1)
    }
    else {
        st_numscalar("__tmo_diagnostic_api",0)
        stata("quietly javacall tmo.tmo prepareTargetScores, jars("+char(34)+libpath+"/tmo.jar"+char(34)+")")
        assert(st_numscalar("__tmo_diagnostic_api")==2 & st_numscalar("__tmo_diagnostic_ready")==1)
        assert(st_numscalar("__tmo_diagnostic_preparations")==1)
    }
    delta=st_numscalar("__tmo_s_sel_thres")
    xxden=st_numscalar("__tmo_xx")^2
    nkeep=nscreen=nsupport=hw=total=abssum=0; pos=1
    if (savepairs) st_addobs(n*(n+1)/2)
    __tmo_hist=J(histogram ? nbins : 0,2,0)
    if (histogram) __tmo_hist[,1]=-1:+((1::nbins):-.5):*(2/nbins)
    for (first=1;first<=n;first=last+1) {
      last=first
      while (last<n & (last+1)*(last+2)/2-first*(first-1)/2<=65536) last++
      nb=last*(last+1)/2-first*(first-1)/2
      if (kernel=="java") {
        __tmo_pair_i=__tmo_pair_j=J(nb,1,.)
        k=1
        for (i=first;i<=last;i++) {
            __tmo_pair_i[|k\k+i-1|]=J(i,1,i)
            __tmo_pair_j[|k\k+i-1|]=(1::i)
            k=k+i
        }
        stata("quietly javacall tmo.tmo targetScores, jars("+char(34)+libpath+"/tmo.jar"+char(34)+")")
      }
      k=1
      for (i=first;i<=last;i++) {
       if (kernel=="java") {
        ci=__tmo_pair_scores[|k\k+i-1|]
        unsupported=1:-__tmo_pair_supported[|k\k+i-1|]
        ss=__tmo_pair_shared[|k\k+i-1|]
       }
       else {
        A=En[|1,1\i,p-1|]
        den=sqrt((A:^2)*M[i,.]')
        unsupported=(den:>=.):|(den:<=0)
        den=sqrt(M[|1,1\i,p-1|]*(En[i,.]:^2)')
        unsupported=unsupported:|(den:>=.):|(den:<=0)
        ss=M[|1,1\i,p-1|]*M[i,.]'
        ci=C[|1,i\i,i|]
        for (j=1;j<i;j++) if (unsupported[j])
            unsupported[j]=missing(tmo_pair_cosine(En[i,],En[j,],M[i,]:&M[j,]))
        ci[i]=1
        unsupported[i]=0
       }
        bprotect=tmo_protected_row(i)
        eligible=1:-bprotect
        screened=eligible:*(abs(ci):>=delta)
        flags=bprotect:|screened
        weight=J(i,1,2); weight[i]=1
        hweight=2:*eligible
        xx=(__tmo_W[|1,1\i,1|]:*__tmo_W[i,1]):*weight:/xxden
        total=total+sum(xx:*flags)
        abssum=abssum+sum(abs(xx))
        nkeep=nkeep+sum(flags)-1
        nscreen=nscreen+sum(screened)
        nsupport=nsupport+sum(unsupported)
        hw=hw+sum(hweight)
        if (histogram) for (j=1;j<=i;j++) {
            h=min((nbins,max((1,floor((ci[j]+1)*nbins/2)+1))))
            __tmo_hist[h,2]=__tmo_hist[h,2]+hweight[j]
        }
        if (savepairs) {
            out=(J(i,1,ids[i]),ids[|1\i|],ci,((1::i):==i),xx,flags,xx:*flags,unsupported,ss,hweight,bprotect,eligible,screened)
            st_store((pos::pos+i-1),tokens("id1 id2 corr diagonal contribution retained selected_contribution missing_support shared_support histogram_weight protected eligible screened_included"),out)
        }
        pos=pos+i
        k=k+i
      }
    }
    assert(pos-1==n*(n+1)/2 & hw==st_numscalar("__tmo_s_sel_calweight"))
    assert(nkeep>=0 & nkeep<=n*(n-1)/2)
    assert(nkeep==nscreen+st_numscalar("__tmo_s_sel_nprotected"))
    assert(nscreen>=0 & nscreen<=st_numscalar("__tmo_s_sel_neligible"))
    if (delta>1) assert(nscreen==0 & nkeep==st_numscalar("__tmo_s_sel_nprotected"))
    assert(abs(total-st_numscalar("__tmo_s_sel_var"))<=1e-10*max((1e-300,abssum,abs(st_numscalar("__tmo_s_sel_var")))))
    if (histogram) assert(sum(__tmo_hist[,2])==hw)
    st_numscalar("__tmo_retained",nkeep)
    st_numscalar("__tmo_screened",nscreen)
    st_numscalar("__tmo_unsupported",nsupport)
    st_numscalar("__tmo_histweight",hw)
    if (kernel=="java") {
        stata("quietly javacall tmo.tmo releaseTargetScores, jars("+char(34)+libpath+"/tmo.jar"+char(34)+")")
        assert(st_numscalar("__tmo_diagnostic_ready")==0)
    }
}
void tmo_selection_anchor_export()
{
    real matrix P,S,O
    real scalar a,g,k
    P=st_matrix("__tmo_anchorpath_r"); S=st_matrix("__tmo_anchorstats_r")
    O=J(rows(P)*cols(P),8,.); k=1
    for (a=1;a<=rows(P);a++) for (g=1;g<=cols(P);g++) {
        O[k,.]=(a,g,P[a,g],S[a,.]); k++
    }
    st_addvar("double",tokens("anchor index variance source_col sumsq phi omega denom"))
    st_addobs(rows(O)); st_store(.,.,O)
}
end

cap program drop tmo_cluster_options
program define tmo_cluster_options, rclass
    syntax [, VCE(string asis) CLuster(string asis) Robust *]
    return local options `"`options'"'
end

capture program drop tmo_fe_options
program define tmo_fe_options, rclass
    * Not `asis': that keeps the caller's compound quotes, and gettoken then
    * binds the whole option string as a single token, so nothing is stripped.
    syntax, TEXT(string)
    local rest `"`text'"'
    local cleaned
    while strtrim(`"`rest'"')!="" {
        gettoken part rest : rest, bind
        local open = strpos(`"`part'"',"(")
        local name = lower(strtrim(`"`part'"'))
        if `open'>0 local name = lower(strtrim(substr(`"`part'"',1,`open'-1)))
        local remove = (strlen("`name'")>=1 & "`name'"==substr("absorb",1,strlen("`name'"))) ///
            | (strlen("`name'")>=3 & "`name'"==substr("residuals",1,strlen("`name'")))
        if `remove' {
            if !`open' & substr(strtrim(`"`rest'"'),1,1)=="(" {
                gettoken discard rest : rest, bind
            }
        }
        else local cleaned `"`cleaned' `part'"'
    }
    return local options `"`cleaned'"'
end


* Position of a command's option comma: the first comma at parenthesis depth
* zero and outside quotes. strpos() finds commas inside if inrange(a,1,9) too.
cap mata: mata drop tmo_optcomma()
mata
    real scalar tmo_optcomma(string scalar s)
    {
        real scalar i, d, q, n
        string scalar c

        d = q = 0
        n = strlen(s)
        for (i=1; i<=n; i++) {
            c = substr(s,i,1)
            if (c==char(34)) q = !q
            else if (q) continue
            else if (c=="(") d++
            else if (c==")") d--
            else if (c=="," & d==0) return(i)
        }
        return(0)
    }
end

* Lower triangle including the diagonal, in vec() order. Column j contributes
* rows j..n, so neither vec() nor an N x N selector has to be formed.
cap mata: mata drop lowtri_vec()
cap mata: mata drop lowtri_outer()
cap mata: mata drop lowtri_ids()
cap mata: mata drop lowtri_samecl()
mata
    real colvector lowtri_vec(real matrix A)
    {
        real scalar n, j, p, m
        real colvector out

        n = rows(A)
        out = J(n*(n+1)/2,1,.)
        p = 1
        for (j=1; j<=n; j++) {
            m = n-j+1
            out[|p,1 \ p+m-1,1|] = A[|j,j \ n,j|]
            p = p+m
        }
        return(out)
    }

    // Lower triangle of g*g' without forming the outer product.
    real colvector lowtri_outer(real colvector g)
    {
        real scalar n, j, p, m
        real colvector out

        n = rows(g)
        out = J(n*(n+1)/2,1,.)
        p = 1
        for (j=1; j<=n; j++) {
            m = n-j+1
            out[|p,1 \ p+m-1,1|] = g[|j,1 \ n,1|]:*g[j]
            p = p+m
        }
        return(out)
    }

    // Same-cluster indicator over the same triangle, from one code per
    // location. Multiway: a pair shares a cluster if ANY dimension matches.
    real colvector lowtri_samecl(real matrix cl)
    {
        real scalar n, j, p, m, k
        real colvector out, blk

        n = rows(cl)
        out = J(n*(n+1)/2,1,0)
        p = 1
        for (j=1; j<=n; j++) {
            m = n-j+1
            blk = J(m,1,0)
            for (k=1; k<=cols(cl); k++) blk = blk :| (cl[|j,k \ n,k|]:==cl[j,k])
            out[|p,1 \ p+m-1,1|] = blk
            p = p+m
        }
        return(out)
    }

    void lowtri_ids(real matrix id, real matrix rowid, real matrix colid)
    {
        real scalar n, j, p, m

        n = rows(id)
        rowid = J(n*(n+1)/2,1,.)
        colid = J(n*(n+1)/2,1,.)
        p = 1
        for (j=1; j<=n; j++) {
            m = n-j+1
            rowid[|p,1 \ p+m-1,1|] = id[|j,1 \ n,1|]
            colid[|p,1 \ p+m-1,1|] = J(m,1,id[j])
            p = p+m
        }
    }
end

capture mata: mata drop tmo_scpc_capture()
capture mata: mata drop tmo_scpc_stream()
mata:
void tmo_scpc_capture(string scalar idvar, string scalar groupvar)
{
    external real matrix Wfin, __tmo_scpc_basis
    external real vector permfin
    external real colvector __tmo_scpc_ids
    real matrix B, info
    real colvector ids, groups, p, ord, first
    real scalar m,k
    m=rows(Wfin); p=vec(permfin)
    ids=st_data(.,idvar); groups=st_data(.,groupvar)
    assert(m>1 & cols(Wfin)>1 & rows(p)==m & !hasmissing(Wfin))
    assert(all(sort(p,1):==(1::m)))
    assert(!hasmissing(ids) & !hasmissing(groups) & min(groups)==1 & max(groups)==m)
    assert(all(groups:==floor(groups)))
    B=J(m,cols(Wfin)-1,.)
    B[p,.]=Wfin[,2..cols(Wfin)]
    ord=order(ids,1); ids=ids[ord]; groups=groups[ord]
    info=panelsetup(ids,1); first=info[,1]
    for (k=1;k<=rows(info);k++) {
        if (any(groups[|info[k,1]\info[k,2]|]:!=groups[info[k,1]])) {
            errprintf("tmo: SCPC must cluster repeated observations at the location level or coarser.\n")
            _error(459)
        }
    }
    __tmo_scpc_ids=ids[first]
    __tmo_scpc_basis=B[groups[first],.]
    st_numscalar("__tmo_scpc_neff",m)
}

void tmo_scpc_stream(real colvector corr, real colvector id,
    real colvector x, real colvector eps, real scalar uncond)
{
    external real matrix __tmo_scpc_basis, __tmo_Bcl
    external real colvector __tmo_scpc_ids
    real matrix info, U, totals
    real colvector ids,g,ys,j,pos,c,base,keep,screen,mult,xx,rr
    real scalar n,q,i,k,cidx,den,normn,delta,orig,vx,vr,vall,variance
    info=panelsetup(id,1); ids=id[info[,1]]; n=rows(ids)
    q=cols(__tmo_scpc_basis); normn=st_numscalar("__tmo_scpc_normN")
    assert(rows(corr)==n*(n+1)/2 & q>0 & rows(__tmo_scpc_ids)==rows(__tmo_scpc_basis))
    den=x'*x
    assert(den>0 & den<.)
    g=panelsum(x:*eps,info)
    ys=panelsum((sqrt(normn)/den):*eps:*x,info)
    U=J(n,q,0); k=1
    for (i=1;i<=n;i++) {
        if (k<=rows(__tmo_scpc_ids)) {
            if (ids[i]==__tmo_scpc_ids[k]) {
                U[i,.]=__tmo_scpc_basis[k,.]
                k++
            }
        }
    }
    assert(k==rows(__tmo_scpc_ids)+1)
    delta=st_numscalar("__tmo_s_sel_thres")
    totals=J(n,5,0)
    for (i=1;i<=n;i++) {
        j=1::i
        pos=(j:-1):*(2*n:-j:+2):/2:+i:-j:+1
        c=corr[pos]
        base=J(i,1,0); base[i]=1
        for (cidx=1;cidx<=cols(__tmo_Bcl);cidx++)
            base=base:|(__tmo_Bcl[|1,cidx\i,cidx|]:==__tmo_Bcl[i,cidx])
        keep=tmo_protected_row(i)
        screen=(c:<.):&(abs(c):>=delta):&!keep
        keep=keep:|screen
        mult=J(i,1,2); mult[i]=1
        xx=(g[|1\i|]:*g[i]):/(den^2):*mult
        rr=((U[|1,1\i,q|]*U[i,.]'):/q):*ys[i]:*ys[|1\i|]:*mult
        totals[i,.]=(quadcolsum(xx:*base),quadcolsum(xx:*keep),
            quadcolsum(rr:*!keep),quadcolsum(rr),quadcolsum(abs(xx)))
    }
    orig=quadcolsum(totals[,1]); vx=quadcolsum(totals[,2])
    vr=quadcolsum(totals[,3]); vall=quadcolsum(totals[,4])
    // The TMO part keeps the selected pairs: it is the selected variance before n/(n-k).
    if (abs(vx*st_numscalar("__tmo_s_dof_adj")-st_numscalar("__tmo_s_sel_var"))>
        1e-10*max((1e-300,quadcolsum(totals[,5])*st_numscalar("__tmo_s_dof_adj"),
        abs(st_numscalar("__tmo_s_sel_var"))))) {
        errprintf("tmo: the SCPC hybrid does not reproduce the selected TMO pairs.\n")
        _error(459)
    }
    if (!uncond & rows(__tmo_scpc_ids)==n & abs(st_numscalar("__tmo_s_scpc_se")-sqrt(vall))>1e-5) {
        errprintf("SCPC SE does not match\n")
        _error(459)
    }
    if (orig<=0 | missing(orig)) {
        errprintf("the original specification's variance estimate is not positive\n")
        _error(459)
    }
    variance=(vx+vr)*st_numscalar("__tmo_s_se")^2/orig
    if (variance<=0) errprintf("warning: the TMO-SCPC hybrid variance is not positive.\n")
    st_numscalar("__tmo_s_tmo_scpc_se",sqrt(variance))
    corr=J(0,1,.)
}
end


cap mata: mata drop sandwich_crosssec()
mata
    void sandwich_crosssec(real matrix id_widerowvec, real matrix id_widecolvec, real matrix xtilde, real matrix Res1, real matrix ResXtildeVec)
	{
		// (xx')*(ee') elementwise is the outer product of x:*e, so only the
		// lower triangle of one outer product is needed
        ResXtildeVec = id_widerowvec, id_widecolvec, lowtri_outer(xtilde:*Res1)
    }
end

cap mata: mata drop sandwich_panel_collapse()
mata
    void sandwich_panel_collapse(real scalar nloc, real matrix id, real matrix xtilde,
                                 real matrix Res1, real matrix id_widerowvec,
                                 real matrix id_widecolvec, real matrix ResXtildeVec)
    {
        real colvector g
        real matrix info

        // Summing over all time-pairs (s,s') within a location-pair factorizes:
        //   sum_{s,s'} (x_is e_is)(x_js' e_js') = g_i * g_j,  g_i = sum_s x_is e_is
        // so only the location-level vector g is needed, not the NT x NT array.
        // Summing by panel boundaries rather than by reshape keeps this correct
        // when the panel is unbalanced. Missing cells contribute zero.
        info = panelsetup(id, 1)
        assert(rows(info)==nloc)
        g = panelsum(editmissing(xtilde:*Res1, 0), info)

        ResXtildeVec = id_widerowvec, id_widecolvec, lowtri_outer(g)
    }
end

* Divide by (xtilde_wgt' * xtilde)^2; double off-diagonals stored only once.
cap mata: mata drop resxtildenorm()
mata
    void resxtildenorm(real matrix ResXtildeVec, real matrix xtilde, real matrix xtilde_wgt)
	{
		real matrix denom, offdiag

        denom = 1/(xtilde_wgt'*xtilde)^2
		offdiag = ResXtildeVec[,1]:!=ResXtildeVec[,2]

        ResXtildeVec[,3] = ResXtildeVec[,3]:*denom
        ResXtildeVec[,3] = ResXtildeVec[,3] + (ResXtildeVec[,3]:*offdiag)
    }
end

cap program drop load_data
program define load_data,
    syntax, nloc(int) [samecl] [distpath(str)] [savedyad] [filesuffix(str)] [scpcwfin(str)] [scpcy1path(str)] [scpcy2path(str)] [CACHEMODE(string)] [STRINGID]

    clear
    char _dta[tmo_outcome_normalized] "`=scalar(__tmo_outcome_normalized)'"

    * Both arrays use the same lower-triangle order, so no merge is needed.
    mata: assert(rows(__tmo_legacy_ResXtildeVec)==rows(__tmo_legacy_CovEpsVec))
    mata: st_local("obsN",strofreal(rows(__tmo_legacy_ResXtildeVec),"%50.0f"))

    gen double id1=.
    gen double id2=.
    gen double xxresxx=.
    gen double corr=.

    qui set obs `obsN'

    * Separate stores avoid allocating an M x 4 double temporary.
    mata: st_store(., ("id1","id2","xxresxx"), __tmo_legacy_ResXtildeVec)
    mata: st_store(., "corr", __tmo_legacy_CovEpsVec)
    mata: __tmo_legacy_ResXtildeVec = __tmo_legacy_CovEpsVec = J(0,0,.)

    * st_store() silently writes missing for values a variable cannot hold
    qui count if missing(id1) | missing(id2)
    if r(N) {
        di as error "location ids could not be stored -- are they larger than 2^53?"
        exit 459
    }

    if _N!=(`nloc'^2+`nloc')/2 {
        di as error "Number of location pairs is incorrect"
        exit 459
    }

    if "`samecl'"!="" {
        gen byte same_cl = 0
        mata: st_store(., "same_cl", lowtri_samecl(__tmo_clmat))
    }

    if "`distpath'"!="" {
        qui merge 1:1 id1 id2 using `distpath', assert(3) nogen
    }

    if "`scpcwfin'"!="" {
        qui merge 1:1 id1 id2 using `scpcwfin', assert(1 3) nogen
        qui merge m:1 id1 using `scpcy1path', assert(1 3) nogen
        qui merge m:1 id2 using `scpcy2path', assert(1 3) nogen
        qui gen double ryyr = (Wfin*y_scpc1*y_scpc2)*(1 + (id1!=id2))
    }

    qui compress

    if "`savedyad'"!="" {
        quietly gen byte diag_schema = 1
        quietly gen str16 diag_method = "legacy_override"
        quietly gen str16 selection_mode = "`cachemode'"
        * Only a custom-threshold cache carries the selection screen in corr.
        if "`cachemode'"=="fixed_threshold" char _dta[tmo_screen] "selection"
        else char _dta[tmo_screen] "none"
        tmo_dyad_save, file(`"`filesuffix'_dyad.dta"') `stringid'
    }
end

capture program drop tmo_dyad_save
program define tmo_dyad_save
    syntax, FILE(string) [STRINGID]
    char _dta[tmo_outcome_normalized] "`=scalar(__tmo_outcome_normalized)'"
    char _dta[tmo_version] "`=scalar(__tmo_version)'"
    if "`stringid'"=="" {
        quietly save `"`file'"', replace
        exit
    }
    preserve
        rename (id1 id2) (id1_code id2_code)
        quietly generate strL id1 = ""
        quietly generate strL id2 = ""
        mata: tmo_dyad_labels(0)
        assert !missing(id1) & !missing(id2)
        assert (id1==id2)==(id1_code==id2_code)
        label variable id1_code "Internal numeric location code for id1"
        label variable id2_code "Internal numeric location code for id2"
        char _dta[tmo_id_encoding] "string_codes_v1"
        order id1 id2 id1_code id2_code
        quietly compress id1 id2
        quietly save `"`file'"', replace
    restore
end

capture program drop tmo_cache_normalization
program define tmo_cache_normalization
    local cached : char _dta[tmo_outcome_normalized]
    if `"`cached'"'=="" local cached 1
    if !inlist(`"`cached'"',"0","1") {
        di as error "tmo: cached outcome-normalization metadata is invalid; regenerate the cache."
        exit 459
    }
    if real(`"`cached'"')!=scalar(__tmo_outcome_normalized) {
        di as error "tmo: cached outcome normalization differs from the requested mode; regenerate the cache."
        exit 459
    }
end

capture program drop tmo_dyad_cache_ids
program define tmo_dyad_cache_ids
    syntax, NLOC(integer) [STRINGID]
    local encoding : char _dta[tmo_id_encoding]
    if "`stringid'"=="" {
        capture confirm numeric variable id1 id2
        if _rc | "`encoding'"!="" {
            di as error "tmo: cached location-ID encoding differs from the current numeric idvar()."
            exit 459
        }
        exit
    }
    if "`encoding'"!="string_codes_v1" {
        di as error "tmo: string idvar() requires a dyad cache with its original string-ID mapping."
        exit 459
    }
    capture confirm string variable id1 id2
    local badids = _rc
    capture confirm numeric variable id1_code id2_code
    if `badids' | _rc {
        di as error "tmo: cached string location IDs or their numeric codes are missing."
        exit 459
    }
    mata: tmo_dyad_labels(1)
    capture assert _N==`nloc'*(`nloc'+1)/2 & id1_code>=id2_code
    local badids = _rc
    capture gisid id1_code id2_code
    if `badids' | _rc {
        di as error "tmo: cached dyads do not contain the complete unique location-pair triangle."
        exit 459
    }
    drop id1 id2
    rename (id1_code id2_code) (id1 id2)
end

capture mata: mata drop tmo_dyad_labels()
mata:
void tmo_dyad_labels(real scalar verify)
{
    external string colvector __tmo_idlabels
    real scalar p, j
    real colvector obs
    real matrix codes
    string scalar name
    string colvector labels

    assert(rows(__tmo_idlabels)>0)
    for(p=1;p<=st_nobs();p=p+50000) {
        obs=p::min((st_nobs(),p+49999))
        codes=st_data(obs,("id1_code","id2_code"))
        if(any(codes:>=.) | any(codes:!=floor(codes)) |
            min(codes)<1 | max(codes)>rows(__tmo_idlabels)) {
            errprintf("tmo: dyad location codes do not match the current string-ID mapping.\n")
            _error(459)
        }
        for(j=1;j<=2;j++) {
            name="id"+strofreal(j)
            labels=__tmo_idlabels[codes[.,j]]
            if(verify) {
                if(any(st_sdata(obs,name):!=labels)) {
                    errprintf("tmo: cached location labels differ from the current string-ID mapping.\n")
                    _error(459)
                }
            }
            else {
                st_sstore(obs,name,labels)
            }
        }
    }
}
end

cap program drop mkkeep
program define mkkeep
    syntax, nloc(int) [dist_cutoff(real 0)]

    cap confirm var same_cl
    if !_rc  local base (same_cl)
    else     local base (id1==id2)

    local keep (`base')
    cap confirm var dist
    if !_rc & `dist_cutoff'>0  local keep (`base' | (dist<=`dist_cutoff'))

    cap drop origpair
    cap drop keeppair
    cap drop offdiag
    cap drop corr_abs

    qui gen byte origpair = (`base')
    qui gen byte keeppair = (`keep')
    * A pair missing either endpoint's coordinates cannot enter the Conley window
    * and remains subject to screening.
    if `dist_cutoff'>0 {
        cap confirm var dist
        if !_rc {
            qui count if missing(dist)
            if r(N) {
                di as txt "note: `r(N)' location pairs have no distance " _c
                di as txt "(missing coordinates); they are screened, not kept"
            }
        }
    }

    qui gen byte offdiag = !missing(corr) & !keeppair
    qui gen double corr_abs = abs(corr)

    qui gstats sum xxresxx if origpair, meanonly
    scalar __tmo_s_v_orig = r(sum)

    qui count if id1!=id2
    scalar __tmo_s_offdN = r(N)
    if scalar(__tmo_s_offdN)!=`nloc'*(`nloc'-1)/2 {
        di as error "Number of off-diagonal location pairs is incorrect"
        exit 459
    }
    * denominator for the share above the threshold: the pairs the threshold is
    * actually applied to, so it excludes the Conley window as the numerator does
    qui count if !keeppair
    scalar __tmo_s_offdNnocl = r(N)
end

cap program drop est_tmo_se
program define est_tmo_se,
    syntax, [custom_thres(real -9)] [thresholdoff] [scpc_uncond]

    if `custom_thres'!=-9 {
        scalar __tmo_s_thres  = `custom_thres'
        scalar __tmo_s_df_iqr = .
        scalar __tmo_s_sd_iqr = .
    }

    if "`thresholdoff'"!="" {
        scalar __tmo_s_thres  = .
        scalar __tmo_s_df_iqr = .
        scalar __tmo_s_sd_iqr = .
    }

    cap confirm var Wfin
    if !_rc {
        local scpc = 1
        * Conditional SCPC reconstructs the native reference variance when coordinates are complete.
        qui count if missing(ryyr)
        if r(N)==0 & "`scpc_uncond'"=="" {
            qui gstats sum ryyr, meanonly
            if abs(scalar(__tmo_s_scpc_se) - sqrt(r(sum)))>1e-5 {
                di as error "SCPC SE does not match"
                exit 459
            }
        }
    }
    else {
        local scpc = 0
    }

    qui gstats sum xxresxx if (offdiag & corr_abs>=scalar(__tmo_s_thres)) | keeppair, meanonly
    local xxresxx_sum = r(sum)
    if `scpc'==1 {
        qui gstats sum ryyr if !(offdiag & corr_abs>=scalar(__tmo_s_thres)) & !keeppair, meanonly
        local ryyr_sum = r(sum)
    }
    else {
        local ryyr_sum = 0
    }
    * Signed contributions can sum to a negative variance; retain the other results
    * and report the resulting missing SE.
    if (`xxresxx_sum'+`ryyr_sum')<=0 {
        di as error "warning: the thresholded variance estimate is not positive (" ///
                    %10.3e (`xxresxx_sum'+`ryyr_sum') ")"
        di as error "         negative variance gives a missing SE; zero variance gives a zero SE"
    }
	scalar __tmo_s_tmo_variance = (`xxresxx_sum'+`ryyr_sum')*scalar(__tmo_s_dof_adj)
	scalar __tmo_s_tmo_se = sqrt(scalar(__tmo_s_tmo_variance))

	qui count if (offdiag & corr_abs>=scalar(__tmo_s_thres)) | (keeppair & id1!=id2)
	scalar __tmo_s_offdP = r(N)/scalar(__tmo_s_offdN)

	qui count if offdiag & corr_abs>=scalar(__tmo_s_thres)
	scalar __tmo_s_offdPnocl = r(N)/scalar(__tmo_s_offdNnocl)
end

cap program drop tmo_save
program define tmo_save,
    syntax, FILEsuffix(str) MODE(string)

    preserve
        clear

        qui set obs 1

        qui gen double beta = scalar(__tmo_s_beta)
        qui gen double orig_se = scalar(__tmo_s_se)
        qui gen double tmo_se = scalar(__tmo_s_tmo_se)
        quietly gen double variance = scalar(__tmo_s_tmo_variance)
        qui gen double lb = scalar(__tmo_s_lb)
        qui gen double ub = scalar(__tmo_s_ub)
        qui gen double threshold = scalar(__tmo_s_thres)
        qui gen double pct_ge_thres = scalar(__tmo_s_offdP)*100
        qui gen double pct_ge_thres_nocl = scalar(__tmo_s_offdPnocl)*100
        qui gen double T = scalar(__tmo_s_T)
        qui gen double N_loc = scalar(__tmo_s_N)
        qui gen double N = scalar(__tmo_s_N_obs)
        qui gen double N_clust = scalar(__tmo_s_N_clust)
        qui gen double N_outcomes = scalar(__tmo_s_D)
        qui gen double N_novar = scalar(__tmo_nnovar)
        qui gen double dof = scalar(__tmo_s_df_iqr)
        qui gen double finite_sample_dof = scalar(__tmo_s_dof_adj)
        qui gen double df_r = scalar(__tmo_s_df_r)
        qui gen double scpc_cv = ${scpc_cv}
        quietly gen double scpc_se = scalar(__tmo_s_scpc_se)
        quietly gen double tmo_scpc_se = scalar(__tmo_s_tmo_scpc_se)
        quietly gen byte diag_schema = 1
        quietly gen str16 diag_method = "legacy_override"
        quietly gen str16 tmo_version = scalar(__tmo_version)
        quietly gen str16 selection_mode = "`mode'"

        qui compress
        char _dta[tmo_outcome_normalized] "`=scalar(__tmo_outcome_normalized)'"
        save "`filesuffix'_est.dta", replace
    restore

end

* Target dyads cannot identify the auxiliary penalty. Each anchor's variance
* path is computed directly from the residual matrix.
* One row per location is required. A panel specification is reduced by
* the same collapse the sandwich uses, w_i = sum_t xtilde_it eps_it, and
* the screen then runs on those collapsed quantities.
cap program drop tmo_estimate
program define tmo_estimate, rclass
    syntax [, ETA(numlist ascending >0 <=1 min=2) PILOT(real -9) ///
        PCLust(string) CLPILOT(name) NONORMalize MISSlimit(real 0.1) ///
        THREADS(int 4) LIBpath(str) DOFADJ(real 1) ///
        MAXELIGIBLEPCT(real 100) PROTECTed NOVAR(name) ]

    if missing(`dofadj') | `dofadj'<=0 {
        di as error "tmo: dofadj() must be finite and strictly positive."
        exit 198
    }
    scalar __tmo_s_dof_adj = `dofadj'
    if missing(`maxeligiblepct') | `maxeligiblepct'<0 | ///
            `maxeligiblepct'>100 {
        di as error "tmo: maxeligiblepct() must be between 0 and 100, inclusive."
        exit 198
    }
    scalar __tmo_max_eligible_pct = `maxeligiblepct'

    if "`nonormalize'"=="" local normalize normalize

    * -9 is the not-supplied sentinel, as for threshold() above. A pilot
    * variance must be positive, so no legitimate value collides with it.
    tempname pilot_value
    scalar `pilot_value' = .
    if `pilot'!=-9 {
        if missing(`pilot') | `pilot'<=0 {
            di as error "tmo: pilot() must be positive"
            exit 198
        }
        if "`pclust'"!="" {
            di as error "tmo: pilot() and pclust() are alternatives"
            exit 198
        }
        scalar `pilot_value' = `pilot'
    }

    if "`eta'"=="" local eta 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1
    local ng : word count `eta'
    local etalast : word `ng' of `eta'
    if `etalast'!=1 {
        di as error "tmo: the eta grid must end at 1"
        exit 198
    }

    local psrc "supplied"
    if `pilot'==-9 & "`pclust'"=="" {
        if scalar(__tmo_s_se)<=0 | missing(scalar(__tmo_s_se)) {
            di as error "tmo: the original standard error is not usable as a pilot; supply pilot()"
            exit 459
        }
        scalar `pilot_value' = scalar(__tmo_s_se)^2
        local psrc "baseline"
    }
    else if "`pclust'"!="" {
        confirm scalar `clpilot'
        scalar `pilot_value' = scalar(`clpilot')
        local psrc "clustered"
    }
    if missing(scalar(`pilot_value')) | scalar(`pilot_value')<=0 {
        di as error "tmo: the pilot variance must be finite and strictly positive."
        exit 459
    }

    scalar __tmo_threads = `threads'
    scalar __tmo_v0hat = .
    scalar __tmo_ndropped = 0
    scalar __tmo_nnovar = 0
    scalar __tmo_hasmiss = 0
    scalar __tmo_nout = .
    scalar __tmo_q25 = .
    scalar __tmo_q75 = .
    scalar __tmo_masked = ("`protected'"!="")
    if `threads'<1 {
        di as error "tmo: threads() must be positive"
        exit 198
    }
    mata: tmo_prepare(st_numscalar("`pilot_value'"), "`normalize'", `misslimit', "`novar'")
    mata: tmo_mask_prepare(st_numscalar("__tmo_nn"), st_numscalar("__tmo_masked"))

    local kern "mata"
    scalar __tmo_java_memory_failure = 0
    scalar __tmo_kernel_api = 0
    scalar __tmo_calibration_error = 0
    scalar __tmo_mask_api = 0
    scalar __tmo_counts_api = 0
    capture confirm scalar __tmo_need_counts
    if _rc scalar __tmo_need_counts = 0
    if "`libpath'"!="" {
        capture confirm file "`libpath'/tmo.jar"
        if !_rc {
            adopath ++ "`libpath'"
            capture javacall tmo.tmo run, jars(tmo.jar)
            if !_rc local kern "java"
            if __tmo_java_memory_failure==1 {
                di as error "tmo: the compiled kernel ran out of memory; reduce threads() or increase available memory."
                exit 909
            }
            if __tmo_kernel_api!=4 {
                local kern "mata"
                di as txt "tmo: this compiled kernel predates the correlation-boundary fixes; using Mata."
            }
            if __tmo_kernel_api==4 & __tmo_calibration_error>0 {
                di as error "tmo: the eligible-pair screening IQR is not positive."
                di as error "     Use less extensive protection or auxiliary outcomes with more screening variation."
                exit 459
            }
            if __tmo_masked & __tmo_mask_api!=2 {
                local kern "mata"
                di as txt "tmo: this compiled kernel does not support protected pairs; using Mata."
            }
            if scalar(__tmo_need_counts) & scalar(__tmo_counts_api)!=1 local kern "mata"
        }
    }
    if "`kern'"=="mata" mata: tmo_paths_mata(st_numscalar("__tmo_masked"))

    if "`kern'"=="java" mata: tmo_select(1)
    else                 mata: tmo_select(0)

    if __tmo_gridmax>1 {
        di as txt "tmo: the largest candidate threshold is " %6.3f __tmo_gridmax ///
            ", above one."
        di as txt "     A correlation screen cannot retain any unprotected off-diagonal pair there."
        if "`nonormalize'"!="" di as txt "     Omit nonormalize to restore default outcome-column scaling."
    }

    return scalar tmo_var   = scalar(__tmo_var)
    return scalar tmo_se    = scalar(__tmo_se)
    return scalar threshold = scalar(__tmo_thres)
    return scalar eta       = scalar(__tmo_eta)
    return scalar grid_max  = scalar(__tmo_gridmax)
    return scalar v0hat     = scalar(__tmo_v0hat)
    return scalar dof       = 1/scalar(__tmo_v0hat)
    return scalar pilot     = scalar(__tmo_pilot)
    return scalar gamma0    = scalar(__tmo_gamma0)
    return scalar logfactor = scalar(__tmo_logfactor)
    return scalar finite_sample_dof = scalar(__tmo_s_dof_adj)
    return scalar interior  = scalar(__tmo_interior)
    return scalar selected_index = scalar(__tmo_selected_index)
    return scalar q25       = scalar(__tmo_q25)
    return scalar q75       = scalar(__tmo_q75)
    return scalar max_eligible_pct = scalar(__tmo_max_eligible_pct)
    return scalar N_loc     = scalar(__tmo_nn)
    return scalar N_protected_pairs = scalar(__tmo_nprotected)
    return scalar N_eligible_pairs = scalar(__tmo_neligible)
    return scalar calibration_weight = scalar(__tmo_calibration_weight)
    return scalar N_dropped = scalar(__tmo_ndropped)
    return scalar N_novar   = scalar(__tmo_nnovar)
    return scalar N_used     = scalar(__tmo_nout)
    return scalar incomplete = scalar(__tmo_hasmiss)
    if scalar(__tmo_need_counts) {
        return matrix screened_counts = __tmo_screened_counts
        return matrix eligiblepctpath = __tmo_eligible_pctpath
        return matrix admissiblepath = __tmo_admissiblepath
        return scalar unsupported_pairs = scalar(__tmo_unsupported_pairs)
    }
    if __tmo_ndropped>0 {
        di as txt "tmo: " __tmo_ndropped " auxiliary outcome(s) dropped for an incomplete residual"
    }
    if __tmo_nnovar>0 {
        di as txt "tmo: " __tmo_nnovar " auxiliary outcome(s) dropped for no variation after the regression"
    }
    return matrix vpath     = __tmo_vpath
    return matrix xipath    = __tmo_xipath
    return matrix dgrid     = __tmo_dgrid_m
    return matrix etagrid   = __tmo_etagrid_m
    matrix colnames __tmo_anchorstats = source_col sumsq phi omega denom
    return matrix anchorpath = __tmo_anchorpath
    return matrix anchorstats = __tmo_anchorstats
    return matrix objectivepath = __tmo_objectivepath
    return local  pilot_source "`psrc'"
    return local  kernel    "`kern'"
end

cap mata: mata drop tmo_prepare()
cap mata: mata drop tmo_screen_inputs()
cap mata: mata drop tmo_paths_mata()
cap mata: mata drop tmo_select()
cap mata: mata drop tmo_tie_index()
cap mata: mata drop tmo_anchor_paths()
cap mata: mata drop tmo_anchor_paths_miss()
cap mata: mata drop tmo_screen_a()
cap mata: mata drop tmo_xi()
cap mata: mata drop tmo_anchor_stats()
cap mata: mata drop tmo_xi_stats()
cap mata: mata drop tmo_iqr_all()
cap mata: mata drop tmo_kth()
cap mata: mata drop tmo_variance_path()
cap mata: mata drop tmo_target_screen()
cap mata: mata drop tmo_bound_corr()
cap mata: mata drop tmo_mask_prepare()
cap mata: mata drop tmo_protected_row()
cap mata: mata drop tmo_protected_matrix()
cap mata: mata drop tmo_iqr_eligible()
cap mata: mata drop tmo_dd_add()
cap mata: mata drop tmo_dd_product()
cap mata: mata drop tmo_baseline()
cap mata: mata drop tmo_pair_cosine()
cap mata: mata drop tmo_project()
cap mata: mata drop tmo_dd_sum()
cap mata: mata drop tmo_mask_ptr()
cap mata: mata drop tmo_baseline_dd()
cap mata: mata drop tmo_complete_setup()
cap mata: mata drop tmo_complete_screen()
cap mata: mata drop tmo_accumulate_pairs()

mata
void tmo_prepare(real scalar pilot, string scalar nrm,
                  real scalar misslim, | string scalar novarname)
{
    external real matrix __tmo_input_R, __tmo_W, __tmo_R
    external real colvector __tmo_input_w, __tmo_input_x, __tmo_input_a
    external real colvector __tmo_etagrid
    real scalar sumw, n, xx, gamma0, lf

    if (args()>=4) tmo_screen_inputs(nrm, misslim, novarname)
    else tmo_screen_inputs(nrm, misslim)

    // Weighted least squares puts the weight once in the denominator and
    // once in each side of the numerator. __tmo_W already carries it once, so
    // the denominator is sum_i a_i xtilde_i^2, matching resxtildenorm above
    // and the point estimate the command reports.
    xx = quadcross(__tmo_input_w, __tmo_input_x)
    sumw = quadcross(__tmo_input_a, J(rows(__tmo_input_a),1,1))
    if (xx<=0 | missing(xx) | sumw<=0 | missing(sumw)) _error(3498)

    // Put the pilot variance on the scale of the residualized treatment.
    // The residual scale makes the endpoint invariant to the units of the
    // outcome; analytic weights enter in the same way as in xx/sumw.
    gamma0 = quadcross(__tmo_input_a, __tmo_input_R[,1]:^2)/sumw
    if (gamma0<=0 | missing(gamma0)) _error(3498)

    // Apply n/(n-k) to every anchor variance before the penalty and objective.
    __tmo_W = __tmo_W:*sqrt(st_numscalar("__tmo_s_dof_adj"))
    n = rows(__tmo_R)

    if (pilot<=0 | missing(pilot)) _error(3498)

    lf = log(gamma0)-log(xx/sumw)-log(pilot)
    if (lf<=0 | missing(lf)) {
        errprintf("tmo: the scale-adjusted endpoint logarithm is not positive.\n")
        errprintf("     The pilot is too large relative to the residual and treatment scales.\n")
        _error(3498)
    }

    st_numscalar("__tmo_xx", xx)
    st_numscalar("__tmo_nn", n)
    st_numscalar("__tmo_pilot", pilot)
    st_numscalar("__tmo_gamma0", gamma0)
    st_numscalar("__tmo_logfactor", lf)
    st_numscalar("__tmo_corrscale", 1)
    // Each anchor pools all other outcomes. Its own projected residual is zero.
    st_numscalar("__tmo_auxonly", 0)
    st_numscalar("__tmo_endpoint", 1)

    __tmo_etagrid = strtoreal(tokens(st_local("eta")))'
    n = rows(__tmo_etagrid)
    if (n<2 | hasmissing(__tmo_etagrid) | min(__tmo_etagrid)<=0 |
        max(__tmo_etagrid)>1) _error(3498)
    if (any(__tmo_etagrid[|2\n|] :<= __tmo_etagrid[|1\n-1|])) _error(3498)
}

// Screening inputs shared by automatic selection and custom thresholds:
// location rows (score sums in panels), retained outcome columns, the
// unscaled sandwich terms __tmo_W, and the normalized screen matrix __tmo_Xs.
void tmo_screen_inputs(string scalar nrm, real scalar misslim, | string scalar novarname)
{
    external real matrix __tmo_input_R, __tmo_Xs, __tmo_W, __tmo_R
    external real colvector __tmo_input_w, __tmo_input_x, __tmo_input_a
    external real colvector __tmo_input_id, __tmo_locid
    external real rowvector __tmo_keepcols
    real matrix R, M, info, fl
    real colvector xt, mj
    real rowvector keepc, nob, novar, miss, rms
    real rowvector sd
    real scalar n, nr, nb, panel, j
    pointer() scalar pnames

    nr = rows(__tmo_input_R)
    if (nr<3 | cols(__tmo_input_R)<3 | rows(__tmo_input_x)!=nr |
        rows(__tmo_input_a)!=nr | rows(__tmo_input_w)!=nr | rows(__tmo_input_id)!=nr) _error(3498)
    if (hasmissing(__tmo_input_x) | hasmissing(__tmo_input_a) |
        hasmissing(__tmo_input_w) | hasmissing(__tmo_input_id) | min(__tmo_input_a)<=0) _error(3498)
    if (any(__tmo_input_id[|2\nr|] :< __tmo_input_id[|1\nr-1|])) {
        errprintf("tmo: TMO input rows must be sorted by location.\n")
        _error(3498)
    }

    info = panelsetup(__tmo_input_id, 1)
    nb = rows(info)
    __tmo_locid = __tmo_input_id[info[,1]]
    panel = (nb < rows(__tmo_input_id))

    if (panel) {
        // w_i^(j) = sum_t a_it xtilde_it eps_it^(j).
        // Screening uses these location-level weighted score sums.
        if (hasmissing(__tmo_input_R))
            M = panelsum(editmissing(__tmo_input_R,0), __tmo_input_w, info)
        else M = panelsum(__tmo_input_R, __tmo_input_w, info)
        R = M
        if (hasmissing(__tmo_input_R)) R=R:/(panelsum(__tmo_input_R:<.,info):>0)
        __tmo_W = M
    }
    else {
        R = __tmo_input_R
        xt = __tmo_input_w
        __tmo_W = xt :* R
    }

    // The outcome of interest must be complete. An auxiliary missing for
    // more than misslimit of the sample is dropped; the rest are kept and
    // each pair of locations is compared over the outcomes both observe.
    if (hasmissing(R[,1])) {
        errprintf("tmo: the outcome of interest has a missing residual.\n")
        _error(3498)
    }

    // Auxiliaries with no variation left after the regression are dropped:
    // flagged by the caller at the residual level, and here, in panels, when
    // the location score sums cancel to rounding noise.
    novar = J(1,cols(R),0)
    if (args()>=3) {
        if (novarname!="") {
            fl = st_matrix(novarname)
            if (rows(fl)!=1 | cols(fl)!=cols(R)) {
                errprintf("tmo: novar() must be a 1 x %g matrix.\n", cols(R))
                _error(3498)
            }
            novar = (fl:!=0)
        }
    }
    if (panel) {
        rms = sqrt(colsum((editmissing(__tmo_input_R,0):*__tmo_input_w):^2) :/
            colsum(__tmo_input_R:<.))
        for (j=2; j<=cols(R); j++) {
            if (novar[j] | !(rms[j]>0)) continue
            mj = select(R[,j], R[,j]:<.)
            if (rows(mj)<2) continue
            if (sqrt(quadvariance(mj)) <= 1e-8*rms[j]) {
                novar[j] = 1
                pnames = findexternal("__tmo_input_names")
                if (pnames!=NULL) {
                    if (cols(*pnames)==cols(R))
                        printf("{txt}tmo: %s has no variation in its location scores; dropped.\n", (*pnames)[j])
                }
            }
        }
    }
    novar[1] = 0
    miss = colsum(R:>=.):/rows(R)
    keepc = 1, (selectindex((miss[2..cols(R)]:<=misslim) :& !novar[2..cols(R)]):+1)
    __tmo_keepcols = keepc
    st_numscalar("__tmo_ndropped", sum((miss[2..cols(R)]:>misslim) :& !novar[2..cols(R)]))
    st_numscalar("__tmo_nnovar", sum(novar[2..cols(R)]))
    if (cols(keepc)<cols(R)) {
        R = R[,keepc]
        __tmo_W = __tmo_W[,keepc]
    }

    n = rows(R)
    if (n<=2 | cols(R)<3) {
        errprintf("tmo: TMO needs three or more locations, and two or\n")
        errprintf("     more auxiliary outcomes inside misslimit().\n")
        _error(3498)
    }
    st_numscalar("__tmo_nout", cols(R)-1)
    st_numscalar("__tmo_hasmiss", hasmissing(R))
    __tmo_W = editmissing(__tmo_W, 0)

    __tmo_R = R
    __tmo_Xs = R
    if (nrm!="") {
        nob = colsum(R:<.)
        sd = (colsum(editmissing((R :- (colsum(editmissing(R,0)):/nob)):^2,0)):/nob):^0.5
        if (min(sd)<=0 | hasmissing(sd)) {
            errprintf("tmo: an outcome has no residual variation.\n")
            _error(3498)
        }
        __tmo_Xs = R :/ sd
    }
}

void tmo_mask_prepare(real scalar n, real scalar masked)
{
    external real matrix __tmo_Bpairs, __tmo_Browinfo
    external real colvector __tmo_Bgroup
    real matrix info
    real colvector ix,sz
    real scalar q, m
    __tmo_Browinfo=J(n,2,0)
    if (!masked) {
        __tmo_Bpairs=J(0,2,.)
        __tmo_Bgroup=J(0,1,.)
    }
    if (findexternal("__tmo_Bgroup")==NULL) __tmo_Bgroup=J(0,1,.)
    if (cols(__tmo_Bpairs)!=2 | hasmissing(__tmo_Bpairs)) {
        errprintf("tmo: protected pairs must be a two-column index matrix.\n")
        _error(3498)
    }
    q=rows(__tmo_Bpairs)
    if (rows(__tmo_Bgroup)>0) {
        if (rows(__tmo_Bgroup)!=n | cols(__tmo_Bgroup)!=1 | q |
            hasmissing(__tmo_Bgroup) | min(__tmo_Bgroup)<1 |
            any(__tmo_Bgroup:!=floor(__tmo_Bgroup))) _error(3498)
        ix=order(__tmo_Bgroup,1); info=panelsetup(__tmo_Bgroup[ix],1)
        sz=info[,2]-info[,1]:+1
        if (max(__tmo_Bgroup)!=rows(info)) _error(3498)
        q=quadcross(sz,sz:-1)/2
    }
    if (rows(__tmo_Bpairs)>0) {
        if (min(__tmo_Bpairs)<1 | max(__tmo_Bpairs)>n |
            any(__tmo_Bpairs:!=floor(__tmo_Bpairs)) |
            any(__tmo_Bpairs[,1]:<=__tmo_Bpairs[,2])) {
            errprintf("tmo: protected pairs must have integer indices 1 <= j < i <= N_loc.\n")
            _error(3498)
        }
        if (q>1) {
            if (any((__tmo_Bpairs[|2,1\q,1|]:<__tmo_Bpairs[|1,1\q-1,1|]):|
                ((__tmo_Bpairs[|2,1\q,1|]:==__tmo_Bpairs[|1,1\q-1,1|]):&
                (__tmo_Bpairs[|2,2\q,2|]:<__tmo_Bpairs[|1,2\q-1,2|]))))
                __tmo_Bpairs=sort(__tmo_Bpairs,(1,2))
        }
        if (q>1) if (any(rowsum(__tmo_Bpairs[|2,1\q,2|]:==__tmo_Bpairs[|1,1\q-1,2|]):==2)) {
            errprintf("tmo: protected pairs contain duplicate entries.\n")
            _error(3498)
        }
        info=panelsetup(__tmo_Bpairs,1)
        __tmo_Browinfo[__tmo_Bpairs[info[,1],1],]=info
    }
    m=n*(n-1)/2-q
    st_numscalar("__tmo_nprotected",q)
    st_numscalar("__tmo_neligible",m)
    st_numscalar("__tmo_calibration_weight",masked ? 2*m : n*n)
    if (masked & m<=0) {
        errprintf("tmo: no unprotected off-diagonal pairs remain for calibration.\n")
        errprintf("     Use less extensive cluster or distance protection.\n")
        _error(459)
    }
}

// Lower-triangle protection indicator for row i; the diagonal is implicit.
real colvector tmo_protected_row(real scalar i)
{
    external real matrix __tmo_Bpairs, __tmo_Browinfo
    external real colvector __tmo_Bgroup
    real colvector b
    real scalar lo, hi
    if (i<1 | i>rows(__tmo_Browinfo) | i!=floor(i)) _error(3498)
    if (rows(__tmo_Bgroup)) return(__tmo_Bgroup[|1\i|]:==__tmo_Bgroup[i])
    b=J(i,1,0); b[i]=1
    lo=__tmo_Browinfo[i,1]; hi=__tmo_Browinfo[i,2]
    if (lo>0) b[__tmo_Bpairs[|lo,2\hi,2|]]=J(hi-lo+1,1,1)
    return(b)
}

// An empty protection set needs no square working matrix.
real matrix tmo_protected_matrix(real scalar n)
{
    external real matrix __tmo_Bpairs
    external real colvector __tmo_Bgroup
    real matrix B
    real scalar k, i, j
    if (rows(__tmo_Bgroup)) {
        B=J(n,n,0)
        for (i=1;i<=n;i++) B[i,]=__tmo_Bgroup':==__tmo_Bgroup[i]
        return(B)
    }
    if (rows(__tmo_Bpairs)==0) return(J(0,0,.))
    B=J(n,n,0)
    for (k=1;k<=rows(__tmo_Bpairs);k++) {
        i=__tmo_Bpairs[k,1]; j=__tmo_Bpairs[k,2]
        B[i,j]=B[j,i]=1
    }
    return(B)
}

// Three-way order-statistic selection; a bounded fallback handles bad pivots.
real scalar tmo_kth(real colvector v, real scalar k)
{
    real scalar lo, hi, mid, pivot, left, right, i, z, steps, budget
    real colvector part
    lo=1; hi=rows(v); steps=0
    if (k<1 | k>hi | k!=floor(k)) _error(3498)
    budget=2*ceil(ln(hi+1)/ln(2))
    while (lo<hi) {
        steps++
        if (steps>budget) {
            part=sort(v[|lo\hi|],1)
            return(part[k-lo+1])
        }
        mid=floor((lo+hi)/2)
        pivot=v[mid]
        if (v[lo]>v[hi]) {
            if (pivot<v[hi]) pivot=v[hi]
            else if (pivot>v[lo]) pivot=v[lo]
        }
        else {
            if (pivot<v[lo]) pivot=v[lo]
            else if (pivot>v[hi]) pivot=v[hi]
        }
        left=lo; right=hi; i=lo
        while (i<=right) {
            if (v[i]<pivot) {
                z=v[left]; v[left]=v[i]; v[i]=z
                left++; i++
            }
            else if (v[i]>pivot) {
                z=v[right]; v[right]=v[i]; v[i]=z
                right--
            }
            else i++
        }
        if (k<left) hi=left-1
        else if (k>right) lo=right+1
        else return(pivot)
    }
    return(v[k])
}

// Inverse empirical CDF over all n^2 pairs, the diagonal included.
real rowvector tmo_iqr_all(real matrix C)
{
    real colvector v
    real scalar m, r25, r75, q25, q75
    if (rows(C)!=cols(C) | hasmissing(C)) _error(3498)
    v = vec(C)
    m = rows(v)
    r25 = ceil(0.25*m); r75 = ceil(0.75*m)
    if (r25<1) r25 = 1
    if (r75<1) r75 = 1
    q25 = tmo_kth(v,r25)
    q75 = tmo_kth(v,r75)
    return((q25,q75))
}

// Signed eligible pairs only; the two orientations have equal multiplicity.
real rowvector tmo_iqr_eligible(real matrix C)
{
    real colvector v, b, rowv
    real scalar n, m, i, k, nr, q25, q75
    n=rows(C)
    if (cols(C)!=n | hasmissing(C)) _error(3498)
    m=st_numscalar("__tmo_neligible")
    if (m<=0 | missing(m)) _error(459)
    v=J(m,1,.); k=1
    for (i=2;i<=n;i++) {
        b=tmo_protected_row(i)
        rowv=select(C[|i,1\i,i-1|]',b[|1\i-1|]:==0)
        nr=rows(rowv)
        if (nr>0) {
            v[|k\k+nr-1|]=rowv
            k=k+nr
        }
    }
    if (k!=m+1 | hasmissing(v)) _error(3498)
    q25=tmo_kth(v,ceil(0.25*m))
    q75=tmo_kth(v,ceil(0.75*m))
    return((q25,q75))
}

// Enforce the numerical range of a correlation without a tolerance rule.
real matrix tmo_bound_corr(real matrix C)
{
    real matrix A
    real scalar i
    A=abs(C)
    A=C:/(A:*(A:>1)+(A:<=1))
    for (i=1;i<=rows(A);i++) A[i,i]=1
    return(A)
}

// The screen for one anchor, formed explicitly. Each auxiliary is
// orthogonalized on the anchor over the locations where both are observed,
// and each pair of locations is normalized over the outcomes both observe.
real scalar tmo_pair_cosine(real rowvector x, real rowvector y, real rowvector support)
{
    real rowvector u,v
    real scalar sx,sy
    u=select(x,support); v=select(y,support)
    if (!length(u)) return(.)
    sx=max(abs(u)); sy=max(abs(v))
    if (sx<=0 | sy<=0 | missing(sx) | missing(sy)) return(.)
    if (length(u)==1) return(sign(u[1])*sign(v[1]))
    u=u:/sx; v=v:/sy
    if (all(u:==v)) return(1)
    if (all(u:==-v)) return(-1)
    return(quadcross(u',v')/sqrt(quadcross(u',u'))/sqrt(quadcross(v',v')))
}

// A projection with no shared anchor variation is omitted for this anchor.
// Exceptional arithmetic is evaluated on the shared support in scaled units.
void tmo_project(real matrix Xs, real scalar a, real matrix En, real matrix Mb)
{
    real matrix E
    real colvector x,ix,u,v,r
    real rowvector den,num,bt
    real scalar b,sx,sy,gamma,minnormal,scaled,overflow
    E=a==1 ? Xs[,2..cols(Xs)] : (a==cols(Xs) ? Xs[,1..a-1] : (Xs[,1..a-1],Xs[,a+1..cols(Xs)]))
    x=editmissing(Xs[,a],0); Mb=(E:<.):*(Xs[,a]:<.)
    En=editmissing(E,0)
    r=x:^2; overflow=hasmissing(r)
    den=quadcross(r,Mb); num=quadcross(x,En); bt=num:/den
    minnormal=2.2250738585072014e-308
    for (b=1;b<=cols(E);b++) {
        ix=vec(selectindex(Mb[,b]))
        if (!rows(ix)) {
            En[,b]=J(rows(x),1,0)
            continue
        }
        r=En[ix,b]-x[ix]*bt[b]
        scaled=0
        if (num[b]==0) {
            sx=max(abs(x[ix])); sy=max(abs(En[ix,b]))
            scaled=sx>0 & sy>0 & (sx*sy<minnormal | missing(sx*sy))
        }
        if (overflow | missing(bt[b]) | den[b]<minnormal | missing(den[b]) |
            (num[b]!=0 & abs(num[b])<minnormal) | scaled | hasmissing(r)) {
            u=x[ix]; v=En[ix,b]; sx=max(abs(u)); sy=max(abs(v))
            if (sx==0) {
                Mb[,b]=J(rows(x),1,0); En[,b]=J(rows(x),1,0)
                continue
            }
            if (sy==0) r=J(rows(ix),1,0)
            else {
                u=u:/sx; v=v:/sy
                if (any((x[ix]:!=0):&(abs(u):<minnormal)) |
                    any((En[ix,b]:!=0):&(abs(v):<minnormal)) |
                    any((u:!=0):&(v:!=0):&(abs(u:*v):<minnormal))) {
                    errprintf("tmo: auxiliary projection loses precision from extreme within-outcome dynamic range.\n")
                    errprintf("     Inspect the observations or transform the affected outcome.\n")
                    _error(459)
                }
                gamma=quadcross(u,v)/quadcross(u,u)
                r=(v-u*gamma)*sy
            }
            if (hasmissing(r)) {
                errprintf("tmo: auxiliary projection overflow; rescale outcome units.\n")
                _error(459)
            }
        }
        En[,b]=J(rows(x),1,0); En[ix,b]=r
    }
}

real matrix tmo_screen_a(real matrix Xs, real scalar a, | real scalar count_support)
{
    real matrix Mb, En, num, SS, C
    real scalar p, i, j, minnormal
    real colvector nobs, single
    p = cols(Xs)
    if (args()<3) count_support=0
    tmo_project(Xs,a,En,Mb)
    num = En*En'
    SS = (En:^2)*Mb'
    C = num :/ sqrt(SS) :/ sqrt(SS')
    minnormal=2.2250738585072014e-308
    if (hasmissing(C) | min(SS)<minnormal |
        any((abs(num):>0):&(abs(num):<minnormal)) | any((num:!=0):&(C:==0))) {
        for (i=1;i<=rows(C);i++) for (j=1;j<i;j++) {
            if (missing(C[i,j]) | SS[i,j]<minnormal | SS[j,i]<minnormal |
                (num[i,j]!=0 & (abs(num[i,j])<minnormal | C[i,j]==0)))
                C[i,j]=C[j,i]=tmo_pair_cosine(En[i,],En[j,],Mb[i,]:&Mb[j,])
        }
    }
    // Recompute near-boundary arithmetic before exact proportionality checks.
    for (i=2;i<=rows(C);i++) {
        single=vec(selectindex((abs(C[|i,1\i,i-1|]'):>=1-sqrt(epsilon(1))):&
            (C[|i,1\i,i-1|]':<.)))
        for (j=1;j<=rows(single);j++)
            C[i,single[j]]=C[single[j],i]=tmo_pair_cosine(En[i,],En[single[j],],Mb[i,]:&Mb[single[j],])
    }
    // A valid one-dimensional cosine is exactly its product's sign.
    nobs=rowsum(Mb)
    nobs=select(nobs,nobs:>0)
    if (rows(nobs)>0) {
        if (2*min(nobs)<=cols(Mb)+1) {
            SS=Mb*Mb'
            for (i=1;i<=rows(C);i++) {
                single=selectindex((SS[,i]:==1):&(C[,i]:<.))
                if (rows(single)>0) C[single,i]=sign(C[single,i])
            }
        }
    }
    for (i=1; i<=rows(C); i++) {
        C[i,i]=1
        if (i>1) C[|1,i\i-1,i|]=C[|i,1\i,i-1|]'
    }
    if (count_support) st_numscalar("__tmo_unsupported_pairs",sum(C:>=.)/2)
    return(editmissing(tmo_bound_corr(C),0))
}

real rowvector tmo_dd_add(real rowvector acc, real scalar x)
{
    real scalar s,e,t,h
    s=acc[1]+x
    e=abs(acc[1])>=abs(x) ? (acc[1]-s)+x : (x-s)+acc[1]
    t=acc[2]+e; h=s+t
    return((h,(abs(s)>=abs(t) ? (s-h)+t : (t-h)+s)))
}
real rowvector tmo_dd_product(real rowvector acc, real scalar x, real scalar y, real scalar factor)
{
    real scalar p,tx,ty,hx,hy,lx,ly,e
    p=x*y; tx=134217729*x; ty=134217729*y
    hx=tx-(tx-x); hy=ty-(ty-y); lx=x-hx; ly=y-hy
    e=((hx*hy-p)+hx*ly+lx*hy)+lx*ly
    return(tmo_dd_add(tmo_dd_add(acc,factor*p),factor*e))
}
// Recover the low word of a quad-accumulated sum without a scalar pair loop.
real rowvector tmo_dd_sum(real colvector x)
{
    real scalar s,r,i
    real rowvector d
    real colvector nz
    if (!rows(x)) return((0,0))
    if (rows(x)==1) return((x[1],0))
    nz=select(abs(x),x:!=0)
    if (length(nz)) {
        // Extreme exponent spreads can exceed even a quad accumulator's span.
        if (max(nz)>min(nz)*2^80) {
            d=(0,0)
            for (i=1;i<=rows(x);i++) d=tmo_dd_add(d,x[i])
            return(d)
        }
    }
    s=quadcross(J(rows(x),1,1),x)
    d=tmo_dd_add((x[1],0),-s)
    r=quadcross(J(rows(x),1,1),(d[1]\x[|2\rows(x)|]))
    return(tmo_dd_add(tmo_dd_add((s,0),r),d[2]))
}

// Vectorize products over a bounded pair block, retaining product roundoff.
void tmo_accumulate_pairs(real matrix buckets, real rowvector countpath,
    real colvector left, real colvector right, real colvector score, real colvector dgrid)
{
    real colvector p,tx,ty,hx,hy,lx,ly,e,bins,ix
    real rowvector s,t
    real scalar g
    p=left:*right; tx=134217729*left; ty=134217729*right
    hx=tx-(tx-left); hy=ty-(ty-right); lx=left-hx; ly=right-hy
    e=((hx:*hy-p)+hx:*ly+lx:*hy)+lx:*ly
    bins=J(rows(score),1,0)
    for (g=1;g<=rows(dgrid);g++) bins=bins+(score:>=dgrid[g])
    for (g=1;g<=rows(dgrid);g++) {
        ix=vec(selectindex(bins:==g)); if (!rows(ix)) continue
        s=tmo_dd_sum(p[ix]); t=tmo_dd_sum(e[ix])
        buckets[g,]=tmo_dd_add(tmo_dd_add(buckets[g,],2*s[1]),2*s[2])
        buckets[g,]=tmo_dd_add(tmo_dd_add(buckets[g,],2*t[1]),2*t[2])
        countpath[g]=countpath[g]+rows(ix)
    }
}

// B is either a private dense mask or a sorted lower-triangle pair list.
real colvector tmo_mask_ptr(real matrix B, real scalar n)
{
    real colvector ptr
    real scalar i,k
    ptr=J(n+1,1,1); k=1
    if (!rows(B) | cols(B)!=2) return(ptr)
    for (i=1;i<=n;i++) {
        ptr[i]=k
        while (k<=rows(B)) {
            if (B[k,1]!=i) break
            if (B[k,2]<1 | B[k,2]>=i) _error(3498)
            if (k>ptr[i]) {
                if (B[k,2]<=B[k-1,2]) _error(3498)
            }
            k++
        }
    }
    ptr[n+1]=k
    if (k!=rows(B)+1) _error(3498)
    return(ptr)
}

real rowvector tmo_baseline_dd(real colvector wa, real matrix B, real colvector groups,
                               real colvector ptr)
{
    real colvector ix,jj
    real matrix info
    real rowvector acc,s
    real scalar i
    acc=(0,0)
    if (rows(groups)) {
        ix=order(groups,1); info=panelsetup(groups[ix],1)
        for (i=1;i<=rows(info);i++) {
            s=tmo_dd_sum(wa[ix[|info[i,1]\info[i,2]|]])
            acc=tmo_dd_product(acc,s[1],s[1],1)
            acc=tmo_dd_product(acc,s[1],s[2],2)
            acc=tmo_dd_product(acc,s[2],s[2],1)
        }
        return(acc)
    }
    for (i=1;i<=rows(wa);i++) {
        acc=tmo_dd_product(acc,wa[i],wa[i],1)
        if (i==1 | !rows(B)) continue
        if (cols(B)==2) jj=ptr[i]<ptr[i+1] ? B[|ptr[i],2\ptr[i+1]-1,2|] : J(0,1,.)
        else jj=vec(selectindex(B[|i,1\i,i-1|]'))
        s=tmo_dd_sum(wa[jj])
        acc=tmo_dd_product(acc,wa[i],s[1],2)
        acc=tmo_dd_product(acc,wa[i],s[2],2)
    }
    return(acc)
}

real scalar tmo_baseline(real colvector wa, real matrix B, real colvector groups)
{
    real rowvector acc
    acc=tmo_baseline_dd(wa,B,groups,tmo_mask_ptr(B,rows(wa)))
    return(acc[1]+acc[2])
}

// Keep high/low words through BOTH protected and screened contributions.
// Bounded pair buckets avoid an n-by-n mask/matrix product at each candidate.
real rowvector tmo_variance_path(real matrix A, real colvector wa,
                                  real colvector dgrid, | real matrix B, real colvector groups, real scalar counts)
{
    real scalar n,G,g,i,scale,power,capacity,used,first,take
    real matrix buckets
    real colvector ptr,eligible,jj,u,left,right,score,part
    real rowvector out,countpath,acc
    n=rows(wa); G=rows(dgrid)
    if (G<1 | (min(dgrid)<=1 & (rows(A)!=n | cols(A)!=n))) _error(3498)
    if (args()<4) B=J(0,0,.)
    if (args()<5) groups=J(0,1,.)
    if (args()<6) counts=0
    if (rows(B)>0 & cols(B)!=2 & (rows(B)!=n | cols(B)!=n)) _error(3498)
    if (rows(groups)>0 & rows(groups)!=n) _error(3498)
    if (G>1) {
        if (any(dgrid[|2\G|]:<dgrid[|1\G-1|])) _error(3498)
    }
    ptr=tmo_mask_ptr(B,n)
    // Binary scaling preserves score bits and protects Dekker splitting.
    scale=1; power=max(abs(wa))
    if (power>1e100 | (power>0 & power<1e-100)) {
        power=max((-1022,min((1023,floor(log(power)/log(2))))))
        scale=2^power
    }
    u=wa:/scale
    acc=tmo_baseline_dd(u,B,groups,ptr)
    out=J(1,G,.); countpath=J(1,G,0); buckets=J(G,2,0)
    capacity=min((65536,n*(n-1)/2)); used=0
    if (min(dgrid)<=1) left=right=score=J(capacity,1,.)
    if (min(dgrid)<=1) for (i=2;i<=n;i++) {
        eligible=J(i-1,1,1)
        if (rows(groups)) eligible=groups[|1\i-1|]:!=groups[i]
        else if (rows(B)) {
            if (cols(B)==2) {
                if (ptr[i]<ptr[i+1]) eligible[B[|ptr[i],2\ptr[i+1]-1,2|]]=J(ptr[i+1]-ptr[i],1,0)
            }
            else eligible=B[|i,1\i,i-1|]':==0
        }
        eligible=eligible:&(A[|i,1\i,i-1|]':>=dgrid[1])
        jj=vec(selectindex(eligible)); if (!rows(jj)) continue
        first=1
        while (first<=rows(jj)) {
            take=min((capacity-used,rows(jj)-first+1))
            part=jj[|first\first+take-1|]
            left[|used+1\used+take|]=J(take,1,u[i])
            right[|used+1\used+take|]=u[part]
            score[|used+1\used+take|]=A[i,part]'
            used=used+take; first=first+take
            if (used==capacity) {
                tmo_accumulate_pairs(buckets,countpath,left,right,score,dgrid)
                used=0
            }
        }
    }
    if (used) tmo_accumulate_pairs(buckets,countpath,left[|1\used|],right[|1\used|],score[|1\used|],dgrid)
    for (g=G;g>=1;g--) {
        acc=tmo_dd_add(tmo_dd_add(acc,buckets[g,1]),buckets[g,2])
        out[g]=((acc[1]+acc[2])*scale)*scale
        if (g<G) countpath[g]=countpath[g]+countpath[g+1]
    }
    if (hasmissing(out)) {
        errprintf("tmo: variance contributions overflow; rescale outcome or treatment units.\n")
        _error(459)
    }
    if (counts) st_matrix("__tmo_screened_counts",countpath)
    return(out)
}

// Use rank-two screens only when every projected row norm is well conditioned.
// This guard changes arithmetic, never the literal threshold comparison.
real scalar tmo_complete_setup(real matrix Xs, real matrix Z, real colvector gd, real matrix GM)
{
    real matrix gram,H,sub,v
    real colvector ss,hh
    gram=quadcross(Xs,Xs)
    ss=diagonal(gram)
    if (min(ss)<=0 | hasmissing(ss)) return(1)
    H=gram:/ss'
    hh=colsum(H:^2)'
    Z=Xs*H-Xs:*(hh'/2)
    gd=rowsum(Xs:^2)
    sub=2*Z:*Xs; v=gd:-sub
    if (hasmissing(v) | any(v:<=sqrt(epsilon(1))*(abs(gd):+abs(sub))) |
        any(v:/(cols(Xs)-1):<=0)) return(1)
    GM=Xs*Xs'
    return(0)
}

real matrix tmo_complete_screen(real matrix Xs, real scalar a, real matrix Z,
                                 real colvector gd, real matrix GM)
{
    real matrix S,L,R,En,Mb,bad
    real colvector xa,rsd
    real scalar i,j,d
    d=cols(Xs)-1; xa=Xs[,a]
    L=Z[,a]*xa'; R=xa*Z[,a]'
    S=GM-L-R
    bad=(abs(S):<=sqrt(epsilon(1))*(abs(GM)+abs(L)+abs(R))):|(S:>=.)
    rsd=1:/sqrt((gd-2*(Z[,a]:*xa)):/d)
    S=((S:/d):*rsd):*rsd'
    bad=bad:|(abs(S):>=1-sqrt(epsilon(1)))
    for (i=1;i<=rows(S);i++) bad[i,i]=0
    if (any(bad)) {
        tmo_project(Xs,a,En,Mb)
        for (i=1;i<=rows(S);i++) for (j=1;j<i;j++) {
            if (bad[i,j]) S[i,j]=S[j,i]=tmo_pair_cosine(En[i,],En[j,],Mb[i,]:&Mb[j,])
        }
    }
    for (i=1;i<=rows(S);i++) {
        S[i,i]=1
        if (i>1) S[|1,i\i-1,i|]=S[|i,1\i,i-1|]'
    }
    return(editmissing(tmo_bound_corr(S),0))
}

// Calibration, candidate paths, and pair exports share one screen evaluation.
real matrix tmo_target_screen(real matrix Xs, real scalar hasmiss)
{
    real matrix Z,GM
    real colvector gd
    if (hasmiss) return(tmo_screen_a(Xs,1))
    if (tmo_complete_setup(Xs,Z,gd,GM)) return(tmo_screen_a(Xs,1))
    return(tmo_complete_screen(Xs,1,Z,gd,GM))
}

// One variance path per anchor. With a common outcome set the p anchors
// share one Gram matrix and each costs two rank-one updates; otherwise the
// screen is formed per anchor so that the shared-outcome normalization is
// exact.
real matrix tmo_anchor_paths_miss(real matrix Xs, real matrix W,
                                   real colvector dgrid, | real matrix B, real colvector groups, real scalar counts,
                                   real rowvector targetpath)
{
    real scalar p, G, a
    real matrix A, OUT
    p = cols(Xs); G = rows(dgrid)
    if (args()<4) B=J(0,0,.)
    if (args()<5) groups=J(0,1,.)
    if (args()<6) counts=0
    if (args()<7) targetpath=J(1,0,.)
    if (cols(targetpath)>0 & cols(targetpath)!=G) _error(3498)
    OUT = J(p, G, .)
    if (cols(targetpath)) OUT[1,]=targetpath
    for (a=1; a<=p; a++) {
        if (a==1 & cols(targetpath)) continue
        A = min(dgrid)>1 ? J(0,0,.) : abs(tmo_screen_a(Xs, a))
        OUT[a,] = tmo_variance_path(A,W[,a],dgrid,B,groups,counts & a==1)
        A = J(0,0,.)
    }
    return(OUT)
}

real matrix tmo_anchor_paths(real matrix Xs, real matrix W,
                              real colvector dgrid, | real matrix B, real colvector groups, real scalar counts,
                              real rowvector targetpath, real matrix Z, real colvector gd, real matrix GM)
{
    real scalar p,G,a
    real matrix A,OUT
    p=cols(Xs); G=rows(dgrid)
    if (args()<4) B=J(0,0,.)
    if (args()<5) groups=J(0,1,.)
    if (args()<6) counts=0
    if (args()<7) targetpath=J(1,0,.)
    if (cols(targetpath)>0 & cols(targetpath)!=G) _error(3498)
    if (min(dgrid)>1) return(tmo_anchor_paths_miss(Xs,W,dgrid,B,groups,counts,targetpath))
    if (args()<10) {
        if (tmo_complete_setup(Xs,Z,gd,GM)) return(tmo_anchor_paths_miss(Xs,W,dgrid,B,groups,counts,targetpath))
    }
    else if (rows(Z)!=rows(Xs) | cols(Z)!=p | rows(gd)!=rows(Xs) |
        rows(GM)!=rows(Xs) | cols(GM)!=rows(Xs)) _error(3498)
    OUT = J(p, G, .)
    if (cols(targetpath)) OUT[1,]=targetpath
    for (a=1; a<=p; a++) {
        if (a==1 & cols(targetpath)) continue
        A=abs(tmo_complete_screen(Xs,a,Z,gd,GM))
        OUT[a,] = tmo_variance_path(A,W[,a],dgrid,B,groups,counts & a==1)
        A = J(0,0,.)
    }
    return(OUT)
}

real matrix tmo_anchor_stats(real matrix R, real rowvector sourcecols)
{
    real scalar p
    real colvector ss, om, den, ph
    real matrix R0, out
    p = cols(R)
    if (p<3 | cols(sourcecols)!=p) _error(3498)
    R0 = editmissing(R,0)
    ss = colsum(R0:^2)'
    if (min(ss)<=0 | hasmissing(ss) | missing(max(abs(R0))^2)) {
        errprintf("tmo: an outcome has zero or unrepresentable residual sum of squares.\n")
        errprintf("     Remove constant auxiliary outcomes or rescale extreme outcome units.\n")
        _error(3498)
    }
    om = quadcross(R0[,1]:/sqrt(ss[1]), R0[|1,2 \ rows(R0),p|]:/sqrt(ss[|2\p|])')'
    om = rowmax((rowmin((om, J(p-1,1,1))), J(p-1,1,-1)))
    den = 1 :- om:^2
    if (min(den) <= 1e-12 | hasmissing(den)) {
        errprintf("tmo: an auxiliary outcome is perfectly or nearly collinear with the target.\n")
        errprintf("     Remove the collinear outcome from ylist(); its penalty denominator is not usable.\n")
        _error(3498)
    }
    ph = ss[1] :/ ss[|2 \ p|]
    if (min(ph)<=0 | hasmissing(ph)) {
        errprintf("tmo: an auxiliary variance-rescaling ratio is outside the numerical range.\n")
        errprintf("     Rescale extreme outcome units before calling tmo.\n")
        _error(3498)
    }
    out=sourcecols', ss, (1\ph), (1\om), (0\den)
    return(out)
}

real rowvector tmo_xi_stats(real matrix paths, real matrix stats)
{
    real scalar p
    real rowvector v0,xi,scale
    real matrix D,D2
    p=rows(paths)
    if (rows(stats)!=p | cols(stats)!=5 | p<3) _error(3498)
    v0 = paths[1,]
    D = stats[|2,3\p,3|]:*paths[|2,1\p,cols(paths)|] :- v0
    if (hasmissing(D) | hasmissing(stats[,5]) | min(stats[|2,5\p,5|])<=0) {
        errprintf("tmo: the auxiliary variance differences or penalty denominators are not representable.\n")
        errprintf("     Rescale extreme outcome units or inspect collinear auxiliary outcomes.\n")
        _error(459)
    }
    D2 = D:^2 :/ stats[|2,5\p,5|]
    xi = sqrt(colsum(D2):/(p-1))
    // Scale exceptional disagreements before squaring; preserve ordinary arithmetic.
    if (hasmissing(D2) | hasmissing(xi) |
        any((abs(D):>0):&(abs(D):<sqrt(2.2250738585072014e-308)))) {
        scale=colmax(abs(D)); scale=scale:+(scale:==0)
        D2=(D:/scale):^2:/stats[|2,5\p,5|]
        if (hasmissing(D2)) _error(459)
        xi=scale:*sqrt(colsum(D2):/(p-1))
    }
    if (hasmissing(xi)) {
        errprintf("tmo: the selection penalty exceeds the numerical range; rescale outcome units.\n")
        _error(459)
    }
    return(xi)
}

real rowvector tmo_xi(real matrix paths, real matrix R)
{
    return(tmo_xi_stats(paths,tmo_anchor_stats(R,1..cols(R))))
}

// Reference path for when the compiled kernel is unavailable. It also
// supplies v0hat, which the kernel computes for itself.
void tmo_paths_mata(| real scalar masked)
{
    external real matrix __tmo_Xs, __tmo_W, __tmo_Bpairs
    external real colvector __tmo_dgrid, __tmo_etagrid, __tmo_Bgroup
    real matrix C,B,Z,GM
    real colvector sdv,groups
    real rowvector q,targetpath
    real scalar p, v0, lf, counts
    if (args()<1) masked=0
    p = cols(__tmo_Xs)
    counts=st_numscalar("__tmo_need_counts")==1
    st_numscalar("__tmo_unsupported_pairs",0)
    if (st_numscalar("__tmo_hasmiss")) st_numscalar("__tmo_explicit_screen",1)
    else st_numscalar("__tmo_explicit_screen",tmo_complete_setup(__tmo_Xs,Z,sdv,GM))
    if (st_numscalar("__tmo_explicit_screen")) C=tmo_screen_a(__tmo_Xs,1,counts)
    else C=tmo_complete_screen(__tmo_Xs,1,Z,sdv,GM)
    q = masked ? tmo_iqr_eligible(C) : tmo_iqr_all(C)
    st_numscalar("__tmo_q25",q[1])
    st_numscalar("__tmo_q75",q[2])
    v0 = ((q[2]-q[1])/(2*invnormal(0.75)))^2
    if (v0<=0 | missing(v0)) {
        if (masked) {
            errprintf("tmo: the eligible-pair screening IQR is not positive.\n")
            errprintf("     Use less extensive protection or auxiliary outcomes with more screening variation.\n")
            _error(459)
        }
        else {
            errprintf("tmo: the interquartile range of the screening statistic is zero.\n")
            _error(3498)
        }
    }
    lf = st_numscalar("__tmo_logfactor")
    __tmo_dgrid = sqrt(2:*__tmo_etagrid:*v0:*lf)
    st_numscalar("__tmo_v0hat", v0)
    B=masked ? __tmo_Bpairs : J(0,0,.)
    groups=masked ? __tmo_Bgroup : J(0,1,.)
    targetpath=tmo_variance_path(min(__tmo_dgrid)>1 ? J(0,0,.) : abs(C),__tmo_W[,1],__tmo_dgrid,B,groups,counts)
    C=J(0,0,.)
    if (st_numscalar("__tmo_explicit_screen")) {
        Z=GM=J(0,0,.)
        st_matrix("__tmo_paths", tmo_anchor_paths_miss(__tmo_Xs, __tmo_W, __tmo_dgrid,B,groups,counts,targetpath))
    }
    else
        st_matrix("__tmo_paths", tmo_anchor_paths(__tmo_Xs, __tmo_W, __tmo_dgrid,B,groups,counts,targetpath,Z,sdv,GM))
    if (counts) st_numscalar("__tmo_counts_api",1)
}

// Objectives within 0.01% of the maximum are ties; the largest tied
// admissible index wins.
real scalar tmo_tie_index(real rowvector obj, real rowvector admissible)
{
    real scalar b
    b = max(obj[selectindex(admissible)])
    return(max(selectindex(admissible :& (obj:>=b-1e-4*abs(b)))))
}

void tmo_select(real scalar usejava)
{
    external real matrix __tmo_R
    external real colvector __tmo_dgrid, __tmo_etagrid
    external real rowvector __tmo_keepcols
    real matrix paths, stats
    real rowvector v0, xi, obj, counts, eligiblepct, admissible
    real scalar xx2, G, i, v, cap, neligible
    if (usejava) __tmo_dgrid = st_matrix("__tmo_dgrid_out")
    xx2 = st_numscalar("__tmo_xx")^2
    if (xx2<=0 | missing(xx2)) _error(3498)
    paths = st_matrix("__tmo_paths") :/ xx2
    G = cols(paths)
    if (G!=rows(__tmo_dgrid) | G!=rows(__tmo_etagrid) |
        rows(paths)!=cols(__tmo_R) | hasmissing(paths) | hasmissing(__tmo_dgrid)) _error(3498)
    stats = tmo_anchor_stats(__tmo_R,__tmo_keepcols)
    v0 = paths[1,]
    xi = tmo_xi_stats(paths, stats)
    obj = v0 - xi
    if (hasmissing(obj)) _error(3498)
    cap = st_numscalar("__tmo_max_eligible_pct")
    if (cap<0 | cap>100 | missing(cap)) _error(3498)
    if (st_numscalar("__tmo_need_counts")==1) {
        counts = st_matrix("__tmo_screened_counts")
        neligible = st_numscalar("__tmo_neligible")
        if (rows(counts)!=1 | cols(counts)!=G | neligible<=0 |
                missing(neligible) | hasmissing(counts) |
                any(counts:<0) | any(counts:>neligible)) _error(3498)
        eligiblepct = 100:*counts:/neligible
        // Cross multiplication makes an exact inclusive comparison for the
        // integer percentage caps.
        admissible = 100:*counts :<= cap*neligible
        st_matrix("__tmo_eligible_pctpath",eligiblepct)
        st_matrix("__tmo_admissiblepath",admissible)
    }
    else admissible = J(1,G,1)
    if (cap<100 & !any(admissible)) {
        errprintf("tmo: no eta candidate satisfies maxeligiblepct(%g).\n",cap)
        errprintf("     Raise maxeligiblepct(); maxeligiblepct(100) removes the restriction.\n")
        _error(459)
    }
    i = tmo_tie_index(obj, admissible)
    v = v0[i]
    st_numscalar("__tmo_var", v)
    st_numscalar("__tmo_se", v>0 ? sqrt(v) : .)
    st_numscalar("__tmo_thres", __tmo_dgrid[i])
    st_numscalar("__tmo_eta", __tmo_etagrid[i])
    st_numscalar("__tmo_gridmax", __tmo_dgrid[G])
    st_numscalar("__tmo_interior", i>1 & i<G)
    st_numscalar("__tmo_selected_index",i)
    st_matrix("__tmo_vpath", v0)
    st_matrix("__tmo_xipath", xi)
    st_matrix("__tmo_dgrid_m", __tmo_dgrid')
    st_matrix("__tmo_etagrid_m", __tmo_etagrid')
    st_matrix("__tmo_anchorpath",paths)
    st_matrix("__tmo_anchorstats",stats)
    st_matrix("__tmo_objectivepath",obj)
}
end
