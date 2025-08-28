
# since GGG2020 provides the integration operator for converting wet mole fraction to column averages, let's simply use int_operator 

get.wgt.tccon.func = function(output, tccon.fn, tccon.species) {

    # grab receptor info and select the particles at first time step back
    receptor = output$receptor
    p = output$particle
	min.time = min(abs(p$time)) * sign(p$time[1])  # MIN time in mins
	sel.p = p[p$time == min.time, ] %>% 
            dplyr::select(indx, zagl, zsfc, pres, xhgt)

	# get TCCON profile first according to lat/lon of receptor, return a list
	if (is.na(tccon.fn)) stop('get.wgt.tccon.func(): missing TCCON file...\n')
	tccon.info = get.tccon.info(tccon.fn, receptor, tccon.species)

    if (is.null(tccon.info)) {
		warning('get.wgt.tccon.func(): NO TCCON info found for this receptor...please check\n')
		return()
	} # end if is.null


    # grab press weighting function, pressures, normalized AK, apriori
	# since aprior profiles are reported in WET mole fraction, conver to dry 
	tccon.df = tccon.info$wgt %>% 
			   mutate(ap_h2o_dry = ap_h2o_wet / (1 - ap_h2o_wet), 
			   		  ap_gas_dry = ap_gas_wet * (1 + ap_h2o_dry))

	if (F) {	# sanity check on wet-dry conversion
		t1 = ggplot(data = tccon.df) + labs(x = 'CO2 prior') + 
			 geom_point(aes(ap_gas_wet, ak_pres, color = 'wet MF')) + 
			 geom_point(aes(ap_gas_dry, ak_pres, color = 'dry MF'))
	}

	tccon.pres.bound = sort(as.numeric(tccon.df$prior_pres))
	tccon.ak.norm = as.numeric(tccon.df$ak_norm)
	tccon.zsfc = tccon.info$zasl_m				# alt above sea level in m
	tccon.psfc = tccon.info$psfc_hpa		    # in hPa
	tccon.xh2o = tccon.info$xh2o_ppm			# in ppm 
	tccon.xgas = tccon.info$x_gas 		

	# --------- correct for model-obs mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RTv_mean * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# Tv_mean is the mean virtual temp 

	# correct for diff in ZSFC between sensor and particle ZSFC
	p.zasl = sel.p$zagl + sel.p$zsfc 	# trajec-level ASL
	p.pres = sel.p$pres 				# trajec-level pressure
	p.zagl.corr.co2 = p.zasl - tccon.zsfc	# corrected trajec-level AGL

	# use satellite surface altitude and pressure to calculate ZAGL that TROPOMI believes
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ tccon.psfc * exp(a * (-p.zagl.corr.co2)), 
					 start = list(a = 1E-4))		
	a.tccon = coef(nls)[[1]]    # a stands for g/RTv_mean based on obs

	# use TCCON sfc pressure to calculate pressure of release heights
	# P = Psfc * exp (g/RT * (Zsfc - Z))
	sel.p$xpres = tccon.psfc * exp(a.tccon * (-sel.p$xhgt))


	# -------------------------------------------------------------------------
	# let's wet STILT footprint for a better comparison with TCCON following 
	# https://tccon-wiki.caltech.edu/Main/AuxiliaryDataGGG2020#Using_pressure_weights, as AK is evaluated around the WET VMR instead of DRY VMR
	# per level: x_ak = x_ap + int_operator * ak * (stilt_wet - ap_wet)
	# because stilt_wet = foot_wet * emission, 
	# xfoot = int_operator * ak * foot_dry * sf.wet
	# remind that int_operator will convert wet mole fraction to dry columns
	# -------------------------------------------------------------------------

	# calculate diff in pressure between particles based on release heights
	dp_par = abs(c(tccon.psfc - max(sel.p$xpres), diff(sel.p$xpres)))
	dp_tccon = abs(c(tccon.psfc - max(tccon.df$prior_pres),
					 diff(tccon.df$prior_pres)))	# in mb
	dp_df = data.frame(prior_pres = tccon.df$prior_pres, dp_tccon)

	# interpolate TCCON profiles to column receptors according to pressure
	rev.p = sel.p %>% 
			mutate(dp_tccon = approx(dp_df$prior_pres, dp_df$dp_tccon, 
									 xout = xpres, rule = 2)$y,
				   dp_par = dp_par, 
				   ak_norm = approx(tccon.df$ak_pres, tccon.df$ak_norm, 
									xout = xpres, rule = 2)$y, 
				   pwf = approx(tccon.df$prior_pres, tccon.df$int_op, 
								xout = xpres, rule = 2)$y * dp_par / dp_tccon,
				   ap_h2o_dry = approx(tccon.df$prior_pres, tccon.df$ap_h2o_dry,
									   xout = xpres, rule = 2)$y, 
				   ap_gas_wet = approx(tccon.df$prior_pres, tccon.df$ap_gas_wet,
									   xout = xpres, rule = 2)$y, 

				   # calculate the scaling factor that convert dry air mole fraction to wet air mole fraction to wet the modeled concentration/footprint 
				   # foot_wet = foot_dry / (1 + dry mole fraction of H2O)
				   # thus sf.wet = 1 / (1 + dry MF of H2O from TCCON prior)
				   sf.wet = 1 / (1 + ap_h2o_dry) ) %>% arrange(desc(xpres))


	# -------------------------------------------------------------------------
	# 3. lastly locate the OCO-2 levels that above the STILT particles 
	min.pres.xstilt = min(sel.p$xpres)
	upper.tccon.df = tccon.df %>% filter(prior_pres < min.pres.xstilt) %>% 
					 arrange(desc(prior_pres)) %>% mutate(sf.wet = 1) %>%
					 dplyr::select(ak_norm, pwf = int_op, sf.wet, 
					 			   pres = prior_pres, ap_gas_wet, ap_h2o_dry)	

	npar = max(rev.p$indx)
	upper.tccon.df$indx = seq(npar + 1, npar + nrow(upper.tccon.df))

	# adjust the int operator for the first level above XSTILT levels
	upper.tccon.df$pwf[1] = sum(tccon.df$int_op) - sum(rev.p$pwf) - 
							sum(upper.tccon.df$pwf[-1])

	# -------------------------------------------------------------------------
    # 4. Combine interpolated and initial TCCON profiles 
	lower.df = rev.p %>% dplyr::select(ak_norm, pwf, sf.wet, pres = xpres, 
									   ap_gas_wet, ap_h2o_dry, indx) 
	combine.prof = rbind(lower.df, upper.tccon.df) %>% 
				   rename(ak.norm = ak_norm) %>% mutate(ak.pwf = ak.norm * pwf)


	# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak,
	#AK*PWF and a priori contribution
	combine.prof$stiltTF = F
	combine.prof[combine.prof$indx <= npar, 'stiltTF'] = T

	combine.prof  	# return weighting profiles
}



if (F) {	# sanity check	
	p1 = ggplot() + theme_bw() + ylim(c(1100, 0)) + 
		 #geom_point(data = rev.p, aes(ak_norm, xpres), color = 'red') + 
		 #geom_point(data = rev.p, aes(int_op * 51, xpres),color = 'blue') +
		 geom_point(data = tccon.df, aes(ak_norm, ak_pres), color = 'red') + 
		 geom_point(data = tccon.df, aes(int_op , ak_pres), color = 'blue')

}

