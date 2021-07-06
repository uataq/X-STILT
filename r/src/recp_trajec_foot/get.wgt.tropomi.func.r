### subroutine to interpolate AK PWF profiles based on TROPOMI's averaging kernel
# Dien Wu, 08/26/2020

# TROPOMI CO has no aprior profiles, update on 09/05/2020

get.wgt.tropomi.func = function(output, tropomi.fn, tropomi.speci, tropo.no2TF = T){
	
	# before weighting trajec-level footprint by AK & PWF
	# get specific humidity and temp profiles that have been extracted via 
	# before_trajec_xstilt()
	qt.prof = output$qt_prof
	if (is.null(qt.prof)) stop('get.wgt.tropomi.func(): no extracted q and T profiles found...\n')

    # grab receptor info and select the particles at first time step back
    receptor = output$receptor
	p = output$particle
	min.time = min(abs(p$time)) * sign(p$time[1])  # MIN time in mins
	sel.p = p[p$time == min.time, ] %>% dplyr::select(indx, zagl, zsfc, pres, xhgt)

	
	# get useful info from TROPOMI column CO and/or NO2 data, DW, 09/05/2020 
	# e.g., surface pressure/height, AK, vertical pressures, retrieved obs...
	tropomi.info = get.tropomi.prof(receptor, tropomi.speci, tropomi.fn = tropomi.fn)

	# use tropospheric AK instead of total AK, DW, 09/23/2020
	tropomi.df = as.data.frame(tropomi.info[c('ak', 'lower.pres', 'upper.pres')]) 
	if (tropomi.speci == 'NO2' & tropo.no2TF == T) 
		tropomi.df = as.data.frame(tropomi.info[c('ak.tropo', 'lower.pres', 'upper.pres')]) 
	
	# overwrite the pressure for the bottom level with sfc pressure
	tropomi.df = tropomi.df %>% mutate(lower.pres = ifelse(lower.pres == max(lower.pres), 
						   					                tropomi.info$tropomi.psfc, 
															lower.pres), 
					  	 				ratio.pres = upper.pres / lower.pres)

	tropomi.pres.bound = c(tropomi.df$upper.pres[1], tropomi.df$lower.pres)
	tropomi.zsfc = tropomi.info$tropomi.zsfc
	tropomi.psfc = tropomi.info$tropomi.psfc
	if (is.na(tropomi.zsfc)) stop('get.wgt.tropomi.func(): TROPOMI-based surface height is NA, please check...\n')
	

	# --------- correct for model-satellite mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RTv_mean * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# Tv_mean is the mean virtual temp 

	# correct for diff in ZSFC between sensor and particle ZSFC
	p.zagl.corr = sel.p$zagl + sel.p$zsfc - tropomi.zsfc	# corrected trajec-level AGL
	p.pres = sel.p$pres 

	# use satellite surface altitude and pressure to calculate ZAGL that TROPOMI believes
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ tropomi.psfc * exp(a * (-p.zagl.corr)), start = list(a = 1E-4))		
	a.tropomi = coef(nls)[[1]]    # a3 stands for g/RTv_mean based on OCO atmos-X

	# ----------------------------------------------------------------------------
	# now start to weight trajec-level footprint for X-STILT
	# ----------------------------------------------------------------------------
	# first bin up specific according to pressure of bottom level
	# +1 is adjust as findInterval() finds the lower end of a pressure interval
	# but pressure declines with height
	qt.bin = qt.prof %>% 
			  mutate(lower.pres = tropomi.pres.bound[findInterval(pres, tropomi.pres.bound) + 1]) %>% 
			  group_by(lower.pres) %>% summarise_all(mean) %>% ungroup() %>% 
			  dplyr::select(-pres) %>% right_join(tropomi.df, by = 'lower.pres') %>% 
			  arrange(lower.pres) %>% mutate_all(~ if_else(is.na(.x), 0, .x))


	# if for CO, we need to normalize AK for TROPOMI CO
	# calculcate the normalized AK, AK_m / dz, 
	# use P2  = P1 * exp (g/RTv_mean * dz) to calculate dz
	if ('CO' %in% tropomi.speci) 
		qt.bin = qt.bin %>% mutate(dz = -log(ratio.pres) / a.tropomi, ak.norm = ak / dz)
	if ('NO2' %in% tropomi.speci) qt.bin = qt.bin %>% rename(ak.norm = ak.tropo)
	if ('CH4' %in% tropomi.speci) qt.bin = qt.bin %>% rename(ak.norm = ak)

	cat(paste('get.wgt.tropomi.func(): tropospheric AK_norm for the surface layer is', 
			  signif(qt.bin[qt.bin$lower.pres == max(qt.bin$lower.pres), 'ak.norm'], 3), '\n'))

	# ----------------------------------------------------------------------------
	# estimate initial release height and pressure based on hyposmetric equation
	zmin = min(receptor$zagl)
	zmax = max(receptor$zagl)
	npar = max(p$indx)
	#xhgt = zmax - zmax / npar * (npar - indx)
	# use OCO sfc pressure to calculate pressure of release heights
	# P = Psfc * exp (g/RT * (Zsfc - Z))
	sel.p = sel.p %>% mutate(lower.pres = tropomi.psfc * exp(a.tropomi * (-xhgt))) 


	# locate the OCO-2 levels that above the STILT particles 
	# and merge them with upper TROPOMI levels
	upper.tropomi.df = tropomi.df %>% filter(lower.pres < min(sel.p$lower.pres)) %>% 
									   arrange(desc(lower.pres)) 
	upper.tropomi.df$indx = seq(npar + 1, npar + nrow(upper.tropomi.df))
	combine.prof = rbind(sel.p[, c('indx', 'lower.pres')], 
						  upper.tropomi.df[, c('indx', 'lower.pres')]) %>% arrange(indx)
	

	# -------------------------------------------------------------------------
	# calculate dry-air column density in mol m-2
	g = 9.8        	      	# m s-2
	Mdry = 29 / 1E3        	# kg mol-1
	
	# interpolate specific humidity, AK_norm to each particle based on its release pressure
	# always use press for lower TROPOMI level as a reference for interpolation
	combine.prof = combine.prof %>% 
					mutate(q = approx(qt.bin$lower.pres, qt.bin$sphu, 
									  xout = lower.pres, rule = 2)$y, 

						   ak.norm = approx(qt.bin$lower.pres, qt.bin$ak.norm, 
						   					xout = lower.pres, rule = 2)$y, 
							
						   # calculate diff in pressure between particles based on release heights
						   dp = abs(c(tropomi.psfc - max(lower.pres), diff(lower.pres))), 

						   # finally calculate dry-air column density per particle/level
						   xdry = (1 - q) / g / Mdry * dp * 100, 		  # mol m-2
						   xdry.tot = sum(xdry), pwf = xdry / xdry.tot)

	# also store TROPOMI info
	if ('CO' %in% tropomi.speci) 
		tropomi.info = tropomi.info[c('tropomi.lat', 'tropomi.lon', 'tropomi.zsfc', 
									   'tropomi.psfc', 'xco', 'xco.uncert')] %>% 
						as.data.frame() %>% 
						mutate(xdry.tot = unique(combine.prof$xdry.tot), 
							   xco.ppb = xco / xdry.tot * 1E9, 
							   xco.uncert.ppb = xco.uncert / xdry.tot * 1E9)

	if ('NO2' %in% tropomi.speci) 
		tropomi.info = tropomi.info[c('tropomi.lat', 'tropomi.lon', 'tropomi.zsfc', 
									   'tropomi.psfc', 'xno2.tropo', 'xno2.tropo.uncert')] %>% 
						as.data.frame() %>% 
						mutate(xdry.tot = unique(combine.prof$xdry.tot), 
							   xno2.tropo.ppb = xno2.tropo / xdry.tot * 1E9, 
							   xno2.tropo.uncert.ppb = xno2.tropo.uncert / xdry.tot * 1E9)

	if ('CH4' %in% tropomi.speci) 
		tropomi.info = tropomi.info[c('tropomi.lat', 'tropomi.lon', 
									  'tropomi.zsfc', 'tropomi.psfc', 
									  'xch4', 'xch4.bc', 'xch4.uncert')] %>% 
						as.data.frame() %>% 
						mutate(xdry.tot = unique(combine.prof$xdry.tot)) %>% 
						rename(xch4.ppb = xch4, xch4.bc.ppb = xch4.bc, 
							   xch4.uncert.ppb = xch4.uncert)


	# adjust the PWF for the first level above XSTILT particles 
	combine.prof$pwf[npar + 1] = 1 - sum(combine.prof$pwf[-(npar + 1)])
	combine.prof$ak.pwf = combine.prof$ak.norm * combine.prof$pwf

	# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak,
	#AK*PWF and a priori contribution
	combine.prof$stiltTF = F
	combine.prof[combine.prof$indx <= npar, 'stiltTF'] = T

	return(list(combine.prof = combine.prof, tropomi.info = tropomi.info))  	
	# return weighting profiles
}

# end of subroutine




# convert TROPOMI CO AK to normalized AK
convert.tropomi.co.ak = function(output, co.fn) {
	
	# before weighting trajec-level footprint by AK & PWF
	# get specific humidity and temp profiles that have been extracted via 
	# before_trajec_xstilt()
	qt.prof = output$qt_prof
	if (is.null(qt.prof)) stop('convert.tropomi.co.ak(): no extracted q and T profiles found...\n')

    # grab receptor info and select the particles at first time step back
    receptor = output$receptor
	p = output$particle
	min.time = min(abs(p$time)) * sign(p$time[1])  # MIN time in mins
	sel.p = p[p$time == min.time, ] %>% dplyr::select(indx, zagl, zsfc, pres, xhgt)
	
	# get useful info from TROPOMI column CO and/or NO2 data, DW, 09/05/2020 
	# e.g., surface pressure/height, AK, vertical pressures, retrieved obs...
	co.info = get.tropomi.prof(receptor, 'CO', co.fn = co.fn)
	co.df = as.data.frame(co.info[c('ak', 'lower.pres', 'upper.pres')]) %>% 
			 mutate(lower.pres = ifelse(lower.pres == max(lower.pres), 
						   				co.info$tropomi.psfc, lower.pres), 
					ratio.pres = upper.pres / lower.pres)

	co.pres.bound = c(co.df$upper.pres[1], co.df$lower.pres)
	co.zsfc = co.info$tropomi.zsfc
	co.psfc = co.info$tropomi.psfc
	if (is.na(co.zsfc)) stop('convert.tropomi.co.ak(): TROPOMI-based surface height is NA, please check...\n')
	

	# --------- correct for model-satellite mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RTv_mean * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# Tv_mean is the mean virtual temp 

	# correct for diff in ZSFC between sensor and particle ZSFC
	p.zagl.corr = sel.p$zagl + sel.p$zsfc - co.zsfc	# corrected trajec-level AGL
	p.pres = sel.p$pres 

	# use satellite surface altitude and pressure to calculate ZAGL that TROPOMI believes
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ co.psfc * exp(a * (-p.zagl.corr)), start = list(a = 1E-4))		
	a.co = coef(nls)[[1]]    # a3 stands for g/RTv_mean based on OCO atmos-X

	# ----------------------------------------------------------------------------
	# now start to weight trajec-level footprint for X-STILT
	# ----------------------------------------------------------------------------
	# first bin up specific according to pressure of bottom level
	# +1 is adjust as findInterval() finds the lower end of a pressure interval
	# but pressure declines with height
	qt.bin = qt.prof %>% 
			  mutate(lower.pres = co.pres.bound[findInterval(pres, co.pres.bound) + 1]) %>% 
			  group_by(lower.pres) %>% summarise_all(mean) %>% ungroup() %>% 
			  dplyr::select(-pres) %>% right_join(co.df, by = 'lower.pres') %>% 
			  arrange(lower.pres) %>% mutate_all(~ if_else(is.na(.x), 0, .x))


	# if for CO, we need to normalize AK for TROPOMI CO
	# calculcate the normalized AK, AK_m / dz, 
	# use P2  = P1 * exp (g/RTv_mean * dz) to calculate dz
	qt.bin = qt.bin %>% mutate(dz = -log(ratio.pres) / a.co, ak.norm = ak / dz)


	cat(paste('get.wgt.co.func(): tropospheric AK_norm for the surface layer is', 
			  signif(qt.bin[qt.bin$lower.pres == max(qt.bin$lower.pres), 'ak.norm'], 3), '\n'))

	qt.bin
}

