### subroutine to interpolate AK from TROPOMI and calc PWF for modeled levels
# Dien Wu, 08/26/2020

# TROPOMI CO has no aprior profiles, update on 09/05/2020
# NOTE THAT the pressure weighting function is calculated based on dp between 
# every two release level (dp) over dp from surface and TOA (Psfc - Ptoa), 
# i.e., PWF = dp_recp / (p_sfc - p_toa)

# for tropospheric PWF, need to swtich to use p_tropopause rather than p_toa
get.wgt.tropomi.func = function(output, tropomi.fn, tropomi.speci, 
								tropo.no2TF = T){
	
	# before weighting trajec-level footprint by AK & PWF
	# get specific humidity and temp profiles that have been extracted via 
	# before_trajec_xstilt()
	qt.prof = output$qt_prof
	if (is.null(qt.prof)) stop('get.wgt.tropomi.func(): no extracted q and T profiles found...\n')

    # grab receptor info and select the particles at first time step back
    receptor = output$receptor
	p = output$particle
	min.time = min(abs(p$time)) * sign(p$time[1])  # MIN time in mins
	sel.p = p[p$time == min.time, ] %>% 
			dplyr::select(indx, zagl, zsfc, pres, xhgt)

	
	# get useful info from TROPOMI column CO and/or NO2 data, DW, 09/05/2020 
	# e.g., surface pressure/height, AK, vertical pressures, retrieved obs...
	tropomi.info = get.tropomi.prof(receptor, tropomi.speci, 
									tropomi.fn = tropomi.fn)

	# use tropospheric AK instead of total AK, DW, 09/23/2020
	tropomi.df = as.data.frame(tropomi.info[c('ak','lower_pres','upper_pres')]) 
	if (tropomi.speci == 'NO2' & tropo.no2TF == T) 
		tropomi.df = as.data.frame(tropomi.info[c('ak_tropo', 'lower_pres', 'upper_pres')]) 
	
	# overwrite the pressure for the bottom level with sfc pressure
	tropomi.df = tropomi.df %>% 
				 mutate(lower_pres = ifelse(lower_pres == max(lower_pres), 
						   					tropomi.info$tropomi_psfc, 
											lower_pres), 
					  	ratio.pres = upper_pres / lower_pres)

	tropomi.pres.bound = c(tropomi.df$upper_pres[1], tropomi.df$lower_pres)
	tropomi.zsfc = tropomi.info$tropomi_zsfc
	tropomi.psfc = tropomi.info$tropomi_psfc
	if (is.na(tropomi.zsfc)) stop('get.wgt.tropomi.func(): TROPOMI-based surface height is NA, please check...\n')
	

	# --------- correct for model-satellite mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RTv_mean * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# Tv_mean is the mean virtual temp 

	# correct for diff in ZSFC between sensor and particle ZSFC
	p.zagl.corr = sel.p$zagl + sel.p$zsfc - tropomi.zsfc	# corrected ZAGL
	p.pres = sel.p$pres 

	# use satellite surface altitude and pres to calc ZAGL that match TROPOMI 
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ tropomi.psfc * exp(a * (-p.zagl.corr)), 
					 start = list(a = 1E-4))		
	a.tropomi = coef(nls)[[1]]    # a3 stands for g/RTv_mean using OCO atmos-X

	# --------------------------------------------------------------------------
	# now start to weight trajec-level footprint for X-STILT
	# --------------------------------------------------------------------------
	# first bin up specific according to pressure of bottom level
	# +1 is adjust as findInterval() finds the lower end of a pressure interval
	# but pressure declines with height
	qt.bin = qt.prof %>% 
			 mutate(lower_pres = tropomi.pres.bound[findInterval(pres, tropomi.pres.bound) + 1]) %>% 
			 group_by(lower_pres) %>% summarise_all(mean) %>% ungroup() %>% 
			 dplyr::select(-pres) %>% 
			 right_join(tropomi.df, by = 'lower_pres') %>% 
			 arrange(lower_pres) %>% mutate_all(~ if_else(is.na(.x), 0, .x))


	# if for CO, we need to normalize AK for TROPOMI CO
	# calculcate the normalized AK, AK_m / dz, 
	# use P2  = P1 * exp (g/RTv_mean * dz) to calculate dz
	if ('CO' %in% tropomi.speci) 
		qt.bin = qt.bin %>% mutate(dz = -log(ratio.pres) / a.tropomi, 
								   ak.norm = ak / dz)
	if ('NO2' %in% tropomi.speci) qt.bin = qt.bin %>% rename(ak.norm = ak_tropo)
	if ('CH4' %in% tropomi.speci) qt.bin = qt.bin %>% rename(ak.norm = ak)

	cat(paste('get.wgt.tropomi.func(): tropospheric AK_norm for the surface layer is', signif(qt.bin[qt.bin$lower_pres == max(qt.bin$lower_pres), 'ak.norm'], 3), '\n'))

	# --------------------------------------------------------------------------
	# estimate initial release height and pressure based on hyposmetric equation
	zmin = min(receptor$zagl)
	zmax = max(receptor$zagl)
	npar = max(p$indx)
	#xhgt = zmax - zmax / npar * (npar - indx)

	# use TROPOMI sfc pressure to calculate pressure of release heights, 
	# as model-obs Psfc could diff, always use satellite Psfc to adjust column 
	# P = Psfc * exp (g/RT * (Zsfc - Z))
	sel.p = sel.p %>% mutate(lower_pres = tropomi.psfc*exp(a.tropomi*(-xhgt))) 


	# identify the TROPOMI levels which locate above the max release levels
	# and merge them with upper TROPOMI levels
	upper.tropomi.df = tropomi.df %>% 
					   filter(lower_pres < min(sel.p$lower_pres)) %>% 
					   arrange(desc(lower_pres)) 
	upper.tropomi.df$indx = seq(npar + 1, npar + nrow(upper.tropomi.df))
	combine.prof = rbind(sel.p[, c('indx', 'lower_pres')], 
						 upper.tropomi.df[, c('indx', 'lower_pres')])
	

	# -------------------------------------------------------------------------
	# calculate dry-air column density in mol m-2
	g = 9.8        	      	# m s-2
	Mdry = 29 / 1E3        	# kg mol-1
	
	# interpolate specific humidity, AK_norm to each particle based on its release pressure
	# always use press for lower TROPOMI level as a reference for interpolation
	combine.prof = combine.prof %>% arrange(indx) %>% 
				   mutate(q = approx(qt.bin$lower_pres, qt.bin$sphu, 
									 xout = lower_pres, rule = 2)$y, 

						  ak.norm = approx(qt.bin$lower_pres, qt.bin$ak.norm, 
						   				   xout = lower_pres, rule = 2)$y, 
							
						  # calculate diff in pressure between particles based on release heights
						  dp = abs(c(tropomi.psfc - max(lower_pres), 
						  			diff(lower_pres))), 

						  # finally calculate dry-air column density per particle/level
						  xdry = (1 - q) / g / Mdry * dp * 100,	# mol m-2
						  xdry.tot = sum(xdry), 	# Xdry for all stilt levels
						  
						  # PWF is normalized by air X-density for total X
						  pwf = xdry / xdry.tot)

	# adjust the PWF for the first level above XSTILT particles 
	combine.prof$pwf[npar + 1] = 1 - sum(combine.prof$pwf[-(npar + 1)])
	combine.prof$ak.pwf = combine.prof$ak.norm * combine.prof$pwf

	# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak,
	#AK*PWF and a priori contribution
	combine.prof$stiltTF = F
	combine.prof[combine.prof$indx <= npar, 'stiltTF'] = T

	return(list(combine.prof = combine.prof))  	
	# return weighting profiles
}

# end of subroutine
