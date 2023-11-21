### subroutine to interpolate AK from TROPOMI and calc PWF for modeled levels
# Dien Wu, 08/26/2020

# TROPOMI CO has no aprior profiles, update on 09/05/2020
# NOTE THAT the pressure weighting function is calculated based on dp between 
# every two release level (dp) over dp from surface and TOA (Psfc - Ptoa), 
# i.e., PWF = dp_recp / (p_sfc - p_toa)

# for tropospheric PWF, need to swtich to use p_tropopause rather than p_toa
# TROPOMI v2 CO has now include the normalized AK profiles, no need for normalization, DW, 10/12/2023

get.wgt.tropomi.func = function(output, tropomi.fn, tropomi.speci, 
								tropo.no2TF = T){
	
	if (tropomi.speci != 'NO2') tropo.no2TF = FALSE 

	# before weighting trajec-level footprint by AK & PWF
	# get specific humidity and temp profiles that have been extracted via 
	# before_trajec_xstilt()
	qt.prof = output$qt_prof
	if (is.null(qt.prof)) 
		stop('get.wgt.tropomi.func(): no extracted q and T profiles found...\n')

    # grab receptor info and select the particles at first time step back
    receptor = output$receptor
	p = output$particle
	min.time = min(abs(p$time)) * sign(p$time[1])  # MIN time in mins
	sel.p = p[p$time == min.time, ] %>% 
			dplyr::select(indx, zagl, zsfc, pres, xhgt)

	# get useful info from TROPOMI column CO and/or NO2 data, DW, 09/05/2020 
	# e.g., surface pressure/height, AK, vertical pressures, retrieved obs...
	trp.info = get.tropomi.prof(receptor, tropomi.speci, tropomi.fn =tropomi.fn)

	# encountered NAs for AKs, AMFs, and NO2 VCDs from L2 NO2 PAL files, 
	# DW, Jan 9, 2023 
	if (tropomi.speci == 'NO2') {
		if (is.na(trp.info$no2_vcd_tropo)) {
			cat('get.wgt.tropomi.func(): error in encountering NA values for observed NO2 or AKs, AMFs...returning NULL\n'); return()
		}
	} else if ( tropomi.speci == 'CO') {
		if (is.na(trp.info$co_vcd)) {
			cat('get.wgt.tropomi.func(): error in encountering NA values for observed CO...returning NULL\n'); return()
		}
	} else if ( tropomi.speci == 'CH4') {
		if ( is.na(trp.info$xch4_bc) ) {
			cat('get.wgt.tropomi.func(): error in encountering NA values for observed CH4...returning NULL\n'); return()
		}
	}

	trp.zsfc = trp.info$tropomi_zsfc
	trp.psfc = trp.info$tropomi_psfc
	if (is.na(trp.zsfc)) stop('get.wgt.tropomi.func(): TROPOMI-based surface height is NA, please check...\n')

	# use tropospheric AK instead of total AK, DW, 09/23/2020
	if ('ak.norm' %in% names(trp.info) ) {
		trp.df = as.data.frame(trp.info[c('ak.norm', 'lower_pres', 
										  'upper_pres')]) 
	} else if ('ak' %in% names(trp.info)) {
		trp.df = as.data.frame(trp.info[c('ak', 'lower_pres', 'upper_pres')]) 
	}	# end if

	
	if ( tropomi.speci == 'NO2' & tropo.no2TF == T ) 
		trp.df = as.data.frame(trp.info[c('ak_tropo', 'lower_pres', 'upper_pres')]) 

	if ( NA %in% trp.df$lower_pres ) {
		
		# occasionally, NA values for TROPOMI XCO pressure levels...
		# create fake 50 pressure layers (51 pressure levels) with even dp using TROPOMI Psfc
		trp.pres.bound = seq(0.01, trp.psfc, length = 51)
		trp.df = data.frame(ak.norm = 1, 
							lower_pres = trp.pres.bound[2:51], 
							upper_pres = trp.pres.bound[1:50]) %>% 
					 mutate(ratio_pres = upper_pres / lower_pres)
	} else {

		# overwrite the pressure for the bottom level with tropomi sfc pressure
		trp.df = trp.df %>% 
				 mutate(lower_pres = ifelse(lower_pres == max(lower_pres), 
											trp.info$tropomi_psfc, lower_pres), 
						ratio_pres = upper_pres / lower_pres)
		trp.pres.bound = c(trp.df$upper_pres[1], trp.df$lower_pres)
	}	# end if


	# --------- correct for model-satellite mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RTv_mean * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# Tv_mean is the mean virtual temp 

	# correct for diff in ZSFC between sensor and particle ZSFC
	p.zagl.corr = sel.p$zagl + sel.p$zsfc - trp.zsfc	# corrected ZAGL
	p.pres = sel.p$pres 

	# use satellite surface altitude and pres to calc ZAGL that match TROPOMI 
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ trp.psfc * exp(a * (-p.zagl.corr)), 
					 start = list(a = 1E-4))		
	a.tropomi = coef(nls)[[1]]    # a3 stands for g/RTv_mean using OCO atmos-X

	# --------------------------------------------------------------------------
	# now start to weight trajec-level footprint for X-STILT
	# --------------------------------------------------------------------------
	# first bin up specific according to pressure of bottom level
	# +1 is to adjust as findInterval() finds the lower end of an interval
	# but pressure declines with height
	qt.bin = qt.prof %>% 
			 mutate(lower_pres = trp.pres.bound[findInterval(pres, trp.pres.bound) + 1]) %>% 
			 group_by(lower_pres) %>% summarise_all(mean) %>% 
			 ungroup() %>% dplyr::select(-pres) %>% 
			 right_join(trp.df, by = 'lower_pres') %>% 
			 arrange(lower_pres) %>%
			 mutate_if(is.numeric, ~ if_else(is.na(.x), 0, .x))

	# if for CO, we need to normalize AK for TROPOMI CO
	# calculcate the normalized AK, AK_m / dz, 
	# use P2  = P1 * exp (g/RTv_mean * dz) to calculate dz
	if ( 'CO' %in% tropomi.speci & !'ak.norm' %in% colnames(qt.bin) ) 
		qt.bin = qt.bin %>% mutate(dz = -log(ratio_pres) / a.tropomi, 
								   ak.norm = ak / dz)
	
	if ('NO2' %in% tropomi.speci) 
		qt.bin = qt.bin %>% rename(ak.norm = ak_tropo)
	
	if ('CH4' %in% tropomi.speci) 
		qt.bin = qt.bin %>% rename(ak.norm = ak)

	cat(paste('get.wgt.tropomi.func(): tropospheric AK_norm for the surface layer is', signif(qt.bin[qt.bin$lower_pres == max(qt.bin$lower_pres), 'ak.norm'], 3), '\n'))

	# correct q and T for bottom layers, when modeled Psfc < satellite Psfc
	if (max(qt.prof$pres) < max(qt.bin$lower_pres)) {
		qt.bin[qt.bin$lower_pres >= max(qt.prof$pres), 'sphu'] = max(qt.bin$sphu)

		qt.bin[qt.bin$lower_pres >= max(qt.prof$pres), 'temz'] = max(qt.bin$temz)
	}	# end if correction


	# --------------------------------------------------------------------------
	# estimate initial release height and pressure based on hyposmetric equation
	zmin = min(receptor$zagl)
	zmax = max(receptor$zagl)
	npar = max(p$indx)
	#xhgt = zmax - zmax / npar * (npar - indx)

	# use TROPOMI sfc pressure to calculate pressure of release heights, 
	# as model-obs Psfc could diff, use satellite Psfc to adjust column 
	# P = Psfc * exp (g/RT * (Zsfc - Z))
	sel.p = sel.p %>% mutate(lower_pres = trp.psfc * exp(a.tropomi*(-xhgt)))

	# -------------------------------------------------------------------------
	# identify the TROPOMI levels which locate above the max release levels
	# and merge them with upper TROPOMI levels
	upper.trp.df = trp.df %>% filter(lower_pres < min(sel.p$lower_pres)) %>% 
				   arrange(desc(lower_pres)) 
	
	# select tropospheric levels if simulating tropo column, DW, 07/08/2022
	if (tropo.no2TF) upper.trp.df = upper.trp.df %>% filter(ak_tropo > 0)
	upper.trp.df$indx = seq(npar + 1, npar + nrow(upper.trp.df))
	
	# total column with STILT for lower atmos and TROPOMI for upper
	combine.prof = rbind(sel.p[, c('indx', 'lower_pres')], 
						 upper.trp.df[, c('indx', 'lower_pres')])


	# ----------------------------------------------------------------------
	# calculate dry-air column density in mol m-2
	g = 9.8        	      	# m s-2
	Mdry = 29 / 1E3        	# 29 g mol-1 -> now in kg mol-1

	# FIRST interpolate specific humidity, AK_norm to each particle level 
	# based on its release pressure
	# always use press for the lower TROPOMI level for interpolation
	combine.prof = combine.prof %>% arrange(indx) %>% 
				   mutate(q = approx(qt.bin$lower_pres, qt.bin$sphu, 
									 xout = combine.prof$lower_pres, 
									 rule = 2)$y,
						
						  ak.norm = approx(qt.bin$lower_pres, 
										   qt.bin$ak.norm, 
										   xout = combine.prof$lower_pres, 
										   rule = 2)$y, 

						# calculate diff in pressure between particles based on release heights
						dp = abs(c(trp.psfc - max(lower_pres), 
								 diff(lower_pres))), 

						# calculate dry-air column density per particle/level
						xdry = (1 - q) / g / Mdry * dp * 100,	# mol m-2

						# Xdry for combined levels, either tropo or total X
						# see Line 137!!! - for tropospheric column
						# upper TROPOMI levels above STILT levels only extend
						# to tropopause --- 
						xdry.tot = sum(xdry), 		# mol m-2
						
						# PWF is normalized by air density for tropo or total X
						pwf = xdry / xdry.tot )
	
	# -------------------------------------------------------------------------
	# adjust the PWF for the first level above XSTILT particles 
	combine.prof$pwf[npar + 1] = 1 - sum(combine.prof$pwf[-(npar + 1)])
	combine.prof$ak.pwf = combine.prof$ak.norm * combine.prof$pwf

	# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak,
	#AK*PWF and a priori contribution
	combine.prof$stiltTF = F
	combine.prof[combine.prof$indx <= npar, 'stiltTF'] = T

	return(list(combine.prof = combine.prof, trp.info = trp.info))  	
	# return weighting profiles
}

# end of subroutine
