
# convert TROPOMI CO AK to normalized AK
convert.tropomi.co.ak = function(output, co.path) {
	
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
	co.info = get.tropomi.prof(receptor, tropomi.speci = 'CO', tropomi.path = co.path)
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

