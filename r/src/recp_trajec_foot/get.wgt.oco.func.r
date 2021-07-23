### subroutine to get new AK PWF profiles based on OCO-2's
# pressure weighting and averaging kernel, for combined levels
# OCO-2 only provides the PWF, AK and a priori at 20 levels,
# use linear interpolation to "approx" values at given STILT releasing levels
# Dien Wu, 08/05/2016

# Bugs fixed, DW:
# fix 1, 11/20/2016, modify the interpolate of PWF from OCO2 to STILT
# fix 2, 11/28/2016, remove storeTF flag, as we always need .RData file stored
#                    for generating footprint using Trajectfoot()

# fix 3, 02/08/2017, add uneven vertical spacing
# fix 4, 04/06/2017, change AK from zero to original OCO2 AK above model level
# add 5, 04/19/2017, store interpolated AK PWF apriori

# version 3 for matching Ben's STILT-R version 2, DW, 05/25/2018
# interpolate ground hgt in this subroutine, DW, 05/25/2018
# output refers to all the content from .rds file using STILTv2, 05/25/2018
# which is the same 'output' from simulation_step()

# *** USE OCO-2 retrieved surface pressure/altitude for air column 
# because modeled Psurf/air column differs from retrieved Psurf/air column
# DW, 08/30/2019 

# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# interpolate vertical profiles onto each particle instead of each release level, 
#    due to changes in HYSPLIT compiler, DW, 07/07/2020

get.wgt.oco.func = function(output, oco.path, oco.fn = NA) {

	# before weighting trajec-level footprint by AK & PWF
	# get specific humidity and temp profiles that have been extracted via 
	# before_trajec_xstilt()
	qt.prof = output$qt_prof
	if (is.null(qt.prof)) stop('get.wgt.oco.func(): no extracted q and T profiles found...\n')

    # grab receptor info and select the particles at first time step back
    receptor = output$receptor; p = output$particle
	min.time = min(abs(p$time)) * sign(p$time[1])  # MIN time in mins
	sel.p = p[p$time == min.time, ] %>% dplyr::select(indx, zagl, zsfc, pres, xhgt)

	# get OCO-2/3 profile first according to lat/lon of receptor, return a list
	if (!is.na(oco.fn)) oco.path = NA 
	oco.info = get.oco.info(oco.path, receptor, oco.fn)
	if (is.null(oco.info)) {
		warning('get.wgt.oco.func(): NO OCO info found for this receptor...please check\n')
		return()
	} # end if is.null


	#### ------------------------ DEALING WITH OCO NOW -------------------- ####
    # grab press weighting function, pressures, normalized AK, apriori
	oco.df = oco.info[c('ak.norm', 'pwf', 'pres', 'ap')] %>% as.data.frame() 
	oco.pres.bound = as.numeric(oco.df$pres)
	oco.ak.norm = as.numeric(oco.df$ak.norm)
	oco.zsfc = oco.info$oco.grdhgt
	oco.psfc = oco.info$oco.psfc
	oco.xh2o = oco.info$oco.xh2o	# already in mol m-2
	if (is.na(oco.zsfc)) stop('get.wgt.oco.func(): satellite-based surface height is NA, please check...\n')
	

	# --------- correct for model-satellite mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RTv_mean * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# Tv_mean is the mean virtual temp 

	# correct for diff in ZSFC between sensor and particle ZSFC
	p.zasl = sel.p$zagl + sel.p$zsfc 	# trajec-level ASL
	p.pres = sel.p$pres 				# trajec-level pressure
	p.zagl.corr.co2 = p.zasl - oco.zsfc	# corrected trajec-level AGL

	# use satellite surface altitude and pressure to calculate ZAGL that TROPOMI believes
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ oco.psfc * exp(a * (-p.zagl.corr.co2)), start = list(a = 1E-4))		
	a.oco = coef(nls)[[1]]    # a3 stands for g/RTv_mean based on OCO atmos-X


	# ----------------------------------------------------------------------------
	# now start to weight trajec-level footprint for X-STILT
	# ----------------------------------------------------------------------------
	# first bin up specific according to pressure of bottom level
	# +1 is adjust as findInterval() finds the lower end of a pressure interval
	# but pressure declines with height
	qt.bin = qt.prof %>% 
			  mutate(xpres = oco.pres.bound[findInterval(pres, oco.pres.bound) + 1]) %>% 
			  group_by(xpres) %>% summarise_all(mean) %>% ungroup() %>% 
			  dplyr::select(-pres) %>% na.omit()
	

	# estimate initial release height and pressure based on hyposmetric equation
	zmin = min(receptor$zagl)
	zmax = max(receptor$zagl)
	npar = max(p$indx)

	# HYSPLITv5 deploys a line source (receptor$zagl) for releasing particles, 
	#xhgt = zmax - zmax / npar * (npar - indx), already incorporated in STILT
	# release height of each particle is roughly distributed between levels of a line source
	sel.p = sel.p %>% 

			 # use OCO sfc pressure to calculate pressure of release heights
			 # P = Psfc * exp (g/RT * (Zsfc - Z))
			 mutate(xpres = oco.psfc * exp(a.oco * (-xhgt)), 
			
					# interpolate specific humidity, AK_norm, apriori to 
					# each particle based on its release pressure
			 		q = approx(qt.bin$xpres, qt.bin$sphu, xout = xpres, rule = 2)$y, 
			        ak.norm = approx(oco.df$pres, oco.df$ak.norm, xout = xpres, rule = 2)$y, 
					ap = approx(oco.df$pres, oco.df$ap, xout = xpres, rule = 2)$y) %>% 

			arrange(desc(xpres))

	# calculate diff in pressure between particles based on release heights
	dp = abs(c(oco.psfc - max(sel.p$xpres), diff(sel.p$xpres)))

	
	# merge all info per particle and calculate dry-air column density in mol m-2
	g = 9.8        	      	# m s-2
	Mdry = 29 / 1E3        	# kg mol-1

	# dry air column density in mol m-2, 100 for converting pres from hPa to Pa
	xdry.tot = 1 / g / Mdry * oco.psfc * 100 - oco.xh2o   
	sel.p = sel.p %>% mutate(xdry = (1 - q) / g / Mdry * dp * 100, pwf = xdry / xdry.tot)

	# lastly locate the OCO-2 levels that above the STILT particles 
	min.pres.xstilt = min(sel.p$xpres)
	upper.oco.df = oco.df %>% filter(pres < min.pres.xstilt) %>% arrange(desc(pres))
	upper.oco.df$indx = seq(npar + 1, npar + nrow(upper.oco.df))

	# adjust the PWF for the first level above XSTILT particles 
	upper.oco.df$pwf[1] = 1 - sum(sel.p$pwf) - sum(upper.oco.df$pwf[-1])
	cat('get.wgt.oco.func(): Done computing PWF...\n')


    ## ----- Combine all interpolated OCO-2 profiles and calculating AK_norm *PWF
	lower.df = sel.p %>% dplyr::select(ak.norm, pwf, xpres, ap, indx) %>% rename(pres = xpres)
	combine.prof = rbind(lower.df, upper.oco.df) %>% mutate(ak.pwf = ak.norm * pwf)


	# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak,
	#AK*PWF and a priori contribution
	combine.prof$stiltTF = F
	combine.prof[combine.prof$indx <= npar, 'stiltTF'] = T

	combine.prof  	# return weighting profiles
}

# end of subroutine




