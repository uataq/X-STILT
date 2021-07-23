### different from get.wgt.oco.func() or get.wgt.tropomi.func(), 
# this subroutine will use the modeled surface height and pressure to 
# perform pressure weighting, no dependence on any satellite sensors

# simply treat normalized AK as 1, DW, 09/15/2020

get.wgt.mod.func = function(output){

	# before weighting trajec-level footprint by AK & PWF
	# get specific humidity and temp profiles that have been extracted via 
	# before_trajec_xstilt()
	qt.prof = output$qt_prof
	if (is.null(qt.prof)) stop('get.wgt.mod.func(): no extracted q and T profiles found...\n')

	# grab receptor info and select the particles at first time step back
	receptor = output$receptor; p = output$particle
	min.time = min(abs(p$time)) * sign(p$time[1])  # MIN time in mins
	sel.p = p[p$time == min.time, ] %>% dplyr::select(indx, zagl, zsfc, pres, xhgt)


	# --------- correct for model-satellite mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RTv_mean * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# Tv_mean is the mean virtual temp 
	# correct for diff in ZSFC between sensor and particle ZSFC
	p.zagl = sel.p$zagl
	p.pres = sel.p$pres 				# trajec-level pressure
	p.psfc = receptor$psfc 

	# use satellite surface altitude and pressure to calculate ZAGL that TROPOMI believes
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ p.psfc * exp(a * (-p.zagl)), start = list(a = 1E-4))		
	a.mod = coef(nls)[[1]]    # a stands for g/RTv_mean based on modeled atmos column

	# check correlation coefficient
	#cor(p.pres, predict(nls))	# 0.99989

	# ----------------------------------------------------------------------------
	# now start to weight trajec-level footprint for X-STILT
	# ----------------------------------------------------------------------------
	# use modeled sfc pressure to calculate corresponding pressure of release heights
	# P = Psfc * exp (g/RT * (Zsfc - Z))
	sel.p$xpres = p.psfc * exp(a.mod * (-sel.p$xhgt))
	pres.bound = sort(unique(sel.p$xpres))     # needs to be in increasing trend
	qt.bin = qt.prof %>% mutate(xpres = pres.bound[findInterval(pres, pres.bound) + 1]) %>% 
						 group_by(xpres) %>% summarise_all(mean) %>% ungroup() %>% 
						 dplyr::select(-pres) %>% na.omit()

	# interpolate specific humidity 
	sel.p = sel.p %>% mutate(q = approx(qt.bin$xpres, qt.bin$sphu, xout = xpres, rule = 2)$y) %>% 
			arrange(desc(xpres))

	# calculate diff in pressure between particles based on release heights
	dp = abs(c(p.psfc - max(sel.p$xpres), diff(sel.p$xpres)))


	# merge all info per particle and calculate dry-air column density in mol m-2
	g = 9.8        	      	# m s-2
	Mdry = 29 / 1E3        	# kg mol-1

	# dry air column density in mol m-2, 100 for converting pres from hPa to Pa
	# need to estimate the total column water vapor, 
	entire.dp = 10		# create pressure levels with even dp of 10 hPa
	entire.pres = seq(0, p.psfc, 10)

	# get q for each one level
	qt.entire = qt.prof %>% mutate(xpres = entire.pres[findInterval(pres, entire.pres) + 1]) %>% 
				group_by(xpres) %>% summarise_all(mean) %>% ungroup() %>% na.omit()	

	# get modeled dry-air column density for the entire column and for each particle
	xdry.tot = sum((1 - qt.entire$sphu) / g / Mdry * entire.dp * 100)
	sel.p = sel.p %>% mutate(xdry = (1 - q) / g / Mdry * dp * 100, pwf = xdry / xdry.tot)

	# all levels are modeled levels, since we removed the dependence from satellites, 
	# no apriori, AK_norm is simply 1. 
	xstilt.prof = sel.p %>% mutate(ak.norm = 1, ak.pwf = ak.norm * pwf, stiltTF = T) %>% 
							dplyr::select(indx, ak.norm, pwf, ak.pwf)
							 
	xstilt.prof  	# return weighting profiles

}	# end of subroutine
