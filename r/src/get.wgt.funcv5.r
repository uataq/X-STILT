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

# add 6, 04/20/2017, add control flags "pwf.wgt" for weighting trajec
# (1) default: return trajec weighted by AK and PWF, ak.wgt == T & pwf.wgt == T
# (2) ak.wgt == F & pwf.wgt == T for only press wgt
# (3) ak.wgt == T & pwf.wgt == F for only normalized AK wgt
# (4) ak.wgt == F & pwf.wgt == F no wgt at all (returns the original trajec)

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


get.wgt.funcv5 <- function(output, oco.info, ak.wgt = T, pwf.wgt = T){

	if (ak.wgt) {
		cat('Turn on weighting OCO-2 simulation using averaging kernel...\n')
	} else cat('NO averaging kernel weighting, set ak to 1...\n') 


    # grab trajectory and receptor info
    receptor <- output$receptor; p <- output$particle
	min.time <- min(abs(p$time)) * sign(p$time[1])  # MIN time in mins

	# select the particles at first time step back
	sel.p <- p[p$time == min.time, ] %>% dplyr::select(indx, zagl, zsfc, pres)

    
	#### ------------------------ DEALING WITH OCO NOW -------------------- ####
    # grab press weighting function, pressures, normalized AK, apriori
	oco.ap      <- oco.info$ap
	oco.pwf     <- oco.info$pwf
	oco.pres    <- oco.info$pres
    oco.ak.norm <- oco.info$ak.norm
	oco.zsfc    <- oco.info$oco.grdhgt
	oco.psfc    <- max(oco.pres)
	if (is.na(oco.zsfc)) stop('get.wgt.funcv5(): satellite-based surface height is NA, please check...\n')
	

	# --------- correct for model-satellite mismatch in sfc P or Z --------- #
	# inferred from trajec, DW, 08/30/2019 
	# P = Psfc * exp (g/RT * (Zsfc - Z)), fit nonlinear regression between P and delta_Z,
	# since T_avg is unknown here
	p.zasl <- sel.p$zagl + sel.p$zsfc 
	p.zagl.rev <- sel.p$zagl + sel.p$zsfc - oco.zsfc	# correct for diff in ZSFC
	p.pres <- sel.p$pres 

	# use OCO surface altitude and pressure to calculate ZAGL that OCO believes
	# use particle-level pressure and ASL to calculate the coeffient
	nls.model <- stats::nls(p.pres ~ oco.psfc * exp(a * (-p.zagl.rev)), 
	                        start = list(a = 1E-4))
	a <- coef(nls.model)[[1]]
	# check correlation coefficient
	#cor(par.pres, predict(nls.model))	# test case: 0.9999208

	# calculate the pressure that OCO believe for each particle 
	sel.p.rev <- sel.p %>% mutate(pres_rev = oco.psfc * exp(a * (-zagl))) %>% 
					   	   arrange(desc(pres_rev)) 
	uni.pres <- unique(sel.p.rev$pres_rev)	
	cat(paste('OCO retrieved Psfc:', signif(oco.psfc, 5), 
			  '; modeled Psfc from met field:', signif(max(p.pres), 5), 
			  '; corrected Psfc:', signif(max(uni.pres), 5)))
	# now the highest pressure values should match the surface pressure retrieved from OCO

	#plot(sel.p.rev$indx, sel.p.rev$pres)
	#points(sel.p.rev$indx, sel.p.rev$pres_rev, col = 'orange')


	### ---------------------- start interpolation ------------------------ ###
	# since STILT display particles from surface-to-space, the opposite as OCO-2
	# product originally, profiles from OCO are from levels 1-20, from TOA to sfc.
	# we need to reverse AK and PWF profiles, as well as renaming attributes as
	# pressure levels. Thus, flipped profiles (now from level 1 to 20) will be
	# from sfc to TOA now...
	oco.pres <- oco.pres[length(oco.pres):1]	# flip pressure levels
	attributes(oco.pres)$names <- attributes(oco.pres)$names[length(attributes(oco.pres)$names):1]

	# flip 20 levels (--> from sfc to TOA) and assign names
	oco.ak.norm <- oco.ak.norm[length(oco.ak.norm):1]; attributes(oco.ak.norm)$names <- oco.pres
	oco.ap  <- oco.ap[length(oco.ap):1]; attributes(oco.ap)$names <- oco.pres
	oco.pwf <- oco.pwf[length(oco.pwf):1]; attributes(oco.pwf)$names <- oco.pres
	# plot(oco.pwf, oco.pres, ylim = c(1013, 0))


	## determine the separate level from STILT to OCO-2, using pressure
	# for model levels, keep OCO2 levels with zero AK above the max STILT level
	# for levels above model levels, use OCO2 prof
	min.par.pres <- min(uni.pres)
	upper.index <- oco.pres <  min.par.pres	  # return T/F
	lower.index <- oco.pres >= min.par.pres

	upper.oco.pres <- oco.pres[upper.index]
	lower.oco.pres <- oco.pres[lower.index]

	# combine LOWER STILT levels and UPPER OCO-2 levels
	combine.pres <- c(uni.pres, upper.oco.pres); names(combine.pres) <- NULL


	### -------------------- FOR a combined AK.norm prof ------------------ ###
	# interpolate AK for STILT levels if ak.wgt is TRUE;
	# OR assign all AK as 1 if ak.wgt is FALSE
	if (ak.wgt) {
		lower.ak.norm <- approx(x = oco.pres, y = oco.ak.norm, xout = uni.pres, rule = 2)$y

		# AK=0 for upper levels --> keep original OCO2 AK profiles, DW, 04/06/2017
		upper.ak.norm   <- oco.ak.norm[upper.index]
		combine.ak.norm <- c(lower.ak.norm, upper.ak.norm)

	} else combine.ak.norm <- rep(1, length(combine.pres))
			# assign 1 to all levels if not performing AK weighting, DW, 02/06/2018
	
	attributes(combine.ak.norm)$names <- NULL   # assign names


	### ------------------- FOR a combined a priori CO2 prof -------------- ###
	# interpolate for lower STILT levels
	# remain the upper OCO-2 apriori profiles for UPPER levels
	lower.ap   <- approx(x = oco.pres, y = oco.ap, xout = uni.pres, rule = 2)$y
	upper.ap   <- oco.ap[upper.index]
	combine.ap <- c(lower.ap, upper.ap)
	attributes(combine.ap)$names <- NULL   # assign names


	### ------------------- FOR a combined PWF prof ------------------------ ###
	# Interpolate and scale PWF for STILT levels; remain OCO2 PWF for OCO2 levels
	# MAY WANT TO TREAT THE VERY BOTTOM LEVEL DIFFERENTLY, SUM(PWF) should be 1

	## Step 1) Directly interpolate PWF for LOWER/STILT levels, need adjustments
	# later treat the bottom layer differently, only use PWF profiles above the
	# 1st layer to interpolate, no weird curve now, DW 09/20/2017

	# get interpolated PWF before scaling:
	lower.stilt.pwf.before <- approx(x = oco.pres[-1], y = oco.pwf[-1],
		                             xout = uni.pres[-1], rule = 2)$y

	# for debug--
	#plot(lower.stilt.pwf.before, uni.pres[-1], ylim = c(1013, 0))

	## Step 2) Calculate dP for STILT levels as well as LOWER/OCO-2 levels
	# diff in pres have one value less than the LEVELS
	lower.stilt.dp <- abs(diff(uni.pres))	        # for LOWER/STILT levels
	lower.oco.dp  <- abs(diff(lower.oco.pres))	  # for LOWER/OCO levels

	# DW 11/20/2016 --
	# !!! also, remember to calculate dp for the interface bwtween upper OCO2
	# levels & lower STILT levels, as dp ratio between two levels is always
	# assigned the upper level

	# get modeled and original OCO-2 dp ratios for the 'interface'
	intf.stilt.dp <- abs(diff(combine.pres))[length(uni.pres)]
	intf.oco.dp <- lower.oco.pres[length(lower.oco.pres)] - upper.oco.pres[1]

	## Step 3) Interpolate dp.oco.lower onto STILT levels using pres (EXCEPT the
	# bottom level, first element) + pres diff at LOWER OCO level

	# bug found, approx needs at least two non-NA values to interpolate
	# it occurred when we have low MAXAGL for bootstrapping
	lower.oco.dp.stilt <- approx(x = lower.oco.pres[-1], y = lower.oco.dp,
		                         xout = uni.pres[-1], rule = 2)$y

	## Step 4) As PWF is a function of pressure difference, larger dp, larger
	# air mass, CO2 or footprint should be weighted more, if ignoring the slight
	# variation in q (moisture). Note that STILT footprint and OCO2 are only
	# simulating/measuring the DRY AIR PROPERTIES

	# Thus, calculate the ratio of dp at STILT levels over interpolated OCO levels
	dp.ratio <- lower.stilt.dp/lower.oco.dp.stilt   	 # dp for LOWER/STILT levels

	# Bug fixed, DW 11/20/2016--
	# always assign new PWF to one upper level
	# !!!! we need to add one dp.ratio for the 1st OCO level above the interface
	intf.dp.ratio <- intf.stilt.dp/intf.oco.dp

	## Step 5) Scale interpolated PWF for STILT levels by multiplying 'dp' ratio
	# remember to leave alone the very bottom STILT level
	lower.pwf <- lower.stilt.pwf.before * dp.ratio
	attributes(lower.pwf)$names <- uni.pres[-1]

	## Step 6) Remain original OCO-2 PWF for UPPER OCO-2 levels
	upper.pwf <- oco.pwf[upper.index]

	# Bug fixed, DW, 11/20/2016--
	## Step 6.5) Replace 'pwf' for the lowest UPPER OCO2 level, calculate based on
	# 'dp' ratio at the interface 'intf.dp.ratio'
	intf.pwf  <- upper.pwf[1] * intf.dp.ratio
	upper.pwf <- c(intf.pwf, upper.pwf[-1])

	## Step 7) Calulate the PWF for the very bottom layer
	bottom.pwf  <- 1 - sum(lower.pwf) - sum(upper.pwf)
	combine.pwf <- c(bottom.pwf, lower.pwf, upper.pwf)
	attributes(combine.pwf)$names <- NULL


    ## ----- Combine all interpolated OCO-2 profiles and calculating AK_norm *PWF
	combine.prof <- data.frame(pres = combine.pres, pwf = combine.pwf,
		                       ak.norm = combine.ak.norm, ap = combine.ap) %>% 
				    mutate(ak.pwf = ak.norm * pwf) %>% 
					left_join(sel.p.rev, by = c('pres' = 'pres_rev')) %>% 
					rename(pres.stilt = pres.y)
					
	# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak,
	#AK*PWF and a priori contribution
	combine.prof$stiltTF <- F
	combine.prof[1:length(uni.pres), 'stiltTF'] <- T

	combine.prof  	# return weighting profiles
}

# end of subroutine



if (F) {

	oco.df <- data.frame(pres = oco.info$pres, ak.norm = oco.info$ak.norm, 
						 pwf = oco.info$pwf, ap = oco.info$ap)
	
	ak1 <- ggplot() + ylim(c(1013, 0)) + theme_bw() + 
		   geom_point(data = combine.prof, aes(ak.norm, pres), size = 0.3, shape = 21) +
		   geom_point(data = oco.df, aes(ak.norm, pres), color = 'red', shape = 21, size = 3)
		  
	ap1 <- ggplot() + ylim(c(1013, 0)) + theme_bw() + 
		   geom_point(data = combine.prof, aes(ap, pres), size = 0.3, shape = 21) +
		   geom_point(data = oco.df, aes(ap, pres), color = 'red', shape = 21, size = 3)

	p1 <- ggplot() + ylim(c(1013, 0)) + theme_bw() + 
		   geom_point(data = combine.prof, aes(pwf, pres), size = 0.3, shape = 21) +
		   geom_point(data = oco.df, aes(pwf, pres), color = 'orange', shape = 21, size = 3)
		  	  
}