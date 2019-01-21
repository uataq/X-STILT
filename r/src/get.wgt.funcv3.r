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

get.wgt.funcv3 <- function(output, oco2.info, ak.wgt = T, pwf.wgt = T){

	if (ak.wgt) {
		cat("Turn on weighting OCO-2 simulation using averaging kernel...\n")
	} else {
		cat("NO averaging kernel weighting, set ak to 1...\n")
	}

    # grab trajectory and receptor info
    receptor  <- output$receptor
	trajdat   <- output$particle
    r_zagl    <- receptor$zagl   # vector form of column release levels, mAGLs
	zsfc      <- receptor$zsfc   # modeled ground heights from get.ground.height()
	recp.nlev <- length(r_zagl)  # numbers of levels for X-STILT

    # grab OCO-2 profiles
	# e.g., press weighting function, pressures, normalized AK, a priori
	oco2.ap       <- oco2.info$ap
	oco2.pwf      <- oco2.info$pwf
	oco2.pres     <- oco2.info$pres
    oco2.ak.norm  <- oco2.info$ak.norm


	#### ------ CONVERT STILT release levels from heights to pressures ----- ####
	# interpolate starting pressure based on starting hgts, by looking at press
	# and altitude of particles at the first timestep back

	# first select the particles at first time step back
	min.time <- min(abs(trajdat$time)) * sign(trajdat$time[1])  # MIN time in mins
	sel.traj <- trajdat[trajdat$time == min.time, ] # trajec near receptor

	# calculate mASL, add ground height (zsfc) to mAGL
	sel.traj$zasl <- sel.traj$zagl + sel.traj$zsfc	# ASL near receptor

	# calculate mASL of each release level (mAGL 'r_zagl' + sfc height 'zsfc')
	recp.zasl <- r_zagl + zsfc

	### interpolate pressures for each release levels given their values in mAGLs
	# rule=2 allows us to interpolate pressure beyond the range by using the data
	#        extreme, e.g., the surface pressure (z=0)
	# ruel=1 return NA values to values beyond data range
	recp.pres <- Hmisc::approxExtrap(x = sel.traj$zasl, y = sel.traj$pres,
		                               xout = recp.zasl, rule = 2)$y

	# for debug--
	#plot(sel.traj$zasl, sel.traj$pres, ylim = c(1013, 0))
	#points(recp.zasl, recp.pres, col="red")


	#### ------------------------ DEALING WITH OCO-2 NOW -------------------- ####
	# since STILT display particles from surface-to-space, the opposite as OCO-2
	# product originally, profiles from OCO are from levels 1-20, from TOA to sfc.
	# we need to reverse AK and PWF profiles, as well as renaming attributes as
	# pressure levels. Thus, flipped profiles (now from level 1 to 20) will be
	# from sfc to TOA now...

	oco2.nlev <- length(oco2.pres)
	oco2.pres <- oco2.pres[length(oco2.pres):1]	# flip pressure levels
	attributes(oco2.pres)$names <-
	            attributes(oco2.pres)$names[length(attributes(oco2.pres)$names):1]

	# flip 20 levels (--> from sfc to TOA) and assign names
	oco2.ak.norm <- oco2.ak.norm[length(oco2.ak.norm):1]
	oco2.ap  <- oco2.ap[length(oco2.ap):1]
	oco2.pwf <- oco2.pwf[length(oco2.pwf):1]

    attributes(oco2.ap)$names      <- oco2.pres
    attributes(oco2.pwf)$names     <- oco2.pres
	attributes(oco2.ak.norm)$names <- oco2.pres

	# for debug--
	# plot(oco2.pwf, oco2.pres, ylim = c(1013, 0))

	## determine the separate level from STILT to OCO-2, using pressure
	# for model levels, keep OCO2 levels with zero AK above the max STILT level
	# for levels above model levels, use OCO2 prof
	min.recp.pres <- min(recp.pres)
	upper.index   <- oco2.pres <  min.recp.pres	  # return T/F
	lower.index   <- oco2.pres >= min.recp.pres


	### -------------------  FOR a combined pressure prof ----------------- ###
	upper.oco2.pres <- oco2.pres[upper.index]
	lower.oco2.pres <- oco2.pres[lower.index]

	# combine LOWER STILT levels and UPPER OCO-2 levels
	combine.pres <- c(recp.pres, upper.oco2.pres)
	combine.nlev <- length(combine.pres)


	### -------------------- FOR a combined AK.norm prof ------------------ ###
	# interpolate AK for STILT levels if ak.wgt is TRUE;
	# OR assign all AK as 1 if ak.wgt is FALSE
	if (ak.wgt) {
	  lower.ak.norm <- approx(x = oco2.pres, y = oco2.ak.norm, xout = recp.pres,
			                      rule = 2)$y

	  # AK=0 for upper levels --> keep original OCO2 AK profiles, DW, 04/06/2017
	  upper.ak.norm   <- oco2.ak.norm[upper.index]
	  combine.ak.norm <- c(lower.ak.norm, upper.ak.norm)

	} else {
	  # assign 1 to all levels if no AK weighting, DW, 02/06/2018
	  combine.ak.norm <- rep(1, combine.nlev)
	}
	attributes(combine.ak.norm)$names <- combine.pres   # assign names


	### ------------------- FOR a combined a priori CO2 prof -------------- ###
	# interpolate for lower STILT levels
	# remain the upper OCO-2 apriori profiles for UPPER levels
	lower.ap   <- approx(x = oco2.pres, y = oco2.ap, xout = recp.pres, rule = 2)$y
	upper.ap   <- oco2.ap[upper.index]
	combine.ap <- c(lower.ap, upper.ap)
	attributes(combine.ap)$names <- combine.pres


	### ------------------- FOR a combined PWF prof ------------------------ ###
	# Interpolate and scale PWF for STILT levels; remain OCO2 PWF for OCO2 levels
	# MAY WANT TO TREAT THE VERY BOTTOM LEVEL DIFFERENTLY, SUM(PWF) should be 1

	## Step 1) Directly interpolate PWF for LOWER/STILT levels, need adjustments
	# later treat the bottom layer differently, only use PWF profiles above the
	# 1st layer to interpolate, no weird curve now, DW 09/20/2017

	# get interpolated PWF before scaling:
	lower.stilt.pwf.before <- approx(x = oco2.pres[-1], y = oco2.pwf[-1],
		                             xout = recp.pres[-1], rule = 2)$y

	# for debug--
	#plot(lower.stilt.pwf.before, recp.pres, ylim = c(1013, 0))

	## Step 2) Calculate dP for STILT levels as well as LOWER/OCO-2 levels
	# diff in pres have one value less than the LEVELS
	lower.stilt.dp <- abs(diff(recp.pres))	        # for LOWER/STILT levels
	lower.oco2.dp  <- abs(diff(lower.oco2.pres))	  # for LOWER/OCO levels

	# DW 11/20/2016 --
	# !!! also, remember to calculate dp for the interface bwtween upper OCO2
	# levels & lower STILT levels, as dp ratio between two levels is always
	# assigned the upper level

	# get modeled and original OCO-2 dp ratios for the 'interface'
	intf.stilt.dp <- abs(diff(combine.pres))[recp.nlev]
	intf.oco2.dp  <- lower.oco2.pres[length(lower.oco2.pres)] - upper.oco2.pres[1]

	## Step 3) Interpolate dp.oco2.lower onto STILT levels using pres (EXCEPT the
	# bottom level, first element) + pres diff at LOWER OCO level

	# bug found, approx needs at least two non-NA values to interpolate
	# it occurred when we have low MAXAGL for bootstrapping
	lower.oco2.dp.stilt <- approx(x = lower.oco2.pres[-1], y = lower.oco2.dp,
		                          xout = recp.pres[-1], rule = 2)$y

	## Step 4) As PWF is a function of pressure difference, larger dp, larger
	# air mass, CO2 or footprint should be weighted more, if ignoring the slight
	# variation in q (moisture). Note that STILT footprint and OCO2 are only
	# simulating/measuring the DRY AIR PROPERTIES

	# Thus, calculate the ratio of dp at STILT levels over interpolated OCO levels
	dp.ratio <- lower.stilt.dp/lower.oco2.dp.stilt   	 # dp for LOWER/STILT levels

	# Bug fixed, DW 11/20/2016--
	# always assign new PWF to one upper level
	# !!!! we need to add one dp.ratio for the 1st OCO level above the interface
	intf.dp.ratio <- intf.stilt.dp/intf.oco2.dp

	## Step 5) Scale interpolated PWF for STILT levels by multiplying 'dp' ratio
	# remember to leave alone the very bottom STILT level
	lower.pwf <- lower.stilt.pwf.before * dp.ratio
	attributes(lower.pwf)$names <- recp.pres[-1]

	## Step 6) Remain original OCO-2 PWF for UPPER OCO-2 levels
	upper.pwf <- oco2.pwf[upper.index]

	# Bug fixed, DW, 11/20/2016--
	## Step 6.5) Replace 'pwf' for the lowest UPPER OCO2 level, calculate based on
	# 'dp' ratio at the interface 'intf.dp.ratio'
	intf.pwf  <- upper.pwf[1] * intf.dp.ratio
	upper.pwf <- c(intf.pwf, upper.pwf[-1])

	## Step 7) Calulate the PWF for the very bottom layer
	bottom.pwf  <- 1 - sum(lower.pwf) - sum(upper.pwf)
	combine.pwf <- c(bottom.pwf, lower.pwf, upper.pwf)
	attributes(combine.pwf)$names <- combine.pres

    ## Step 8) Combine all interpolated OCO-2 profiles
	combine.prof <- data.frame(pres = combine.pres, pwf = combine.pwf,
		                       ak.norm = combine.ak.norm, ap = combine.ap)
	rownames(combine.prof) <- seq(1, combine.nlev, 1)

	# calculating the AK*PWF prof, and store back into "stilt.prof"
	combine.prof$ak.pwf <- combine.prof$ak.norm * combine.prof$pwf

	# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak,
	#AK*PWF and a priori contribution
	combine.prof$stiltTF <- F
	combine.prof[1:length(recp.pres), "stiltTF"] <- T

	combine.prof  	# return weighting profiles
}

# end of subroutine
