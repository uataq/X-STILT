### subroutine to get new AK PW profiles based on OCO-2's
# pressure weighting and averaging kernel, for combined levels
# OCO-2 only provides the PW, AK and a priori at 20 levels,
# use linear interpolation to "approx" values at given STILT releasing levels
# Dien Wu, 08/05/2016

# Bugs fixed, DW:
# fix 1, 11/20/2016, modify the interpolate of PWF from OCO2 to STILT
# fix 2, 11/28/2016, remove storeTF flag, as we always need .RData file stored
#                    for generating footprint using Trajectfoot()

# fix 3, 02/08/2017, add uneven vertical spacing
# fix 4, 04/06/2017, change AK from zero to original OCO2 AK above model level
# add 5, 04/19/2017, store interpolated AK PW apriori

# add 6, 04/20/2017, add control flags "pw.weight" for weighting trajec
# (1) default is to return profiles using AK and PW, ak.weight=T & pw.weight=T;
# (2) ak.weight=F & pw.weight=T for only pres wgt (when no need to use apriori);
# (3) ak.weight=T & pw.weight=F for only AK wgt
# (4) ak.weight=F & pw.weight=F no weighting (which returns the original trajec)

get.weight.funcv2 <- function(ident, recp.info, agl.info, trajdat, recp.grdhgt,
														  oco.ak.norm, oco.pw, oco.pres, oco.apriori,
														  ak.weight=TRUE, pw.weight=TRUE){

library(Hmisc)

if(ak.weight){
	cat("Turn on weighting OCO-2 simulation using averaging kernel...\n")
}else{
	cat("NO averaging kernel weighting, set ak to 1...\n")
}

#### ---------- DEALING WITH STILT TRAJ and STILT MULTIPLE HGTS ----------- ####
stilt.agl <- agl.info
stilt.nlevel <- recp.info$stilt.nlevel

#### -----------  CONVERT STILT altitude to pressure ---------------------- ####
# interpolate starting pressure based on starting hgts,
# by looking at press and altitude of particles at the first timestep backwards
# first select the particles at first delt backwards
cat("get.weight.funcv2(): inter/extrapolating OCO-2 profiles to STILT levels...\n")
min.time <- max(trajdat[, "time"])
sel.traj <- trajdat[trajdat[, "time"] == min.time, ]

# if interpolating using AGL
asl <- sel.traj[,"agl"] + sel.traj[,"grdht"]	# ASL from RData traj
sel.traj <- cbind(sel.traj, asl)

### do linear interpolation on the relationship between pressure and hgts FIRST,
# and output the starting pressure based on starting hgt

# rule=2 allows me to interpolate pressure beyond the range
#                   by using the data extreme, e.g., the surface pressure (z=0)
# ruel=1 return NA values to values beyond data range

# sel.traj$agl and sel.traj$pres for AGL heights and pressure of particles;
# also stilt.agl for stilt AGL hgts => we wanna output stilt.pres
stilt.asl <- stilt.agl + recp.grdhgt

# since there are values beyond the range which "approx" cannot predict, use "approxExtrap" function INSTEAD
stilt.pres.extrap <- approxExtrap(sel.traj[, "asl"], sel.traj[, "pres"], stilt.asl, method="linear", rule=2)
stilt.pres <- stilt.pres.extrap$y

#plot(sel.traj[,"asl"], sel.traj[,"pres"],ylim=c(0,1013))
#points(stilt.pres.extrap$x,stilt.pres,col="red")

#### ------------------------ DEALING WITH OCO-2 NOW ---------------------- ####
# since STILT display particles from surface-to-space, the opposite as OCO-2 product
# originally, profiles from OCO are from levels 1-20, from TOA to sfc.
# we need to reverse AK and PW profiles, as well as renaming attributes as pressure levels
# thus, flipped profiles (now from level 1 to 20) will be from sfc to TOA now...

oco.nlevel <- length(oco.pres)
oco.pres   <- oco.pres[length(oco.pres):1]	# flip pressure levels
attributes(oco.pres)$names <- attributes(oco.pres)$names[length(attributes(oco.pres)$names):1]

# now larger weighting as starting from the surface
oco.ak.norm <- oco.ak.norm[length(oco.ak.norm):1]
attributes(oco.ak.norm)$names <- oco.pres

# flip PW
oco.pw <- oco.pw[length(oco.pw):1]
attributes(oco.pw)$names <- oco.pres

# flip prior profile
oco.apriori <- oco.apriori[length(oco.apriori):1]
attributes(oco.apriori)$names <- oco.pres

# plot(oco.pw, oco.pres,ylim=c(1013,0))

# determine the separate level from STILT to OCO-2, using pressure
# for only STILT levels, keep OCO2 levels with zero AK above the max STILT level
# when pressure is smaller than the min of STILT pressure, use OCO2 profile
min.stilt.pres <- min(stilt.pres)
upper.index <- oco.pres < min.stilt.pres	# return T/F
lower.index <- oco.pres >= min.stilt.pres

### -------------------  FOR a combined pressure profile ------------------ ###
upper.oco.pres <- oco.pres[upper.index]
lower.oco.pres <- oco.pres[lower.index]

# combine LOWER STILT levels and UPPER OCO-2 levels
combine.pres <- c(stilt.pres, upper.oco.pres)
combine.nlevel <- length(combine.pres)

### -------------------- FOR a combined AK.norm profile ------------------- ###
# interpolate for LOWER STILT levels if ak.weight==TRUE;
# OR set all AK to 1 if ak.weight==FALSE

# DW, 04/06/2017, AK=0 for upper levels --> now keep original OCO2 AK profiles
if(ak.weight){
	lower.ak.norm <- approx(oco.pres, oco.ak.norm, stilt.pres, method="linear", rule=2)$y
	upper.ak.norm <- oco.ak.norm[upper.index]
	combine.ak.norm <- c(lower.ak.norm, upper.ak.norm)
}else{

	# DW, 02/06/2018, also assign 1 to upper levels if no AK weighting
	combine.ak.norm<-rep(1,combine.nlevel)
}

attributes(combine.ak.norm)$names <- combine.pres

### ------------------- FOR a combined a priori CO2 profile ---------------- ###
# interpolate for lower STILT levels
# remain the upper OCO-2 apriori profiles for UPPER levels
lower.apriori <- approx(oco.pres, oco.apriori, stilt.pres, method="linear", rule=2)$y
upper.apriori <- oco.apriori[upper.index]
combine.apriori <- c(lower.apriori, upper.apriori)
attributes(combine.apriori)$names <- combine.pres

### ------------------- FOR a combined PW profile ------------------------- ###
# Interpolate and scale PW for LOWER/STILT levels
# Remain PW for UPPER OCO-2 levels
# MAY WANT TO TREAT THE VERY BOTTOM LEVEL DIFFERENTLY, BY USING XSUM(PW)=1

# simply using dp/p_surface, do not turn on this
#combine.pw<-c(1-sum(abs(diff(combine.pres))/combine.pres[1]),abs(diff(combine.pres))/combine.pres[1])

# 1) directly interpolate the PW for LOWER/STILT levels, need adjustments later
# treat the bottom layer differently, only use PWF profiles above the first
#                       layer to interpolate, no weird curve now, DW 09/20/2017
# "lower.stilt.pw.before"--interpolated PW before scaling
lower.stilt.pw.before <- approx(oco.pres[-1], oco.pw[-1], stilt.pres[-1], method="linear", rule=2)$y
#plot(lower.stilt.pw.before, stilt.pres,ylim=c(1013,0))

# 2) calculate dP for STILT levels as well as LOWER/OCO-2 levels
# diff in pres have one value less than the LEVELS
#entire.oco.dp<-abs(diff(oco.pres))	      # for all original OCO-2 levels
lower.stilt.dp <- abs(diff(stilt.pres))	    # for LOWER/STILT levels
lower.oco.dp   <- abs(diff(lower.oco.pres))	  # for LOWER/OCO levels

# DW 11/20/2016 --
# !!! also, remember to calculate dp for the OCO/STILT interface
# because dp scaling factor between two levels is always assigned for the upper one level
interface.stilt.dp <- abs(diff(combine.pres))[stilt.nlevel]
interface.oco.dp   <- lower.oco.pres[length(lower.oco.pres)] - upper.oco.pres[1]

# 3) interpolate dp.oco.lower onto STILT levels
# using pres (EXCEPT the bottom level, first element) + pres diff at LOWER OCO level
# bug found DW, approx needs at least two non-NA values to interpolate
# bug occurs when we have small MAXAGL for bootstrapping
lower.oco.dp.stilt <- approx(lower.oco.pres[-1], lower.oco.dp, stilt.pres[-1], method="linear", rule=2)$y
#x11();plot(lower.oco.dp, lower.oco.pres[-1], ylim=c(1000,500), pch=19)
#points(lower.oco.dp.stilt, stilt.pres[-1], col="blue",pch=20)

# 4) since PW is a function of pressure difference,
# larger dp, larger air mass, should be weighted more
# if ignoring the slight variation in q (moisture), note that STILT footprint
#                 and OCO2 are only simulating/measuring the DRY AIR PROPERTIES
# thus, calculate the ratio of dp at LOWER/STILT levels
#                               over interpolated OCO dp for LOWER/STILT levels
dp.ratio <- lower.stilt.dp/lower.oco.dp.stilt

# DW 11/20/2016--
# always assign new pw to one upper level
# !!!! THUS, we need to add one more dp.ratio for the first OCO level above the interface
interface.dp.ratio <- interface.stilt.dp/interface.oco.dp

# 5) scale interpolated PW for STILT levels by multiplying the "dp ratio"
# remember to put aside the very bottom STILT level
lower.pw <- lower.stilt.pw.before * dp.ratio
attributes(lower.pw)$names <- stilt.pres[-1]

# 6) remain the PW for UPPER OCO-2 levels, from oco.pw and upper.oco.pres
upper.pw <- oco.pw[upper.index]	        # from the dividing level to TOA

# DW, 11/20/2016--
# 6.5) CHANGE the lowest UPPER OCO2 level using DP.RATIO at interface
interface.pw <- upper.pw[1] * interface.dp.ratio
upper.pw <- c(interface.pw, upper.pw[-1])

# 7) calulate the PW for the very bottom layer
bottom.pw  <- 1 - sum(lower.pw) - sum(upper.pw)
combine.pw <- c(bottom.pw, lower.pw, upper.pw)

attributes(combine.pw)$names <- combine.pres

combine.profile <- data.frame(combine.pres, combine.pw, combine.ak.norm, combine.apriori)
rownames(combine.profile) <- seq(1, combine.nlevel, 1)

# plotting
if(F){
	cex<-4.5
	png(filename="int_PWF.tiff",width=2000, height=2400)
	par(mar=c(12,12, 10, 5),mgp=c(8,2,0))
	plot(oco.pw,oco.pres,ylim=c(1000,0),xlim=c(0,0.053),pch=19,cex=cex*0.6,
	     main="Pressure Weighting Function (PWF)\nat OCO-2 and model levels",
			 xlab="Pressure Weighting Function (PWF)",ylab="Pressure [hpa]",
			 cex.main=cex,cex.lab=cex,cex.axis=cex)

	#lines(oco.pw,oco.pres,lty=2)
	points(lower.stilt.pw.before,stilt.pres[-1],col="darkorange",cex=cex,lwd=cex)
	points(oco.pw[1],oco.pres[1],col="darkorange",cex=cex,lwd=cex)
	points(upper.pw,upper.oco.pres,col="blue",cex=cex,lwd=cex)
	points(lower.pw,as.numeric(attributes(lower.pw)$names),col="red",cex=cex,lwd=cex)
	text(0.008, 475,"MAXAGL",col="red",cex=cex)
	abline(h=c(upper.oco.pres[1],upper.oco.pres[length(upper.oco.pres)]),col="blue",lty=2,lwd=cex)
	abline(h=c(stilt.pres[1],stilt.pres[length(stilt.pres)]),col="red",lty=2,lwd=cex)
	legend(0.0015, 25,c("Initial OCO-2 PW","OCO-2 PW at upper model levels",
	                   "Linearly interpolated PW at lower model levels",
										 "Scaled PW at lower model levels"),
										 col=c("black","blue","darkorange","red"),
										 text.col=c("black","blue","darkorange","red"),
										 bty="n", pch=c(19,1,1,1), cex=cex)
	dev.off()


	### plot AK.norm
	png(filename="int_AK.tiff", width=2000, height=2400)
	par(mar=c(12, 12, 10, 5), mgp=c(8,2,0))
	plot(oco.ak.norm, oco.pres, ylim=c(1000,0), xlim=c(0,1.2), pch=19, cex=cex*0.6,
	     main="Normalized averaging Kernel Profiles (AK)\nat OCO-2 and model levels",
			 xlab="Normalized Averaging Kernel (AK)", ylab="Pressure [hpa]",
			 cex.main=cex, cex.lab=cex, cex.axis=cex)

	lines(oco.ak.norm, oco.pres, lty=2, lwd=cex-1)
	points(upper.ak.norm, upper.oco.pres, col="blue", cex=cex, lwd=cex)
	points(lower.ak.norm, stilt.pres, col="red", cex=cex, lwd=cex)
	text(0.2, 475, "MAXAGL", col="red", cex=cex)
	abline(h=c(upper.oco.pres[1], upper.oco.pres[length(upper.oco.pres)]), col="blue", lty=2, lwd=cex)
	abline(h=c(stilt.pres[1], stilt.pres[length(stilt.pres)]), col="red", lty=2, lwd=cex)
	legend(0.001, 750, c("Initial OCO-2 AK","OCO-2 AK at upper model levels",
	                     "Linearly interpolated AK\nat lower model levels"),
											 col=c("black","blue","red"), text.col=c("black","blue","red"),
											 bty="n", lty=c(2,NA,NA), pch=c(19,1,1), cex=cex, lwd=c(cex,NA,NA))
	dev.off()
}


# calculating the AK*PW profile, and store back into "stilt.profile"
#stilt.profile$stilt.ak.pw <- stilt.profile$stilt.ak.norm * stilt.profile$stilt.pw
combine.profile$combine.ak.pw <- combine.profile$combine.ak.norm * combine.profile$combine.pw
# NOW ALL PROFILES CONTAIN--pressure, pressure weighting, normalized ak, AK*PW and a priori contribution

combine.profile$stiltTF <- F
combine.profile[1:length(stilt.pres), "stiltTF"] <- T
colnames(combine.profile) <- c("pres", "pw", "ak.norm", "oco2.prior", "ak.pw", "stiltTF")

return(combine.profile)	# return both weighting profiles and weighted trajec

}




# end of subroutine
