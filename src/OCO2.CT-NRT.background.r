# subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# use CO2_total mole fractions 3D fields to get the CO2 concentration for STILT particle endpoints
# for current purpose (Riyadh, Middle East), only use global 3x2 files, may want to add a nested NAM 1x1 for future sites
# ALSO, NEED TO WEIGHT THE ENDPOINT CONCENTRATION USING AVERAGING KERNEL AND PRESSURE WEIGHTING FUNCTIONS !!!
# Written by Dien Wu, 09/12/2016

# Previously, we only consider background CO2 for model levels using trajec endpoints, and AK aloft is ZERO
# but, what we should do is grabbing CT for upper atmos, and keep original OCO-2 AK aloft
# last fix on 04/06/2017, DW

# separate STILT level and CT-NRT levels, and output two XCO2, 04/24/2017, DW

# input variables needed--
#adjusted ak*pw profiles for STILT levels, .RData traj file, receptor time,

debugTF<-FALSE
if(debugTF){
	ident=ident
	trajdat=new.trajdat
	combine.profile=combine.profile
}

ctnrt.background<-function(ident=ident, trajdat=trajdat, combine.profile=combine.profile){

library(ncdf4)

# -------------------- PART 1 ---------------------------------- DEALING WITH ENDPOINTS OF TRAJ for model levels ------------------------------------------------------ #
# for grabbing the CO2 concentration at endpoints, we need the lon, lat, asl height and absolute time of STILT particles
cat("ctnrt.background.r(): Grabbing the ENDPOINTS of trajs...It takes time...");cat("\n")

ct.version<-"2016-1"
if(substr(ident,1,10)>="2016x01x01")ct.version<-"2017"
ctpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/group_data/CT-NRT_v",ct.version,"/molefractions/co2_total/",sep="")

# then, grab the information for particle that dropped off for each index, (before some of them drop off when crossing ZERO lon)
# change in 09/27/2016, using John's advice
# find out the min time according to index # as *ENDPOINTS* before DROPPING OFF, using TAPPLY commend
min.time<-tapply(trajdat[,"time"], trajdat[,"index"], min)	# return the min time [mins] back with attributes name as index

# NOW, grab the endtraj for those min.time, using STRING MATCHING METHOD
traj.time.index<-paste(trajdat[,"time"], trajdat[,"index"])
end.time.index<-paste(min.time, attributes(min.time)$dimnames[[1]])
row.index<-match(end.time.index,traj.time.index)	# which returns row number for each index
endtraj<-trajdat[row.index,]	# finally grab the endpoints

# checking the time...
#diff.time<-unique(endtraj[,"time"]-as.numeric(min.time))	# should be zero, meaning grabbing endpoint correctly

# THUS, endtraj data frame should always have numpar rows
# compute the ASL heights for STILT particles
asl<-endtraj[,"agl"] +endtraj[,"grdht"]
endtraj<-cbind(endtraj,asl)

# use CO2_total for STILT background
# locate which file to use, based on the receptor time and nhrs backwards
recp.year<-as.numeric(substr(ident,1,4))
recp.mon<-as.numeric(substr(ident,6,7))
recp.day<-as.numeric(substr(ident,9,10))
recp.hour<-as.numeric(substr(ident,12,13))

# compute the end yr, mon, day, hour for nhrs back to open a CT file
endp.time<-weekdayhr(yr=recp.year, mon=recp.mon, day=recp.day, hr=recp.hour, runtt=endtraj[,"time"])	# runtt in mins
uni.endp.day<-sort(unique(endp.time[,"day"]))	# for example for 26, 27, 28
uni.endp.year<-sort(unique(endp.time[,"yr"]))
uni.endp.mon<-sort(unique(endp.time[,"mon"]))


# since CT files contain 3D CO2 fields every 3 hours in UTC, assign a endpoint time index to each particle
# CT time--centered time on 1:30, 4:30, 7:30, 10:30, 13:30, 16:30, 19:30 and 22:30
# stands for time range from 00-03, 03-06, 06-09, 09-12, 12-15, 15-18, 18-21 and 21-24.

# Open all CT files into a list before assigning background CO2
all.ct.gph<-rep(list(NULL), length(uni.endp.day))
all.ct.co2<-all.ct.gph
ctfile<-paste(ctpath, "CT-NRT.v",ct.version,".molefrac_glb3x2_", formatC(uni.endp.year,width=2,flag=0), "-", formatC(uni.endp.mon,width=2,flag=0), "-", formatC(uni.endp.day,width=2,flag=0), ".nc", sep="")
cat(paste("Reading from file",ctfile,"\n"))

for(f in 1:length(uni.endp.day)){	# from smallest endp, meaning most backward time

		# find and open the CT file
		ctdat<-nc_open(ctfile[f])

		# list variable names for CO2_components,
		# bg (background, component due to initial condition), bio (terrestrial biopshere exchange), ff (fossil fuel burning), fires (direct fire emissions), ocean (air-sea exchange), latitude, longitude, level, time (since 2000-01-01 00:00:00 UTC)
		# list variable names for CO2_total,
		# air_mass, blh (PBL thickness), co2, decimal_date, gph (geoponential_height), orography (sfc geopotential), pressure, specific_humidity, temperature, time_components

		# ------------------------------------------------------ readin CT-NRT for getting molefractions ------------------------------------------------------ #
		# spatial resolution, 2 deg in N-S, 3 deg in E-W, move centered lat, lon to lower left
		ct.lat<-ncvar_get(ctdat, "latitude")-1
		ct.lon<-ncvar_get(ctdat, "longitude")-1.5
		ct.time<-ncvar_get(ctdat,"time_components")	# UTC time components
		rownames(ct.time)<-c("year", "month", "day", "hour", "minute", "second")
		# convert to YYYYMMDDHHmmss
		ct.timestr<-paste(ct.time["year",], formatC(ct.time["month",], width=2, flag=0), formatC(ct.time["day",], width=2, flag=0), formatC(ct.time["hour",], width=2, flag=0), formatC(ct.time["minute",], width=2, flag=0), formatC(ct.time["second",], width=2, flag=0), sep="")

		ct.level<-ncvar_get(ctdat,"level")	# 25 levels
		ct.boundary<-ncvar_get(ctdat,"boundary")	# 26 boundaries

		# grab geopotentail height, 26 BOUNDS FOR HGT
		ct.gph<-ncvar_get(ctdat,"gph")	# [LON, LAT, BOUND, TIME]
		dimnames(ct.gph)<-list(ct.lon, ct.lat, ct.boundary, ct.timestr)

		# grab CO2 fields, 25 LEVELS FOR CO2
		ct.co2<-ncvar_get(ctdat,"co2")	# [LON, LAT, LEVEL, TIME]
		dimnames(ct.co2)<-list(ct.lon, ct.lat, ct.level, ct.timestr)

		# combine all arrays into list
		all.ct.gph[[f]]<-ct.gph
		all.ct.co2[[f]]<-ct.co2

		nc_close(ctdat)

} # end reading in ct files
names(all.ct.gph)<-uni.endp.day
names(all.ct.co2)<-uni.endp.day


# use lat, lon, asl to locate the CT grid (all CT files have same lat lon grid)
lat.index<-trunc((endtraj[,"lat"]-ct.lat[1])/2+1)
lon.index<-trunc((endtraj[,"lon"]-ct.lon[1])/3+1)
# CHECKING... 2 degs in lat, 3 degs in lon
#diff.lat<-endtraj[,"lat"]-ct.lat[lat.index]
#diff.lon<-endtraj[,"lon"]-ct.lon[lon.index]

# locate which day of this endpoint locates
last.day<-endp.time[,"day"]

# find which CT hour range that the last STILT hour (endp.time[,"hr"]) falls in
# CT hours are 0130 (for 0000-0300), 0430, 0730, 1030, 1330, 1630, 1930, and 2230
# calculate the CT hour index, and grab the CO2 fields at this hour
last.hour<-endp.time[,"hr"]
ct.hour.index<-trunc(last.hour/3+1)	# which ct time to look for CO2

# -------------------------------------- DEALING WITH CT 26 BOUNDARY AND 25 LEVEL for hgt.index ----------------------------------- #
# assign "hgt.index" to each particles, referring which CT levels to grab, comparing ASL from STILT and GPH from CT
# our goal is to relate CO2 concentration with the hgts, not levels
#Determines the location, i.e., index of the (first) minimum or maximum of a numeric (or logical) vector, using "which.min" function
cat("ctnrt.background.r(): Grabbing background CO2 for each particle");cat("\n")

# initial ct.bg.co2 with NA
ct.bg.co2<-rep(NA, nrow(endtraj))
endtraj<-cbind(endtraj, ct.bg.co2)

for (p in 1:nrow(endtraj)){

		# simplify on 04/12/2017, DW
		bound.index<-findInterval(endtraj[p,"asl"],ct.gph[lon.index[p],lat.index[p],,ct.hour.index[p]])

		if(f){
			# for each particles, compare the STILT ASL and CT GPH and return a index for which CT level that this particle falls in
			# first CT.GPH must be smaller than STILT ASL, thus find the LOWER BOUNDARY index for which diff.hgt is the largest negative number
			diff.hgt<-ct.gph[lon.index[p],lat.index[p],,ct.hour.index[p]] - endtraj[p,"asl"]

			# if all CT GPHs are above the STILT ASL, just assign CO2 at the very bottom level to this particle
			neg.diff.hgt.bound<-as.numeric(attributes(diff.hgt[diff.hgt<=0])$names)
			if(length(neg.diff.hgt.bound)==0){
				bound.index<-1
			}else{
				# if there are CT GPHs <= STILT ASL, then--
				bound.index<-max(as.numeric(attributes(diff.hgt[diff.hgt<=0])$names))	# lower boundary index
			}
		}

		# convert to level index, e.g., level 1 CO2 is the concentration between BOUND 1 and BOUND 2
		# thus, level index == lower bound index
		# grab co2 according to day.index, lon.index, lat.index, bound.index, and hour index
		endtraj[p,"ct.bg.co2"]<-all.ct.co2[names(all.ct.gph)==last.day[p]][[1]][lon.index[p],lat.index[p],bound.index,ct.hour.index[p]]

} # end assigning CO2 background

# for checking whether CT-NRT is reasonable
#plot(endtraj[,"ct.bg.co2"], endtraj[,"agl"])
#library(scatterplot3d)
#cex<-1.2
#trajplot<-scatterplot3d(x=endtraj[,"lon"],y=endtraj[,"lat"],z=endtraj[,"agl"],pch=20,angle=60,color="red", xlim=c(1, 51), ylim=c(7, 29), zlim=c(0, 8000), xlab="Longitude",ylab="Latitude",zlab="Altitude [m AGL]",box=FALSE, cex.main=cex, cex.axis=cex, cex.lab=cex, cex=cex)

# checking--whether ct.bg.co2 all have values, if zero, good
#check<-length(which(is.na(endtraj[,"ct.bg.co2"])==TRUE))

# take the mean of ensemble at each STILT releasing levels
co2.bg.stilt<-tapply(endtraj[,"ct.bg.co2"], endtraj[,"level"], mean)
# because ct.bg.co2 is not normal distributed, try 50th percentile, instead of mean, DW, 11/16/2017
#co2.bg.stilt<-tapply(endtraj[,"ct.bg.co2"], endtraj[,"level"], median)

# --------------------- PART 2 --------------------------------- DEALING WITH ENDPOINTS OF TRAJ for model levels ------------------------------------------------------ #
# grab CT pressure and interpolate CT at receptors for high up
cat("ctnrt.background.r(): Grabbing CT CO2 above max model level");cat("\n")
id<-oco2.info[[2]]$find.id
recp.yr<-substr(id,1,4);recp.mon<-substr(id,5,6);recp.day<-substr(id,7,8);recp.hr<-substr(id,9,10)
recp.lat<-oco2.info[[2]]$find.lat;recp.lon<-oco2.info[[2]]$find.lon

# find and open the CT file
#ctpath<-"/uufs/chpc.utah.edu/common/home/lin-group1/group_data/CT-NRT_v2016-1/molefractions/co2_total/"
ctfile<-paste(ctpath, "CT-NRT.v",ct.version,".molefrac_glb3x2_", formatC(recp.yr,width=2,flag=0), "-", formatC(recp.mon,width=2,flag=0), "-", formatC(recp.day,width=2,flag=0), ".nc", sep="")
ctdat<-nc_open(ctfile)

# find the closest lat lon hour grid (to receptor) and grab the CT profiles
# and then interpolate onto combine levels
# spatial resolution, 2 deg in N-S, 3 deg in E-W, move centered lat, lon to lower left
ct.lat<-ncvar_get(ctdat, "latitude")-1
ct.lon<-ncvar_get(ctdat, "longitude")-1.5
ct.time<-ncvar_get(ctdat,"time_components")	# UTC time components
rownames(ct.time)<-c("year", "month", "day", "hour", "minute", "second")
ct.time["hour",]<-ct.time["hour",]-1;ct.time["minute",]<-ct.time["minute",]-30  # move time from center time to the beginning of each time period
ct.hour<-ct.time["hour",]	# starting hours

# compare recptor hour with CT hour, and find the time index
# similar for lat, lon
hr.index<-findInterval(as.numeric(recp.hr),ct.time["hour",])
lat.index<-findInterval(recp.lat,ct.lat)
lon.index<-findInterval(recp.lon,ct.lon)

# grab CO2 fields, 25 LEVELS FOR CO2
ct.co2<-ncvar_get(ctdat,"co2")	# [LON, LAT, LEVEL, TIME]
ct.level<-ncvar_get(ctdat,"level")	# 25 levels
dimnames(ct.co2)<-list(ct.lon, ct.lat, ct.level, ct.hour)

# grab pressure level instead of geopotentail height, at 26 boundaries
ct.pres<-ncvar_get(ctdat, "pressure")/100 	# convert Pa, to hPa, [LON, LAT, BOUND, TIME]
ct.boundary<-ncvar_get(ctdat,"boundary")	# 26 boundaries
dimnames(ct.pres)<-list(ct.lon, ct.lat, ct.boundary, ct.hour)

# NOW, grab the CT CO2 profiles at the recptor
sel.co2<-ct.co2[lon.index, lat.index,, hr.index]  # 25 CT levels
sel.pres<-ct.pres[lon.index, lat.index,, hr.index]  # pressure in hPa

if(F){
	# compute the Pressure height relation
	# compute geopotential height for combine levels
	g0<-9.80665;re<-6371*1E3  # in meters
	standard.asl<-seq(0,81000, 100)
	gh<-g0*(re/(re+standard.asl))
	standard.pres<-1013.25*exp(-gh*0.0289*standard.asl/8.31445/288.15)
	#standard.pres<-oco2.info[[2]]$psfc*exp(-gh*0.0289*standard.asl/8.31445/288.15)

	# due to exponential relation between pressure and asl, apply log to asl and exp() after the extrapolating
	asl.extrap1<-approxExtrap(pres.hgt$stilt.pres,log(pres.hgt$stilt.asl),combine.profile$combine.pres, method="linear", rule=2)
	asl.extrap1$y<-exp(asl.extrap1$y) # x for pressure, y for asl height
	asl.extrap2<-approxExtrap(standard.pres,standard.asl,combine.profile$combine.pres[38:47], method="linear", rule=2)
	upper.pres.asl<-data.frame(asl.extrap2$x, asl.extrap2$y);colnames(upper.pres.asl)<-colnames(pres.hgt)
	combine.pres.asl<-as.data.frame(rbind(pres.hgt, upper.pres.asl))

	plot(pres.hgt$stilt.asl,pres.hgt$stilt.pres,ylim=c(1000,0),xlim=c(0,50000),pch=19, xlab="ASL height [m]", ylab="Pressure")
	lines(standard.asl, standard.pres,lty=2,col="red",lwd=2)
	points(asl.extrap1$y[38:47],asl.extrap1$x[38:47],col="blue",pch=2,type="o",lwd=2,lty=2)
	points(asl.extrap2$y,asl.extrap2$x,col="red",lwd=2)
	legend("bottomright",c("STILT pres-hgt relation for Modeled Levels","Direct Exponential Extrapolation","Using Simple Press-Hgt Equation"),col=c("black","blue","red"),pch=c(19,2,1),lwd=c(NA,2,2),lty=c(NA,2,2))

	# interpolating CT CO2 profiles onto combine levels
	ctco2.extrap1<-approxExtrap(sel.gph[-26], sel.co2, combine.pres.asl$stilt.asl, method="linear", rule=2)
	layer.index<-findInterval(combine.pres.asl$stilt.asl, sel.gph)
	ctco2.extrap2<-sel.co2[layer.index]

	x11();plot(sel.co2,sel.gph[-26], ylim=c(0, 80000),xlim=c(385,405),pch=19,type="o",lty=2)
	points(ctco2.extrap1$y,ctco2.extrap1$x,col="blue",pch=2,type="o",lty=2)
	points(ctco2.extrap2,sel.gph[layer.index],col="red",type="o",lty=2)
}

# interpolating CT CO2 profiles onto combine levels, only above STILT levels
model.pres<-combine.profile$pres
upper.pres<-model.pres[model.pres < model.pres[max(which(combine.profile$stiltTF==T))]]

# assigning CT CO2 values to each model pressure, go through each
co2.bg.upper<-NULL
for (l in 1:length(upper.pres)){

		# locate which layer that model level falls in
		diff.pres<-sel.pres - upper.pres[l]

		# differ from GPH, pressure is opposite
		pres.bound<-as.numeric(attributes(diff.pres[diff.pres>= 0])$names)
		if(length(pres.bound)==0){bound.index<-1}else{bound.index<-max(pres.bound)}	# lower boundary index
		co2.bg.upper<-c(co2.bg.upper,sel.co2[names(sel.co2)==bound.index])

} # end assigning CO2 background for upper levels


# --------------------- PART 3 --------------------------------- COMBINING UPPER AND LOWER bg CO2 ------------------------------------------------------ #
co2.bg.combine<-c(co2.bg.stilt, co2.bg.upper)
combine.profile$ctnrt.bg<-co2.bg.combine

# apply ak pw weighting for each level mean
xco2.bg<-sum(combine.profile$ak.pw*combine.profile$ctnrt.bg, na.rm=T)

return(list(combine.profile,xco2.bg))

} # end of subroutine
