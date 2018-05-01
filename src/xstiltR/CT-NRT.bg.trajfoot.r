# subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# and derive background CO2 for each trajec (with/without transport errors)
# by Dien Wu, 04/12/2017

# fix 1, 04/19/2017, DW, add AK & PW weighting to background CO2
# fix 2, remove fix1, 04/20/2017, DW, not to weight background CO2 for now

if(F){
  trajdat=upper.traj
  ident=sel.err.outname[2]
}

ctnrt.bg.trajfoot<-function(trajdat=NULL,ident=NULL){

library(ncdf4)
ct.version<-"2016-1"
if(substr(ident,1,10)>="2016x01x01")ct.version<-"2017"
ctpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/group_data/CT-NRT_v",ct.version,"/molefractions/co2_total/",sep="")

#### NOW grab a CO2.bio fluxes for each selected trajec ###
# match 3-houly footprint and 3-houly biosperhic fluxes
# first locate the date and hour for each backwards hours of STILT footprint
recp.yr<-as.numeric(substr(ident,1,4));recp.mon<-as.numeric(substr(ident,6,7));recp.day<-as.numeric(substr(ident,9,10))
recp.hr<-substr(ident,12,13)
if(substr(recp.hr,1,1)=="x"){recp.hr<-10}else{recp.hr<-as.numeric(recp.hr)}
endpoint.min<-tapply(trajdat[,"time"],trajdat[,"index"],min)  # in mins

# NOW, grab the endtraj for those min.time, using STRING MATCHING METHOD
traj.time.index<-paste(trajdat[,"time"], trajdat[,"index"])
end.time.index<-paste(endpoint.min, attributes(endpoint.min)$dimnames[[1]])
row.index<-match(end.time.index,traj.time.index)	# which returns row number for each index
endtraj<-trajdat[row.index,]	# finally grab the endpoints

# open all three-day-back CT-NRT files
# use weekdayhr() to return the acutal date, weekdayhr<-function(yr,mon,day,hr,runtt,diffGMT=NA)
endtraj.time<-weekdayhr(recp.yr, recp.mon, recp.day, recp.hr, runtt=endtraj[,"time"])
traj.hh<-endtraj.time[,"hr"]
traj.dd<-paste(endtraj.time[,"yr"],formatC(endtraj.time[,"mon"],width=2,flag=0),formatC(endtraj.time[,"day"],width=2,flag=0),sep="")
uni.day<-sort(unique(traj.dd))  # increasing trend

# grab CT-NRT files for all backwards hours, and then put in a 3D fluxes array
ctfile<-paste(ctpath, "CT-NRT.v",ct.version,".molefrac_glb3x2_", substr(uni.day,1,4), "-", substr(uni.day,5,6), "-", substr(uni.day,7,8), ".nc", sep="")
tmp<-nc_open(ctfile[1])  # for fluxes

# grab variables
ct.lat<-ncvar_get(tmp, "latitude")-1
ct.lon<-ncvar_get(tmp, "longitude")-1.5
ct.time<-ncvar_get(tmp,"time_components")	# UTC time components
ct.level<-ncvar_get(tmp, "level")	# 25 levels
ct.pres<-ncvar_get(tmp, "pressure")/100 # 26 boundaries for pressure
ct.bound<-ncvar_get(tmp,"boundary")	# 26 boundaries
rownames(ct.time)<-c("year", "month", "day", "hour", "minute", "second")
ct.hh<-ct.time["hour",]-1 # every 3 hours

# create a big 5D array for storing CO2 concentration, [lon, lat, 25layer, day, hour]
ctco2.all<-array(0, dim=c(length(ct.lon), length(ct.lat), length(ct.level), length(uni.day), length(ct.hh)), dimnames=list(ct.lon, ct.lat, ct.level, uni.day, ct.hh))
ctgph.all<-array(0, dim=c(length(ct.lon), length(ct.lat), length(ct.bound), length(uni.day), length(ct.hh)), dimnames=list(ct.lon, ct.lat, ct.bound, uni.day, ct.hh))

for (f in 1:length(uni.day)){

  # open the daily file just once
  ctdat<-nc_open(ctfile[f])

	# grab CO2 fields, 25 LEVELS FOR CO2
	ct.co2<-ncvar_get(ctdat,"co2")	# [LON, LAT, LEVEL, TIME]
	dimnames(ct.co2)<-list(ct.lon, ct.lat, ct.level, ct.hh)

	# grab geopotentail height, 26 BOUNDS FOR HGT
	ct.gph<-ncvar_get(ctdat,"gph")	# [LON, LAT, BOUND, TIME]
	dimnames(ct.gph)<-list(ct.lon, ct.lat, ct.bound, ct.hh)

  # put into the big array for storing
  ctco2.all[,,,f,]<-ct.co2  # ppm
  ctgph.all[,,,f,]<-ct.gph  # in meter

  nc_close(ctdat)
}# end for f

# to save time, subset CO2 concentration and gph
cut.lat<- ct.lat >= floor(min(endtraj[,"lat"])-1) & ct.lat <= ceiling(max(endtraj[,"lat"]))
cut.lon<- ct.lon >= floor(min(endtraj[,"lon"])-2) & ct.lon <= ceiling(max(endtraj[,"lon"]))

sel.ct.lat<-ct.lat[cut.lat] # N
sel.ct.lon<-ct.lon[cut.lon] # E
sel.ct.co2<-ctco2.all[cut.lon,cut.lat,,,]
sel.ct.gph<-ctgph.all[cut.lon,cut.lat,,,]

# after reading all bio fluxes, assign lon.index, lat.index, day.index and hour.index for each selected trajec
lat.index<-findInterval(trunc(endtraj[,"lat"]),sel.ct.lat) # should not be zero
lon.index<-findInterval(trunc(endtraj[,"lon"]),sel.ct.lon)
day.index<-match(traj.dd,uni.day)
hour.index<-findInterval(traj.hh,ct.hh)

bg.co2<-rep(NA,nrow(endtraj))
endtraj<-cbind(endtraj,bg.co2)

for(e in 1:nrow(endtraj)){

  asl<-endtraj[e,"agl"]+endtraj[e,"grdht"]

  ### !!! if all endtraj times are the same, previous 5D array will become 4D, which causes a dimension error
  if(length(dim(sel.ct.co2))==4)ct.hgt<-sel.ct.gph[lon.index[e], lat.index[e], , hour.index[e]]
  if(length(dim(sel.ct.co2))==5)ct.hgt<-sel.ct.gph[lon.index[e], lat.index[e], , day.index[e], hour.index[e]]

  hgt.index<-findInterval(asl,ct.hgt)
  if(hgt.index==0)hgt.index<-1  # if below the lowest CT levels, use the first level

  # assign CO2 concentration from CT-NRT to endtraj
  if(length(dim(sel.ct.co2))==4)endtraj[e,"bg.co2"]<-sel.ct.co2[lon.index[e], lat.index[e], hgt.index, hour.index[e]]
  if(length(dim(sel.ct.co2))==5)endtraj[e,"bg.co2"]<-sel.ct.co2[lon.index[e], lat.index[e], hgt.index, day.index[e], hour.index[e]]

}

# fix 2, comment out the weighting ..., DW, 04/20/2017
# further weighted through AK and PW profiles
#wgt.bg.co2<-rep(NA,nrow(endtraj))
#endtraj<-cbind(endtraj,wgt.bg.co2)
#for (n in 1:max(endtraj[,"level"])){
#  endtraj[endtraj[,"level"]==n,"wgt.bg.co2"]<-endtraj[endtraj[,"level"]==n,"bg.co2"]*prof[n,"ak.pw"]
#}

# also, compute total dCO2 for each traj over all backwards hours
edp.co2<-data.frame(index=endtraj[,"index"],edp=endtraj[,"bg.co2"])

return(edp.co2)

} # end of subroutine
