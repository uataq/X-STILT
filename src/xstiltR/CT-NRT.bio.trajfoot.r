# subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# and calculate the dCO2 for each trajec for both biospheric signal)
# by Dien Wu, 04/11/2017
# fix a bug, (all footprint columns are zero), 05/11/2017
# add dmassTF for weighting footprint if violate mass conservation, DW, 10/20/2017

if(F){
  trajdat=orig.traj
  ident=sel.orig.outname
  trajdat=lower.traj
  ident=sel.err.outname[2]
}

ctnrt.bio.trajfoot<-function(trajdat=NULL,ident=NULL,dmassTF=T){

library(ncdf4)
library(reshape)
ct.version<-"2016-1"
if(substr(ident,1,10)>="2016x01x01")ct.version<-"2017"
ctpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/group_data/CT-NRT_v",ct.version,"/fluxes/optimized/",sep="")

max.index<-max(trajdat[,"index"])
full.index<-seq(1,max.index,1)

# mass conservation correction, grabbed from Trajecfoot()
if (dmassTF){

  #remove trajdaticles with too strong dmass violation
  ind<-unique(trajdat[trajdat[,"dmass"]>1E3|trajdat[,"dmass"]<1/1E3,"index"])

  if (length(ind) >= length(unique(trajdat[, "index"]))/2){
    message("Trajecvprm(): ", length(ind), ' of ', length(unique(trajdat[, "index"])), ' trajdaticles have mass defect; returning NA')
    return(NA)
  }
  trajdat<-trajdat[!trajdat[,"index"]%in%ind,]
  #trajdat<-cbind(trajdat,hour=trajdat[,"time"]/60)

  # get average dmass to "correct correction" (allow multiplication w/ dmass without changing total mass)
  # i.e. correction for average mass loss of trajdaticles, since they get attracted to areas of mass destruction
  mean.dmass<-rev(tapply(trajdat[, "dmass"], trajdat[, "time"], mean)) # this gives for each btime a mean dmass

  # DMM, To account for situations where mean.dmass is zero (mass violation total), need to avoid division by zero to avoid downstream problems.
  mean.dmass[which(mean.dmass == 0)] <- 0.00001

  # need to "merge" this with trajdat; can't use array since not all paticles are left at very large time
  nparleft<-rle(trajdat[, "time"])$length # number of times the same time is repeated
  mean.dmass<-rep(mean.dmass, nparleft) # long vector of mean dmass

  # need to link this info to each trajdaticles dmass: normalize individual dmass by mean.dmass
  trajdat[, "dmass"]<-trajdat[, "dmass"]/mean.dmass             # Dan Matross gets problems with that

  # also multiplied by dmass (accumulated weight of trajdaticles due to mass violation, normalized by average dmass to conserve total mass over time)
  trajdat[,"foot"]<-trajdat[,"foot"]*trajdat[,"dmass"]
}

# select traj where foot is non-zero
foot.flag<-trajdat[,"foot"]>1E-10
sel.trajdat<-trajdat[foot.flag,]

# initialize with NA for emissions
bio.flux<-rep(NA,nrow(trajdat))
trajdat<-cbind(trajdat,bio.flux)

if(length(sel.trajdat)!=0){
  #### NOW grab a CO2.bio fluxes for each selected trajec ###
  # match 3-houly footprint and 3-houly biosperhic fluxes
  # first locate the date and hour for each backwards hours of STILT footprint
  # use weekdayhr() to return the acutal date, weekdayhr<-function(yr,mon,day,hr,runtt,diffGMT=NA)

  recp.yr<-as.numeric(substr(ident,1,4));recp.mon<-as.numeric(substr(ident,6,7));recp.day<-as.numeric(substr(ident,9,10))
  recp.hr<-substr(ident,12,13)
  if(substr(recp.hr,1,1)=="x"){recp.hr<-10}else{recp.hr<-as.numeric(recp.hr)}
  back.min<-trajdat[,"time"]
  back.time<-weekdayhr(recp.yr, recp.mon, recp.day, recp.hr, runtt=back.min)
  traj.dd<-paste(back.time[,"yr"],formatC(back.time[,"mon"],width=2,flag=0),formatC(back.time[,"day"],width=2,flag=0),sep="") # increasing trend
  traj.hh<-back.time[,"hr"]

  sel.back.min<-sel.trajdat[,"time"]  # in mins
  sel.back.time<-weekdayhr(recp.yr, recp.mon, recp.day, recp.hr, runtt=sel.back.min)
  sel.traj.dd<-paste(sel.back.time[,"yr"],formatC(sel.back.time[,"mon"],width=2,flag=0),formatC(sel.back.time[,"day"],width=2,flag=0),sep="") # increasing trend
  sel.traj.hh<-sel.back.time[,"hr"]

  uni.day<-unique(traj.dd)

  # grab CT-NRT files for all backwards hours, and then put in a 3D fluxes array
  tmpfile<-list.files(path=ctpath, pattern=paste("CT-NRT.v",ct.version,".flux1x1.",uni.day[1],sep=""))
  tmp<-nc_open(paste(ctpath,tmpfile,sep=""))  # for fluxes
  if(ct.version=="2016-1"){ct.lat<-ncvar_get(tmp,"lat")-0.5;ct.lon<-ncvar_get(tmp,"lon")-0.5;ct.time<-ncvar_get(tmp,"date_components")}
  if(ct.version=="2017"){ct.lat<-ncvar_get(tmp,"latitude")-0.5;ct.lon<-ncvar_get(tmp,"longitude")-0.5;ct.time<-ncvar_get(tmp,"time_components")}
  rownames(ct.time)<-c("year", "month", "day", "hour", "minute", "second")
  ct.hh<-ct.time["hour",]-1

  # create a big 3D array for storing bio fluxes, only covers the regions in footprint [lon, lat, time]
  ctbio.all<-array(0, dim=c(length(ct.lon), length(ct.lat), length(uni.day), length(ct.hh)), dimnames=list(ct.lon, ct.lat, uni.day, ct.hh))
  #store.length<-rep(0, length(unique(foot.dd)))

  for (f in 1:length(uni.day)){

    # open the daily file just once
    ctfile<-list.files(path=ctpath, pattern=paste("CT-NRT.v",ct.version,".flux1x1.",uni.day[f],sep=""))
    ctdat<-nc_open(paste(ctpath,ctfile,sep=""))

    # grab all biosperhic fluxes
    ct.bio<-ncvar_get(ctdat, "bio_flux_opt")	# unit in mol/m2/s
    dimnames(ct.bio)<-list(ct.lon, ct.lat, ct.hh)
    ct.bio<-ct.bio*1E6	# convert to umol/m2/s

    # put into the big array for storing
    ctbio.all[,,f,]<-ct.bio  # umol/m2/s

    nc_close(ctdat)
  }# end for f

  # to save time, subset CO2 fluxes
  sel.ct.lat<-ct.lat[ct.lat >= floor(min(trajdat[,"lat"])) & ct.lat <= ceiling(max(trajdat[,"lat"]))] # N
  sel.ct.lon<-ct.lon[ct.lon >= floor(min(trajdat[,"lon"])) & ct.lon <= ceiling(max(trajdat[,"lon"]))] # E
  sel.bio<-ctbio.all[ct.lon >= floor(min(trajdat[,"lon"])) & ct.lon <= ceiling(max(trajdat[,"lon"])),ct.lat >= floor(min(trajdat[,"lat"])) & ct.lat <= ceiling(max(trajdat[,"lat"])),,]

  # after reading all bio fluxes, assign lon.index, lat.index, day.index and hour.index for each trajec/sel.trajec
  sel.lat.index<-match(trunc(sel.trajdat[,"lat"]),sel.ct.lat)
  sel.lon.index<-match(trunc(sel.trajdat[,"lon"]),sel.ct.lon)
  sel.day.index<-match(sel.traj.dd,uni.day)
  sel.hour.index<-findInterval(sel.traj.hh,ct.hh)

  lat.index<-match(trunc(trajdat[,"lat"]),sel.ct.lat)
  lon.index<-match(trunc(trajdat[,"lon"]),sel.ct.lon)
  day.index<-match(traj.dd,uni.day)
  hour.index<-findInterval(traj.hh,ct.hh)

  sel.combine.index<-paste(sel.lon.index, sel.lat.index, sel.day.index, sel.hour.index)
  combine.index<-paste(lon.index, lat.index, day.index, hour.index)

  # find the unique gridcells
  uni.index<-unique(sel.combine.index)
  uni.str<-matrix(as.numeric(unlist(strsplit(uni.index, " "))), ncol=4, byrow=TRUE)
  colnames(uni.str)<-c("uni.lon.index","uni.lat.index","uni.day.index","uni.hour.index")

  # now find whether all combine.index are within the uni.index
  # if yes, return the number in "uni.index"; if no, return NA
  match.index<-match(combine.index, uni.index)  # will have the same nrow as trajdat

  # in order to save time, only grab CO2 from ODIAC ONLY ONCE for all uni.index
  uni.bio<-NULL
  cat("ctnrt.trajfoot(): Grabbing CT-NRT bio fluxes. It takes time...\n")
  for(i in 1:length(uni.index)){
    if(length(dim(sel.bio))==3){
      uni.bio<-c(uni.bio,sel.bio[uni.str[i,"uni.lon.index"],uni.str[i,"uni.lat.index"],uni.str[i,"uni.hour.index"]])
    }else{
      uni.bio<-c(uni.bio,sel.bio[uni.str[i,"uni.lon.index"],uni.str[i,"uni.lat.index"],uni.str[i,"uni.day.index"],uni.str[i,"uni.hour.index"]])
    }
  }

  # now put fluxes in this large matrix
  # if match.index==NA, no emission assigned; if not, find the emission from vector "uni.co2"
  assign.flag<-!is.na(match.index)  # assign only when match.index is non-NA
  nonzero.match.index<-match.index[assign.flag]
  trajdat[assign.flag,"bio.flux"]<-uni.bio[nonzero.match.index]

} # end if

# if emission is NA or footprints are all zero, replace with zero
trajdat[is.na(trajdat[,"bio.flux"]),"bio.flux"]<-0

# NOW, foot and emiss have the same dim, multiple to get a dCO2
dco2<-trajdat[,"bio.flux"]*trajdat[,"foot"]
trajdat<-cbind(trajdat,dco2)

# also, compute total dCO2 for each traj over all backwards hours
traj.dco2.bio<-tapply(trajdat[,"dco2"], trajdat[,"index"], sum)
melt.dco2.bio<-melt(traj.dco2.bio);colnames(melt.dco2.bio)<-c("index","dco2.bio")

# fill with all index, if dmassTF==T, particles being removed, fill with 0
if(dmassTF){
  sel.index<-sort(unique(trajdat[,"index"]))
  rm.index<-full.index[full.index%in%sel.index==FALSE]
  rm.dco2.bio<-data.frame(index=rm.index,dco2.bio=rep(0,length(rm.index)))

  all.dco2.bio<-rbind(melt.dco2.bio,rm.dco2.bio)
  sort.dco2.bio<-all.dco2.bio[order(all.dco2.bio$index),]
}else{
  sort.dco2.bio<-melt.dco2.bio[order(melt.dco2.bio$index),]
}

return(sort.dco2.bio)

} # end of subroutine
