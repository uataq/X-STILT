# subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# use CO2 optmized fluxes to get the biosperhic co2 exchange
# match weighted footprint with biosperhic co2 exchange
# Written by Dien Wu, 02/10/2017

#sourcepath<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/stiltR/";source(paste(sourcepath,"sourceall.r",sep=""))

# input variables needed--
# ident, footprint matrix, receptor time,

debugTF<-FALSE
if(debugTF){

  ident=ident
  foot=foot.bio
  #ident=boot.outname
  #foot=foot.bio2
  storeTF=FALSE
  res=lon.res
  ncdfpath=new.ncdfpath

}

ctnrt.bio<-function(ident=ident, foot=foot.bio, storeTF=FALSE, res=lon.res, ncdfpath=NULL){

library(ncdf4)
ctpath<-"/uufs/chpc.utah.edu/common/home/lin-group1/group_data/CT-NRT_v2016-1/fluxes/optimized/"

# match 3-houly footprint and 3-houly biosperhic fluxes
# first locate the date and hour for each backwards hours of STILT footprint
recp.yr<-as.numeric(substr(ident,1,4));recp.mon<-as.numeric(substr(ident,6,7))
recp.day<-as.numeric(substr(ident,9,10));recp.hr<-as.numeric(substr(ident,12,13))
foot.hour<-as.numeric(dimnames(foot)[[3]])

# checking if shifted two hours
#foot.hour<-seq(0,69,3)
#foot.hour<-c(0,seq(1,70,3),71)
#foot.hour<-c(0,1,seq(2,71,3))

# attention-- "0" actually means footprint during -3 to -0 hours;"3" actually means footprint during -6 to -3 hours
# thus, add 3 to dimnames(foot)[[3]], and make it negative due to backwards traj
diff.hour<-diff(foot.hour)

# try another way to get back.hour, 08/29/2017, DW
# if we have 0, 1, 4, 7, 10..., which means foot1 is for 0 to -1 hrs, foot2 is for -1 to -4 hrs, ...
# will fix those acutal hours later
max.hr<-max(foot.hour)+diff.hour[1]
back.hour<-c(-1*foot.hour[-1],-1*max.hr)

# check if dimensions match
if(length(back.hour)!=length(foot.hour))cat("ctnrt.bio():dimensions NOT match...\n")

if(F){  # old code
  if(length(unique(diff.hour))==1){ # for even footprint backwards hours
    back.hour<-0-(foot.hour+unique(diff.hour))
  }else{  # for shifted uneven footprint backwards hours
    # e.g., in this case (recp=1000UTC), footprint hours are shifted by an hour
    sel.diff.hour<-diff.hour[1:3] # only select the first 3 differences
    shift.hour<-length(which(sel.diff.hour!=3))
    if(shift.hour==1){  # shift by one hour
      #back.hour<-0-seq(0, 72, 3)-shift.hour # normal case is seq(0, 72, 3)
      #back.hour<-0-seq(0,max.hr,3)-shift.hour # normal case is seq(0, 72, 3)
      #back.hour<-c(back.hour, back.hour[length(back.hour)])
    }else if(shift.hour==2){  # shift by two hours
      #back.hour<-0-seq(0,max.hr, 3)-shift.hour # normal case is seq(0, 72, 3)
      #back.hour<-c(back.hour[1], back.hour)
    }
  }
}

# use weekdayhr() to return the acutal date, weekdayhr<-function(yr,mon,day,hr,runtt,diffGMT=NA)
# runtt needs to be in minutes
back.times<-weekdayhr(recp.yr, recp.mon, recp.day, recp.hr, runtt=back.hour*60)

# fix actual hours, if back.hours%%!=0, 08/29/2017, DW
if(recp.hr%%3!=0){
  stand.hours<-seq(0,21,3)
  back.times[back.times[,"hr"]%%3!=0,"hr"]<-stand.hours[findInterval(back.times[back.times[,"hr"]%%3!=0,"hr"],stand.hours)]
}

foot.hh<-paste(back.times[,"yr"],formatC(back.times[,"mon"],width=2,flag=0),formatC(back.times[,"day"],width=2,flag=0),formatC(back.times[,"hr"],width=2,flag=0),sep="")
foot.dd<-unique(paste(back.times[,"yr"],formatC(back.times[,"mon"],width=2,flag=0),formatC(back.times[,"day"],width=2,flag=0),sep=""))

# if foot hours differ from CT hours, the subroutine cannot find the CT times
# break...
if(unique(back.times[,"hr"]%%3!=0)){
  cat("oco2.ctnrt.bio(): foot hours differ from CT hours, NO CT-NRT files found...\n")
}else{  # if not, we can find CT files

  # grab CT-NRT files for all backwards hours, and then put in a 3D fluxes array
  # create a big 3D array for storing bio fluxes, only covers the regions in footprint [lat, lon, hours]
  foot.lon<-as.numeric(dimnames(foot)[[2]]) # all lower left corner for footprint
  foot.lat<-as.numeric(dimnames(foot)[[1]])

  # use the same dimension as CT-NRT, NOT footprint,
  # remember to flip footprint dimension later
  ctbio.all<-array(0, dim=c(length(foot.lon), length(foot.lat), length(foot.hh)), dimnames=list(foot.lon, foot.lat, foot.hh))
  store.length<-rep(0, length(unique(foot.dd)))

  for (f in 1:length(foot.dd)){

    tmp.foot.hh<-foot.hh[substr(foot.hh,1,8)==foot.dd[f]]

    # open the daily file just once
    ctfile<-list.files(path=ctpath, pattern=paste("CT-NRT.v2016-1.flux1x1.",foot.dd[f],sep=""))
    ctdat<-nc_open(paste(ctpath,ctfile,sep=""))

    # 2D map for each 3 hours
    # spatial resolution, 1 deg in N-S, 1 deg in E-W, move centered lat, lon to lower left
    ct.lat<-ncvar_get(ctdat, "lat")-0.5
    ct.lon<-ncvar_get(ctdat, "lon")-0.5
    ct.time<-ncvar_get(ctdat,"date_components")	# UTC time components
    rownames(ct.time)<-c("year", "month", "day", "hour", "minute", "second")

    # convert to YYYYMMDDHHmm, centered hours, move to the beginning hours
    # 1:30 --> 00:00 to 03:00, use the begining of each 3 hours
    ct.YYYYMMDDHHmmss<-paste(ct.time["year",], formatC(ct.time["month",], width=2, flag=0), formatC(ct.time["day",], width=2, flag=0), formatC(ct.time["hour",]-1, width=2, flag=0), formatC(ct.time["minute",]-30, width=2, flag=0), formatC(ct.time["second",], width=2, flag=0),sep="")
    ct.hh<-substr(ct.YYYYMMDDHHmmss, 1, 10)

    # match footprint hour string with ct hour string, using match for hour index
    # also, select the lat, lon used in footprint, sel.ct.bio[lon, lat, hour]
    match.hh<-match(tmp.foot.hh, ct.hh)
    match.lat<-match(foot.lat, ct.lat)
    match.lon<-match(foot.lon, ct.lon)

    # grab all biosperhic fluxes
    ct.bio<-ncvar_get(ctdat, "bio_flux_opt")	# unit in mol/m2/s, avg fluxes
    dimnames(ct.bio)<-list(ct.lon, ct.lat, ct.hh)
    ct.bio<-ct.bio*1E6	# convert to umol/m2/s

    # select bio for certain lat, lon and time
    sel.ct.bio<-ct.bio[match.lon,match.lat,match.hh]

    # checking, plotting
    plotTF<-FALSE
    if(plotTF){
      cex<-1
      library(fields);library(maps)
      title<-paste("Biospheric CO2 fluxes [umol/m2/s] from 1x1 CT-NRT\nfor Middle East on 12/29/2014 0900-1200UTC\nRiyadh Local Time = UTC+3\nCairo Local Time = UTC+2")
      image.plot(ct.lon[match.lon], ct.lat[match.lat], sel.ct.bio[,,1], zlim=c(-8, 8), col=colorRampPalette(brewer.pal(11,"PiYG"))(100), main=title,xlab="Longitude (degE)",ylab="Latitude (degN)", cex=cex, cex.axis=cex, cex.main=cex, cex.lab=cex)
      text(46.6753,24.7136,"Riyadh", cex=cex)
      text(31.2357,30.0444,"Cairo", cex=cex)
      map("world",add=T,boundary=T)
    }

    # create indices for storing
    store.length[f]<-length(tmp.foot.hh)  # store current files number
    if(f==1){
      min.store<-f; max.store<-store.length[f]
    }else{
      min.store<-sum(store.length[1:(f-1)])+1; max.store<-sum(store.length)
    }
    #print(c(min.store,max.store))

    # put into the big array for storing
    ctbio.all[,,seq(min.store, max.store)]<-sel.ct.bio  # umol/m2/s
    nc_close(ctdat)
  }# end for f

  # Finally, we can match footprint with fluxes
  # flip footprint array first
  flip.foot<-aperm(foot, c(2,1,3))  # now [lon, lat, hr]
  dxco2.bio<-sum(flip.foot*ctbio.all)

  # store CT-NRT fluxes * 3D foot = CO2 contribution map into .nc file
  if(storeTF){

    ident2<-gsub( "&","+",ident)
    cat("ctnrt.bio(): Storing foot x fluxes into ncdf files...");cat("\n")
    netcdf.name<- paste("foot_bio_", ident2, ".nc", sep="")

    # foot.anthro and xfoot.anthro have dims of [LAT, LON]
    dco2.bio<-apply(flip.foot*ctbio.all, c(1,2),sum)
    contri.lat<-as.numeric(rownames(dco2.bio))
    contri.lon<-as.numeric(colnames(dco2.bio))

    x<-ncdim_def("Lon", "degreesE", contri.lon)                 #Set equal to our lat lon vectors we created earlier
    y<-ncdim_def("Lat", "degreesN", contri.lat)

    # flip 2D foot and store footprint in [LAT, LON]
    contri.var<-ncvar_def(name="foot_bio", units="PPM", list(y,x),longname="XCO2 enhancement due to biospheric exchange")
    ncnew<-nc_create(filename=netcdf.name, vars=contri.var)
    ncvar_put(nc=ncnew, varid=contri.var, vals=dco2.bio)            #puts our variable into our netcdf file
    nc_close(ncnew)                                               #Closes our netcdf4 file

    #Move the output file name to our model output directory
    system(paste("mv", netcdf.name, ncdfpath))

  }  # end store nc file

  return(dxco2.bio)
  } # end if, checking whether ct files can be found

} # end of subroutine
