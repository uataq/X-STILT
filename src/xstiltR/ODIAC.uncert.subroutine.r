# subroutine to readin ODIAC emissions uncertainty, and then couple with STILT footprint
# originated from ODIAC.subroutine.r, DW, 09/04/2017

debugTF<-FALSE
# input from main script
if(debugTF){

  foot=foot.anthro
  res=lon.res
  ncdfpath=new.ncdfpath
  storeTF=FALSE
  ident=ident
  site=site
  odiac.uncert=anthro.uncert
}

odiac.anthro.uncert<-function(foot=foot, odiac.uncert=anthro.uncert, storeTF=FALSE, res=1/10, ident=ident, ncdfpath=ncdfpath){

library(ncdf4)

# grab year and month from ident
yr<-substr(ident,1,4);mon<-substr(ident,6,7)

# dealing with footprint
foot.lat<-as.numeric(rownames(foot))
foot.lon<-as.numeric(colnames(foot))

# determine whether to use hourly emissions, based on footprint that passed on
hourlyTF<-FALSE
if(length(dim(foot))==3)hourlyTF<-TRUE

if (res == 1/10){

  # get lat lon of uncertainty
  uncert.lon<-as.numeric(dimnames(odiac.uncert)[[1]])
  uncert.lat<-as.numeric(dimnames(odiac.uncert)[[2]])
  if(hourlyTF)uncert.hr<-as.numeric(dimnames(odiac.uncert)[[3]])

  # subset the max overlapping region between foot and emission
  # foot lat lon fall into the range of emission lat lon
  minlat<-max(min(uncert.lat), min(foot.lat));maxlat<-min(max(uncert.lat), max(foot.lat))
  minlon<-max(min(uncert.lon), min(foot.lon));maxlon<-min(max(uncert.lon), max(foot.lon))

  lat.index<-uncert.lat<= maxlat & uncert.lat>= minlat
  lon.index<-uncert.lon<= maxlon & uncert.lon>= minlon

  if(hourlyTF){   # 3D hourly file
    # flip uncertainty
    f.uncert<-aperm(odiac.uncert,c(2,1,3))
    sel.foot<-foot[foot.lat<= maxlat & foot.lat>= minlat, foot.lon<= maxlon & foot.lon>= minlon,]
    sel.uncert<-f.uncert[lat.index,lon.index,]
  }else{
    # flip uncertainty
    f.uncert<-t(odiac.uncert)
    sel.foot<-foot[foot.lat<= maxlat & foot.lat>= minlat, foot.lon<= maxlon & foot.lon>= minlon]
    sel.uncert<-f.uncert[lat.index,lon.index]
  }

  # NOW, foot and absolute uncertainty have the same dim, multipling...
  dxco2.anthro.uncert<-sel.foot*sel.uncert
  sum.dxco2.uncert<-sum(dxco2.anthro.uncert)

} # end using 1/10 deg ODIAC uncertainty

# store ODIAC emission * integrated footprint = CO2 contribution map into .nc file
if(storeTF){

  ident2<-gsub( "&","+",ident)
  cat("ODIAC.anthro(): Storing footxuncert into ncdf files...");cat("\n")
  netcdf.name<- paste("foot_anthro_uncert_", ident2, ".nc", sep="")

  # foot.anthro and xfoot.anthro have dims of [LAT, LON]
  contri.lat<-as.numeric(rownames(dxco2.anthro.uncert))
  contri.lon<-as.numeric(colnames(dxco2.anthro.uncert))

  x<-ncdim_def("Lon", "degreesE", contri.lon)                 #Set equal to our lat lon vectors we created earlier
  y<-ncdim_def("Lat", "degreesN", contri.lat)

  # flip 2D foot and store footprint in [LAT, LON]
  var<-ncvar_def(name="foot_anthro_uncert", units="PPM", list(y,x),longname="dXCO2 due to absolute ODIAC emission uncertainty")
  ncnew<-nc_create(filename=netcdf.name, vars=var)
  ncvar_put(nc=ncnew, varid=var, vals=dxco2.anthro.uncert)            #puts our variable into our netcdf file
  nc_close(ncnew)                                               #Closes our netcdf4 file

  #Move the output file name to our model output directory
  system(paste("mv", netcdf.name, ncdfpath))
}  # end store nc file

return(sum.dxco2.uncert)
} # end of subroutine
