# subroutine to readin ODIAC emissions, v2015a and then couple with STILT footprint
# FOR NOW, use integrated NEWLY WEIGHTED footprint due to lack of diurnal cycle of anthro emissions
# since footprints have already been weighted by AK and PW, just multiple emission with 2D footprint map
# written by DIEN WU, 09/13/2016

# ADD ODIACv2016, flag 1 for using v2015a, flag 2 for v2016, DW 02/14/2017
# note that v2015a does not have emission for year 2015
# ADD PRD ODIAC emissions, DW 03/08/2017
# ADD TIMES hourly scaling factors for ODIACv2016, DW 03/08/2017

debugTF<-FALSE
# input from main script
if(debugTF){

foot=foot.anthro
res=lon.res
ncdfpath=new.ncdfpath
storeTF=FALSE
ident=ident
odiac.ver=2
site=site
odiac.co2=odiac.co2
}

odiac.anthro<-function(foot=foot, odiac.co2=NULL, odiac.domain=odiac.domain, odiac.ver=2, storeTF=FALSE, res=res, ident=ident, ncdfpath=ncdfpath, site=site){

library(ncdf4)

# grab year and month from ident
yr<-substr(ident,1,4);mon<-substr(ident,6,7)

# dealing with footprint
foot.lat<-as.numeric(rownames(foot))
foot.lon<-as.numeric(colnames(foot))

# determine whether to use hourly emissions, based on footprint that passed on
hourlyTF<-FALSE
if(length(dim(foot))==3)hourlyTF<-TRUE

# define ODIAC versions and domains
if(odiac.ver==1)odiac.vname<-"2015a"
if(odiac.ver==2)odiac.vname<-"2016"

path<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_input/anthro_inventories/ODIAC/v",odiac.vname,"/",sep="")	# stored ODIAC
if(hourlyTF)path<-paste(path,"hourly/",site,"/",sep="")

if (res == 1.0){

  file<-paste("odiac",odiac.vname,"_1.0x1.0d_",yr,"_rev1.nc",sep="")
  odiac.dat<-nc_open(paste(path,file,sep=""))

  # ODIAC contains two portions, ff emissions over land and emissions from international aviation and marine bunker
  # CO2 have unit of g-C/m2/day, [LON, LAT, month], centered lat, lon, 1deg*1deg
  odiac.co2<-ncvar_get(odiac.dat,"land")
  odiac.lat<-ncvar_get(odiac.dat,"lat")-0.5
  odiac.lon<-ncvar_get(odiac.dat,"lon")-0.5
  odiac.mon<-ncvar_get(odiac.dat,"month")
  dimnames(odiac.co2)<-list(odiac.lon,odiac.lat,odiac.mon)

  # determine the month and grab the emission grid
  stilt.mon<-as.numeric(substr(ident, 6, 7))
  sel.odiac.co2<-odiac.co2[,,stilt.mon]

  # convert unit g-C/m2/d to umol-CO2/m2/s,
  # transpose the matrix to make it [LAT, LON]
  sel.odiac.co2<-sel.odiac.co2/12*1E6/24/60/60
  t.sel.odiac.co2<-t(sel.odiac.co2)

  # since this version is globally, grab the emission where footprint fall into
  sub.odiac.co2<-t.sel.odiac.co2[odiac.lat<=max(foot.lat) & odiac.lat>=min(foot.lat), odiac.lon<=max(foot.lon) & odiac.lon>=min(foot.lon)]

  # couple with footprint
  dco2.anthro<-foot*sub.odiac.co2
  dxco2.anthro<-sum(dco2.anthro)

  nc_close(odiac.dat)

}  # end using 1deg ODIA

if (res == 1/120){

  # readin ODIAC if emission is NULL
  if(length(odiac.co2)==0){
    if(hourlyTF){file<-paste("odiac",odiac.vname,"_1kmx1km_",yr,mon,"_",odiac.domain,"_hrly.nc",sep="")}else{file<-paste("odiac",odiac.vname,"_1kmx1km_",yr,mon,"_",odiac.domain,".nc",sep="")}
    odiac.dat<-nc_open(paste(path,file,sep=""))

    # grab emission grid, [LAT, LON, HR.back], lat, lon have already converted to lower left for this version
    # hours back are negative,
    cat("ODIAC.anthro(): take a while to readin hourly ODIAC CO2 emissions...\n")
    odiac.co2<-ncvar_get(odiac.dat, "odiac_co2_emiss")
    odiac.lat<-ncvar_get(odiac.dat,"lat");odiac.lon<-ncvar_get(odiac.dat,"lon")
    if(hourlyTF){odiac.hr<-ncvar_get(odiac.dat,"hr");dimnames(odiac.co2)<-list(odiac.lat, odiac.lon, odiac.hr)}else{dimnames(odiac.co2)<-list(odiac.lat, odiac.lon)} # endif hourlyTF
    nc_close(odiac.dat)

  }else{
    odiac.lat<-as.numeric(dimnames(odiac.co2)[[1]]);odiac.lon<-as.numeric(dimnames(odiac.co2)[[2]])
    if(hourlyTF)odiac.hr<-as.numeric(dimnames(odiac.co2)[[3]])
  } # endif length()

  # subset the max overlapping region between foot and emission
  # foot lat lon fall into the range of emission lat lon
  minlat<-max(min(odiac.lat), min(foot.lat));maxlat<-min(max(odiac.lat), max(foot.lat))
  minlon<-max(min(odiac.lon), min(foot.lon));maxlon<-min(max(odiac.lon), max(foot.lon))

  if(hourlyTF){   # 3D hourly file
    sel.foot<-foot[foot.lat<= maxlat & foot.lat>= minlat, foot.lon<= maxlon & foot.lon>= minlon,]
    sel.odiac.co2<-odiac.co2[odiac.lat<= maxlat & odiac.lat>= minlat, odiac.lon<= maxlon & odiac.lon>= minlon,]
  }else{
    sel.foot<-foot[foot.lat<= maxlat & foot.lat>= minlat, foot.lon<= maxlon & foot.lon>= minlon]
    sel.odiac.co2<-odiac.co2[odiac.lat<= maxlat & odiac.lat>= minlat, odiac.lon<= maxlon & odiac.lon>= minlon]
  }

  # NOW, foot and emiss have the same dim, multiple to get a dCO2
  dco2.anthro<-sel.foot*sel.odiac.co2
  dxco2.anthro<-sum(dco2.anthro)

} # end using 1/120 deg ODIAC

# store ODIAC emission * integrated footprint = CO2 contribution map into .nc file
if(storeTF){

  ident2<-gsub( "&","+",ident)
  cat("ODIAC.anthro(): Storing footxemission into ncdf files...");cat("\n")
  netcdf.name<- paste("foot_anthro_", ident2, ".nc", sep="")

  # foot.anthro and xfoot.anthro have dims of [LAT, LON]
  contri.lat<-as.numeric(rownames(dco2.anthro))
  contri.lon<-as.numeric(colnames(dco2.anthro))

  x<-ncdim_def("Lon", "degreesE", contri.lon)                 #Set equal to our lat lon vectors we created earlier
  y<-ncdim_def("Lat", "degreesN", contri.lat)

  # flip 2D foot and store footprint in [LAT, LON]
  foot.var<-ncvar_def(name="foot_anthro", units="PPM", list(y,x),longname="dCO2 due to ODIAC emission")
  ncnew<-nc_create(filename=netcdf.name, vars=foot.var)
  ncvar_put(nc=ncnew, varid=foot.var, vals=dco2.anthro)            #puts our variable into our netcdf file
  nc_close(ncnew)                                               #Closes our netcdf4 file

  #Move the output file name to our model output directory
  system(paste("mv", netcdf.name, ncdfpath))
}  # end store nc file

return(dxco2.anthro)
} # end of subroutine
