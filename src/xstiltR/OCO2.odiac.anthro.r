# subroutine to readin ODIAC emissions, v2015a and then couple with STILT footprint
# FOR NOW, use integrated NEWLY WEIGHTED footprint due to lack of diurnal cycle of anthro emissions
# since footprints have already been weighted by AK and PW, just multiple emission with 2D footprint map
# written by Dien Wu, 09/13/2016

# Updates:
# ADD ODIACv2016, flag 1 for using v2015a, flag 2 for v2016, DW 02/14/2017
# note that v2015a does not have emission for year 2015
# ADD PRD ODIAC emissions, DW 03/08/2017
# ADD TIMES hourly scaling factors for ODIACv2016, DW 03/08/2017
# Get rid of variable "odiac.vname", always preprocess and read ODIAC emission before call this function...

odiac.anthro<-function(foot,odiac.co2,ident,storeTF=F,ncdfpath){

  library(ncdf4)

  # grab year and month from ident
  yr<-substr(ident,1,4);mon<-substr(ident,6,7)

  # grab footprint dimensions
  foot.lat<-as.numeric(rownames(foot))
  foot.lon<-as.numeric(colnames(foot))

  # determine whether to use hourly emissions, based on footprint that passed on
  hourlyTF<-FALSE;if(length(dim(foot))==3)hourlyTF<-TRUE

  # grab emission dimensions
  odiac.lat<-as.numeric(dimnames(odiac.co2)[[1]]);odiac.lon<-as.numeric(dimnames(odiac.co2)[[2]])
  if(hourlyTF)odiac.hr<-as.numeric(dimnames(odiac.co2)[[3]])

  ### find the overlap region between foot and emission
  # foot lat lon fall into the range of emission lat lon
  minlat<-max(min(odiac.lat), min(foot.lat));maxlat<-min(max(odiac.lat), max(foot.lat))
  minlon<-max(min(odiac.lon), min(foot.lon));maxlon<-min(max(odiac.lon), max(foot.lon))

  # select footprint and emission fields
  if(hourlyTF){   # 3D hourly file
    sel.foot<-foot[foot.lat<= maxlat & foot.lat>= minlat, foot.lon<= maxlon & foot.lon>= minlon,]
    sel.odiac.co2<-odiac.co2[odiac.lat<= maxlat & odiac.lat>= minlat, odiac.lon<= maxlon & odiac.lon>= minlon,]
  }else{
    sel.foot<-foot[foot.lat<= maxlat & foot.lat>= minlat, foot.lon<= maxlon & foot.lon>= minlon]
    sel.odiac.co2<-odiac.co2[odiac.lat<= maxlat & odiac.lat>= minlat, odiac.lon<= maxlon & odiac.lon>= minlon]
  }

  # NOW, foot and emiss should have the same dimension, multiple them to get contribution map of CO2 enhancements
  # sum the map to get the XCO2 enhancements, note that AK and PW have been incorporated in footprint
  dco2.anthro<-sel.foot*sel.odiac.co2
  dxco2.anthro<-sum(dco2.anthro)

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
    ncvar_put(nc=ncnew, varid=foot.var, vals=dco2.anthro)   #puts our variable into our netcdf file
    nc_close(ncnew)   # close our netcdf4 file

    #Move the output file name to our model output directory
    system(paste("mv", netcdf.name, ncdfpath))
  }  # end store nc file

  return(dxco2.anthro)

} # end of subroutine
