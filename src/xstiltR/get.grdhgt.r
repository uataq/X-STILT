# subroutine to get ground height using trajwind()
# DW, 10/20/2017

# Inputs variables:
# recp.info: output from ident.to.info(), including receptor time, lat, lon (recp.info) and release levels (agl.info)

# Amendments --
# add 0.5 deg GDAS, DW, 01/25/2018
# add another approach to interpolate ground height, using Fortran function "profile", call profileARL() written by John, DW, 01/29/2018
# interpolate ground heights of multiple receptors, vector forms of lat/lon/agl, DW, 04/18/2018

get.grdhgt<-function(recp.info,nummodel,agl=10,nhrs=-1,rundir="/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/STILT_Exe/",
                     grdhgt.path="/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/trajwind/Riyadh/",
										 metpath="/uufs/chpc.utah.edu/common/home/u0947337/",metfile){

	# METHOD 1--call fortran code profile() for directly extracting ground height
	#profile.info<-profileARL(LAT=recp.info$recp.lat,LON=recp.info$recp.lon,metdir=metpath,metfile=metfile,ttoff=0)

	# METHOD 2--use trajwind() for interpolating ground height
	# backward for 1hour, if nhrs is negative, u, v winds are backward wind speed and directions
	# but, we are interpolating the ground height, so no need to flip it
  # feed lat,lon,agl as vectors, same dimension; still one outname !!!
	trajwind.info<-trajwind(yr=unique(recp.info$recp.year-2000),mon=unique(recp.info$recp.mon),day=unique(recp.info$recp.day),hr=unique(recp.info$recp.hour),
                          lat=recp.info$recp.lat,lon=recp.info$recp.lon,agl=rep(10,nrow(recp.info)),nhrs=nhrs,
													metlib=metpath,metfile=metfile,rundir=rundir,outpath=grdhgt.path,nummodel=nummodel,
													metd=c("fnl","awrf"),varsout=c("time","index","lon","lat","agl","grdht","zi","temp","pres"))

  recp.info$recp.grdhgt<-trajwind.info[,"grdht"]  # ground height in meters
  recp.info$recp.ubar<-trajwind.info[,"ubar"]*nhrs/abs(nhrs)  # if nhrs<0, flip sign of winds
  recp.info$recp.vbar<-trajwind.info[,"vbar"]*nhrs/abs(nhrs)

  recp.info$recp.zi<-trajwind.info[,"zi"]  # PBL heights in meters
  recp.info$recp.temp<-trajwind.info[,"temp"]  # temp in K
  recp.info$recp.pres<-trajwind.info[,"pres"]  # pressure in mb
	return(recp.info)
}
