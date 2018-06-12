# subroutine to get ground height using trajwind()
# DW, 10/20/2017

# Inputs variables:
# 1) "recp.info" from ident.to.info(),
# including receptor time, lat, lon (recp.info) and release levels (agl.info)
# 2) "nummodel" for copy number
# 3) "outpath, metpath, metfile" for trajwind()...

# Updates --
# add 0.5 deg GDAS, DW, 01/25/2018
# interpolate ground heights from multiple receptors,
# add vector forms of lat/lon/agl, DW, 05/02/2018

get.grdhgt <- function(recp.info, nummodel, agl=10, nhrs=-1, rundir, outpath,
                       site, metpath, metfile){

	# METHOD 1--call fortran code profile() for directly extracting ground height
	# profile.info<-profileARL(LAT=recp.info$recp.lat, LON=recp.info$recp.lon,
  #                          metdir=metpath,metfile=metfile,ttoff=0)

	# METHOD 2--use trajwind() for interpolating ground height
  # if nhrs is negative, u, v winds are backward wind speed and directions
	# but, we are interpolating the ground height, so no need to flip it
  # feed lat,lon,agl as vectors, same dimension; still one outname !!!
  varsout <- c("time","index","lon","lat","agl","grdht","zi","temp","pres")

  # vector of AGLs
  agl <- rep(10, nrow(recp.info))  # if multiple receptor lat/lon

  outname <- paste(unique(recp.info$timestr), "_", site, sep="")

  # if outpath did not end with "/"
  if(substr(metpath, nchar(metpath), nchar(metpath))!="/"){
    metpath <- paste(metpath, "/",sep="")
  }

  # if outpath did not end with "/"
  if(substr(rundir, nchar(rundir), nchar(rundir))!="/"){
    rundir <- paste(rundir, "/",sep="")
  }

  # yr, mon, day, hr should be the same,
  # trajwind() can return matrix if lat/lon/agl is in vector form
	trajwind.info <- trajwind(yr =unique(recp.info$recp.year - 2000),
                            mon=unique(recp.info$recp.mon),
                            day=unique(recp.info$recp.day),
                            hr =unique(recp.info$recp.hour), outname = outname,
                            lat=recp.info$recp.lat, lon=recp.info$recp.lon,
                            agl=agl, nhrs=nhrs, metlib=metpath, metfile=metfile,
                            rundir=rundir, outpath=outpath, varsout=varsout,
                            nummodel=nummodel, metd=c("fnl","awrf"))

  # add trajwind() output onto original "recp.info"
  recp.info$recp.grdhgt <- trajwind.info[, "grdht"]  # ground height in meters

  # if nhrs<0, flip sign of winds, m/s
  recp.info$recp.ubar <- trajwind.info[, "ubar"] * sign(nhrs)
  recp.info$recp.vbar <- trajwind.info[, "vbar"] * sign(nhrs)

  recp.info$recp.zi   <- trajwind.info[, "zi"]    # PBL heights in meters
  recp.info$recp.temp <- trajwind.info[, "temp"]  # temp in K
  recp.info$recp.pres <- trajwind.info[, "pres"]  # pressure in mb

	return(recp.info)
}
