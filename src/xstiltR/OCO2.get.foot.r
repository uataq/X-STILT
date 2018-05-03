# subroutine for obtaining the time-integrated footprint,
# allows you to not regenerate footprint if there is one existed
# two cases for regenerating--
# 1. no footprint ever generated before, where "find.flag==FALSE";
# 2. overwrite==TRUE, overwrite previous footprint anyway, use cautiously

# IF you want to save time by setting overwrite==FALSE, this function will go
#                                and find whether the footprint really exists
# yes -- open existing nc file and return the int.foot
# no -- the function will take care by generating foot

# written by Dien Wu, 09/14/2016

# debug weird footprints over 29.1N, 6E and 14E, DW, 10/16/2017
# add a new flag "sumTF" for returning whether the time-integrated 2D footprint
#              (sumTF=T) or time-varying 3D footprint (sumTF=T), DW, 02/10/2017

get.foot <- function(ident, foot.overwrite=FALSE, trajdat, trajpath, footpath,
	                   foottimes, zlim=c(zbot,ztop), dmassTF, coarse=1,
									   fluxweighting=NULL, numpix.x, numpix.y, storeTF=FALSE,
									   lon.ll, lat.ll,lon.res,lat.res){

library(ncdf4)

# if foottimes has only two components,
# meaning this run only needs the time-integrated foot
if (length(foottimes)==2){
	sumTF<-TRUE
}else{
	sumTF<-FALSE
}

# only store integrated footprint, if sumTF=TRUE
# if storeTF is T, store int foot
intfoot.name <- paste("intfoot", ident, ".nc", sep="")

# only for uneven ident, change "&" to "+"
if (grepl("&", ident) == TRUE){
	ident2 <- gsub("&", "+", ident)
	intfoot.name <- paste("intfoot", ident2,".nc", sep="")
}

## if integrated footprint can be found in the the corresponding directory,
# just use the existing int foot from nc file
find.flag <- TRUE  # initialize by TRUE

if (foot.overwrite == FALSE){
	find.file <- list.files(path=footpath, pattern=intfoot.name)
	if(length(find.file)==0)find.flag <- FALSE
}

# #if we don't want to regenerate the foot and there is an existing nc file
# then we open the nc file for intfoot
if (foot.overwrite==FALSE & find.flag==TRUE){

	cat("get.foot(): FOUND the footprint, open the file and grab footprint...\n")

	# since we find the existing file, find.file is the filename
	xfootdat <- nc_open(paste(footpath, find.file, sep=""))
	xfoot    <- ncvar_get(xfootdat, "footprint")
	dimnames(xfoot) <- list(ncvar_get(xfootdat,"Lat"), ncvar_get(xfootdat,"Lon"))
	nc_close(xfootdat)
}  # end not overwritting AND file existing


## condition for generating foot, no existing foot found OR foot.overwrite=TRUE
if (foot.overwrite==TRUE | find.flag==FALSE){

	cat("get.foot(): Generating footprint, call Trajecfoot.r...\n")

	# pathname in trajecfoot.r is the path for storing .RData
	# !!! Trajecfoot() needs the .RData file stored in dir,
	# thus never turn storeTF to FALSE when calling weight.trajecfoot()
	# can just pass on the trajdat to "part" for trajecfoot()
	foot <- Trajecfoot(ident=ident, part=trajdat, pathname=trajpath, zlim=zlim,
		                 foottimes=foottimes, dmassTF=dmassTF, fluxweighting=NULL,
										 coarse=1, lon.ll=lon.ll, lat.ll=lat.ll, numpix.x=numpix.x,
										 numpix.y=numpix.y, lon.res=lon.res, lat.res=lat.res)

  # if need to store time-integrated footprint
	if(sumTF){

		if (length(foottimes) == 2){
			xfoot <- foot[,,1]
		}else{
			xfoot <- apply(foot, c(1, 2), sum)
		}

		if(storeTF){

			cat("get.foot(): Saving time-integrated footprint as .nc file...\n\n")

			# foot and xfoot have dims of [LAT, LON]
			foot.lat <- as.numeric(rownames(xfoot))
			foot.lon <- as.numeric(colnames(xfoot))

			#Set equal to our lat lon vectors we created earlier
			x <- ncdim_def("Lon", "degreesE", foot.lon)
			y <- ncdim_def("Lat", "degreesN", foot.lat)

			# flip 2D foot and store footprint in [LAT, LON]
			foot.var <- ncvar_def(name="footprint", units="PPM/(umoles/m2/s)",
			                      list(y,x), longname="time-integrated 2D footprint")

			ncnew <- nc_create(filename=intfoot.name, vars=foot.var)

			# puts our variable into our netcdf file
			ncvar_put(nc=ncnew, varid=foot.var, vals=xfoot)
			nc_close(ncnew)  # Closes our netcdf4 file

			#Move the output file name to our model output directory
			system(paste("mv", intfoot.name, footpath))
		} # end storeTF

	}	# end sumTF
}	# end (re)generating and storing int foot

if(sumTF){
	# return time-int footprint due to monthly temporal resolution of ODIAC
	return(xfoot)
}else{
	# return time-varying footprint due to hourly ODIAC
	return(foot)
}  # end if sumTF

}



# end of subroutine
