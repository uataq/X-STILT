# subroutine for obtaining the time-integrated footprint,
# allows you to not regenerate footprint if there is one existed
# two cases for regenerating--
# 1. no footprint ever generated before, where "find.flag==FALSE";
# 2. overwrite==TRUE, overwrite previous footprint anyway, use cautiously

# IF you want to save time by setting overwrite==FALSE, this function will go and find whether the footprint really exists
# yes -- open existing nc file and return the int.foot
# no -- the function will take care by generating foot
# Written by DIEN WU, 09/14/2016

# add a new flag "sumTF" for returning whether the time-integrated 2D footprint (sumTF=T) or time-varying 3D footprint (sumTF=T), 02/10/2017

debugTF<-FALSE
if (debugTF){
	ident=ident
	foot.overwrite=TRUE
	trajpath=new.trajpath
	footpath=new.intpath
	foottimes=c(0,72)
	zlim=c(0,0)
	fluxweighting=NULL
	coarse=1
	numpix.x=3000
	numpix.y=2400
	lon.ll=10
	lat.ll=20
	lon.res=1/120
	lat.res=1/120
	storeTF=FALSE
	part=new.trajdat
	sumTF=FALSE

	# debug DW, 10/16/2017
	plot.foot<-function(foot=foot){
		library(reshape)
		melt.foot<-melt(foot[,,1])
		colnames(melt.foot)<-c("lat","lon","foot")
		signal<-1E-8
		melt.foot<-melt.foot[melt.foot$foot>signal,]
		if(site=="Riyadh"){minlat<-17;maxlat<-25;minlon<-0;maxlon<-50}
		if(site=="Cairo"){minlat<-20;maxlat<-35;minlon<-5;maxlon<-33}
		m1<-ggplot.map(minlat=minlat,maxlat=maxlat,minlon=minlon,maxlon=maxlon,ocean.col="lightsteelblue2",land.col="white",land.outline="gray30")
		shift.lon<-0;shift.lat<-0;font.size=rel(1.2)

		if(met=="1km")met<-"wrf"
		title<-paste(toupper(met),"-STILT footprint (>",toupper(signal),"umol/m2/s)\nlon.ll=",lon.ll,"numpix.x=",numpix.x)
		p1<-m1+labs(title=title,x="Longitude [E]",y="Latitude [N]")+theme_bw()#+coord_equal(ratio=1.1)
		p1<-p1+geom_raster(data=melt.foot,aes(x=lon,y=lat,fill=foot),alpha=0.8,interpolate=F)#+facet_grid(.~fac)
		p1<-p1+scale_fill_gradient(limits=c(signal,4E-4),name="STILT Footprint",low="gray90",high="black",trans="log10",breaks<-c(signal,1E-4))+theme(legend.key.width=unit(0.5, "cm"),legend.key.height=unit(2, "cm"))
		recp.loc<-data.frame(lon=31.0144,lat=29.1483)
		p2<-p1+geom_point(data=recp.loc,aes(x=lon,y=lat),colour="purple",size=3)
		return(p2)
	}
}

get.foot<-function(ident=ident, foot.overwrite=FALSE, part=new.trajdat, trajpath=trajpath, footpath=footpath, foottimes=foottimes,zlim=c(zbot,ztop),dmassTF=dmassTF,fluxweighting=NULL,coarse=1,vegpath=vegpath,numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res, storeTF=FALSE){

# if foottimes has only two components, meaning this run only needs the time-integrated foot
if(length(foottimes)==2){sumTF<-TRUE}else{sumTF<-FALSE}

# only store integrated footprint, if sumTF=TRUE
intfoot.name<- paste("intfoot", ident,".nc", sep="")	# if storeTF is T, store int foot

# only for uneven ident, change "&"
if(grepl("&",ident)==TRUE){ident2<-gsub("&","+",ident);intfoot.name<- paste("intfoot", ident2,".nc", sep="")}

# if integrated footprint have already been found in the the corresponding directory, just use the existing int foot from nc file
find.flag<-TRUE	# initialize by TRUE
if (foot.overwrite==FALSE){
	find.file<-list.files(path=footpath, pattern=intfoot.name)
	if(length(find.file)==0)find.flag<-FALSE
}

# if we dont want to regenerate the foot and there is an existing nc file
# then we open the nc file for intfoot
if (foot.overwrite==FALSE & find.flag==TRUE){

	cat("get.foot(): FOUND the footprint, open the file and grab footprint...");cat("\n")
	# since we find the existing file, find.file is the filename
	library(ncdf4)
	xfootdat<-nc_open(paste(footpath, find.file, sep=""))
	xfoot<-ncvar_get(xfootdat, "footprint")
	dimnames(xfoot)<-list(ncvar_get(xfootdat,"Lat"), ncvar_get(xfootdat,"Lon"))
	nc_close(xfootdat)
}  # end not overwritting AND file existing


# condition for generating foot, no existing foot found OR foot.overwrite=TRUE
if (foot.overwrite==TRUE | find.flag==FALSE){

	cat("get.foot(): Regenerate footprint(overwrite or no foot found), calling Trajecfoot.r...");cat("\n")

	# pathname in trajecfoot.r is the path for storing .RData
	# !!! Trajecfoot() needs the .RData file stored in dir, thus never turn storeTF to FALSE when calling weight.trajecfoot()
	# can just pass on the trajdat to "part" for trajecfoot()
	foot<-Trajecfoot(ident=ident,part=part, pathname=trajpath,foottimes=foottimes,zlim=zlim,dmassTF=dmassTF,fluxweighting=NULL,coarse=1,vegpath=vegpath,numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res)
	#foot<-Trajecfoot(ident=ident, pathname=trajpath,foottimes=foottimes,zlim=c(zbot,ztop),fluxweighting=NULL,coarse=1,vegpath=vegpath,numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res)

	#p<-plot.foot(foot)

	# debug weird footprints over 29.1N, 6E and 14E, DW, 10/16/2017
  #lat<-rownames(foot)
  #lon<-colnames(foot)
  #selfoot1<-foot[lat<29.2&lat>29,lon>6&lon<6.1,]
  #selfoot2<-foot[lat<29.4&lat>29.3,lon>14.2&lon<14.4,]

	#As an additional step lets save the time-integrated footprint as a netcdf file as well!
	# if we let Trajectfoot integrates over every certain hours, do apply for obtaining 2D foot
	# if we've already let Trajectfoot integrates over all nhrs hours (foottime=0, 72), just one for the 3rd dimention, just grab the first 2 dims now
	    #  e.g., > str(foot)
		# num [1:6000, 1:7200, 1] 0 0 0 0 0 0 0 0 0 0 ...
		# - attr(*, "dimnames")=List of 3
		# ..$ : chr [1:6000] "0" "0.00833333333333333" "0.0166666666666667" "0.025" ...
		# ..$ : chr [1:7200] "0" "0.00833333333333333" "0.0166666666666667" "0.025" ...
		# ..$ : chr "0"

	if(sumTF){	# if need time-integrated footprint

		if(length(foottimes)==2){
			xfoot<-foot[,,1]
		}else{
			xfoot<-apply(foot, c(1,2), sum)
		}

		if(storeTF){
			library(ncdf4)
			cat("get.foot(): Saving footprint grid to netcdf output file..."); cat("\n"); cat("\n")
			# foot and xfoot have dims of [LAT, LON]
			foot.lat<-as.numeric(rownames(xfoot))
			foot.lon<-as.numeric(colnames(xfoot))
			x<-ncdim_def("Lon", "degreesE", foot.lon)                 #Set equal to our lat lon vectors we created earlier
			y<-ncdim_def("Lat", "degreesN", foot.lat)

			# flip 2D foot and store footprint in [LAT, LON]
			foot.var<-ncvar_def(name="footprint", units="PPM/(umoles/m2/s)", list(y,x),longname="time-integrated 2D footprint")
			ncnew<-nc_create(filename=intfoot.name, vars=foot.var)
			ncvar_put(nc=ncnew, varid=foot.var, vals=xfoot)            #puts our variable into our netcdf file
			nc_close(ncnew)                                               #Closes our netcdf4 file

			#Move the output file name to our model output directory
			system(paste("mv", intfoot.name, footpath))
		} # end storeTF

	}	# end sumTF
}	# end (re)generating and storing int foot

if(sumTF){
	return(xfoot)	# return time-int footprint due to monthly temporal resolution of ODIAC
}else{
	return(foot)	# return time-varying footprint due to hourly ODIAC
}

} # end of subroutine
