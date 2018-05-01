# create a namelist for running X-STILT by DW, 04/18/2018
# originated from "STILT_OCO2_sim.r"

# required functions: ident.to.info()

##### source all functions and load all libraries
# change working directory ***
workpath <- "/uufs/chpc.utah.edu/common/home/lin-group1/wde/github/XSTILT"
setwd(workpath)
source(file.path(workpath,"src/sourceall.r"))  # source all functions
library(ncdf4);library(ggplot2)

#------------------------------ STEP 1 --------------------------------------- #
#### CHOOSE CITIES AND USED MET RESOLUTION ***
index <- 1
site  <- c("Riyadh","Cairo","PRD")[index]
met   <- c("1km","gdas1","gdas0p5")[3] # customized WRF, 1 or 0.5deg GDAS

# spatial domains for generating footprints or grabbing ODIAC emissions,
# these variables will determine the numpix.x, numpix.y, lon.ll, lat.ll;
# e.g., for Riyadh, 	# 0-50N, 0-60E,  ***
if(site == "Riyadh"){min.lon<-0;  max.lon<-60; min.lat<-0;  max.lat<-50}
if(site == "PRD")   {min.lon<-10; max.lon<-50; min.lat<-20; max.lat<-130}
if(site == "Cairo") {min.lon<-0;  max.lon<-60; min.lat<-0;  max.lat<-50}

#### CHOOSE OVERPASSED TIMESTR ***
tt <- 2  # which track to model

# input all track numbers to be modeled, can be YYYYMMDD OR YYYYMMDDHH ***
riyadh.timestr <- c("2014122710", "2014122910", "2015012810", "2015081710",
                    "2015111210", "2015121610", "2016011510", "2016021610",
                    "2016072510")
cairo.timestr  <- c("20150228", "20150318", "20160224")
prd.timestr    <- "2015011505"

# final timestr, YYYYMMDD and file strings for trajec
track.timestr <- c(riyadh.timestr[tt],cairo.timestr[tt],prd.timestr)[index]
filestr <- paste(substr(track.timestr,1,4), "x", substr(track.timestr,5,6), "x",
                 substr(track.timestr,7,8), sep="")

#### CHOOSE THE CO2 SOURCES/SINKS TO BE MODELED ***
# total 5 components, no wildfire component for now, but can be easily added
ss <- c(1,2,3,4,5)
sel <- c("anthro", "bio", "ocean", "edp", "apriori")[ss]

#------------------------------ STEP 2 --------------------------------------- #
#### TURN ON/OFF FLAGS for XSTILT setups ***
ppTF <- F       # whether use ODIAC with addtional PP10+PP14
selTF <-F       # true for only use 1-day trajs.
dmassTF <-T     # dmassTF when generating trajecfoot().
storeTF <- F    # whether store model results;
                # e.g., footprint, co2 contribution maps, weighted trajec.

trajecTF <- T   # whether running trajec
nummodel <- 10	# copy for running trajwind() or trajec()
grdhgtTF <- !trajecTF  # whether run ground heights only when trajecTF==FALSE

hourlyTF <- F   # true for using hourly ODIAC emissions;
                # false for using monthly mean emissions.
columnTF <- T   # whether a column receptor or fixed receptor
forwardTF <- F  # forward or backward traj

uncertTF <- F   # sim using trajec with transport error components
biascorrTF <- F # whether use bias-corrected trajec

ak.weight <- T  # whether weighting footprint by averaging kernels
pw.weight <- T  # whether weighting footprint by pressure weighting functions
if(columnTF)pw.weight <-T  # if release from a column, always XCO2

#### INPUT/CHANGE PATHS ***
# outpath for storing model intermediate and results
# inpath for OCO2, ODIAC, CT..etc.
outpath <- "/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2"
inpath <- "/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_input"

### OCO-2 path and version ***
oco2.version <- c("b7rb","b8r")[1]
ocopath  <- file.path(inpath , paste("OCO-2/OCO2_lite_",oco2.version,sep=""))

## STILT paths
trajpath <- file.path(outpath, "Traj")  # trajec path
ncdfpath <- file.path(outpath, "NetCDF")  # footprint and contribution paths
windpath <- file.path(outpath, "trajwind",site) # path for storing trajwind()

# Where to run STILT, where Copy lies ***
rundir <- "/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/STILT_Exe/"

# path for the ARL format of WRF and GDAS, CANNOT BE TOO LONG ***
metpath <- "/uufs/chpc.utah.edu/common/home/u0947337/"
#if(met=="gdas0p5")metpath <- file.path(metpath, "GDAS0p5")

# compute string for original trajec path
str<-NULL
if(uncertTF==TRUE)str <-"_uncert"
if(forwardTF==TRUE)str <-"_forward"
if(biascorrTF==TRUE)str <-"_bc"
if(columnTF==TRUE)aglstr <-"multi_agl"
if(columnTF==FALSE)aglstr <-"fixed_agl"

## GET ORIGINAL TRAJEC (NO AK WEIGHTING YET) PATH ***
orig.trajpath <- file.path(trajpath, site, met, aglstr, "orig_traj")
if(length(str) != 0)orig.trajpath <- file.path(orig.trajpath, str)
if(length(str) == 0)orig.trajpath <- file.path(orig.trajpath, track.timestr)

### paths for storing ak.pw weighted traj (new.trajpath), ***
# weighted column footprint (new.intpath),
# XCO2 contribution maps for each source/sink (new.ncdfpath)
new.trajpath <- file.path(trajpath, site, met, aglstr, "akpw_traj")
new.intpath  <- file.path(ncdfpath, site, met, aglstr, "akpw_intfoot")
new.ncdfpath <- file.path(ncdfpath, site, met, aglstr, "akpw_foot_emiss")
if(length(str) != 0) new.trajpath <- file.path(new.trajpath, str)
if(length(str) != 0) new.intpath  <- file.path(new.intpath, str)
if(length(str) != 0) new.ncdfpath <- file.path(new.ncdfpath, str)

#------------------------------ STEP 3 --------------------------------------- #
#### METFILES and X-STILT configurations ***
if(trajecTF == TRUE){

  nhrs <- -72 # nhrs for generating trajec, when trajecTF=T
  if(selTF) nhrs <- -24 # 1day if trajecTF=T and selTF
  back.info <- run.backward.trajec(columnTF=columnTF, ocopath=ocopath,
                                   YYYYMMDD= substr(track.timestr,1,8), )
  orig.outname <-
  recp.info <-
  agl.info <-
  timestr <- unique(recp.info$timestr)

  # get metfile for generating backward trajectories
  metfile <- find.metfile(timestr=track.timestr, nhrs=nhrs, metpath=metpath,
                          met=met, trajecTF=trajecTF)

}else{

  ### if trajec exists, just for modeling XCO2
  # get trajec files, names and info, using ident.to.info()
  orig.outname <- list.files(pattern=filestr, path=orig.trajpath)
  recp.info <- ident.to.info(ident=orig.outname)[[1]]	# for receptor info
  agl.info <- ident.to.info(ident=orig.outname)[[2]]	# for release level info
  timestr <- unique(recp.info$timestr)

  # call get.grdhgt() to interpolate modeled ground heights from metfile,
  # add to original "recp.info"
  nhrs <- -1    # nhrs for trajwind() when trajecTF=F
  cat("Interpolating Ground Heights & winds @ receptors from met fields...\n")
  recp.info<-get.grdhgt(recp.info=recp.info,nummodel=nummodel,agl=10,nhrs=nhrs,
                        rundir=rundir, grdhgt.path=windpath, metpath=metpath,
                        metfile=metfile)

}

#------------------------------ STEP 4 --------------------------------------- #
#### INPUT X-STILT PARAMETERS FOR GENERATING FOOTPRINT
zbot <- 0
ztop <- 0			# Set to 0 if we want the surface influence volume
max.hour <- -nhrs       # change if not 3 days backward

# horizontal resolution of ODIAC
lon.res.anthro <- 1/120
lat.res.anthro <- lon.res.anthro

# horizontal resolution of CT-NRT
lon.res.bio <- 1
lat.res.bio <- lon.res.bio

# horizontal resolution of CT-NRT
lon.res.ocean <- 1
lat.res.ocean <- lon.res.ocean

#------------------------------ STEP 5 --------------------------------------- #
#### choose ODIAC version
odiac.ver   <- 3		# odiac version, 2 for v2016; 3 for v2017
odiac.vname <- c("2015a","2016","2017")[odiac.ver]
odiac.path  <- file.path(inpath,"anthro_inventories/ODIAC")
odiac.domain <- paste(min.lat,"-",max.lat,"N_",min.lon,"-",max.lon,"E",sep="")

## read ODIAC emissions to save time
if(hourlyTF){
	odiac.path <- file.path(odiac.path, odiac.vname, "hourly", site)
	odiac.file <- list.files(path=odiac.path, pattern=odiac.domain)
}else{
	odiac.path <- file.path(odiac.path, odiac.vname)
  pattern <- paste(substr(timestr,1,6),"_",odiac.domain,".nc",sep="")
  if(ppTF)pattern <- paste(substr(timestr,1,6),"_",odiac.domain,"_PP.nc",sep="")
	odiac.file <- list.files(path=odiac.path, pattern=pattern)
}

## create txt files for results and
vv<-"v5"  # avoid overwrite
txtpath <- workpath
txtname <- paste("STILT_OCO2_XCO2_",site,"_",timestr,"_",odiac.vname,"_",met,str,vv,".txt",sep="")
if(hourlyTF)txtname <- paste("STILT_OCO2_XCO2_",site,"_",timestr,"_",odiac.vname,"_hrly_",met,str,vv,".txt",sep="")
headers <- paste("xco2.", sel, sep="")  # headers for txtfile

## read in emissions, grab emission grid, [LAT, LON, (HR.back)]
# lat, lon should have already been converted to LOWER LEFT !!!
cat("take a while to readin ODIAC CO2 emissions...\n")
odiac.dat <- nc_open(file.path(odiac.path, odiac.file))
odiac.co2 <- ncvar_get(odiac.dat, "odiac_co2_emiss") # [lat,lon] or [lat,lon,hr]
odiac.lat <- ncvar_get(odiac.dat, "lat")
odiac.lon <- ncvar_get(odiac.dat, "lon")

if(hourlyTF){
	odiac.hr <- ncvar_get(odiac.dat, "hr")  # hours back are negative
	dimnames(odiac.co2) <- list(odiac.lat, odiac.lon, odiac.hr)
}else{
	dimnames(odiac.co2) <- list(odiac.lat, odiac.lon)
}

#### CHOOSE CT-NRT version
ct.version  <- c("v2016-1", "v2017")[2]
ctflux.path <- file.path(inpath, ct.version, "fluxes", "optimized")
ctmole.path <- file.path(inpath, ct.version, "molefractions")

#------------------------------ STEP 6 --------------------------------------- #
#### create a namelist (data.frame) including all variables
# and pass to subroutines
namelist<-list(met=met,site=site,timestr=timestr,filestr=filestr,sel=sel,
               min.lat=min.lat,max.lat=max.lat,min.lon=min.lon,max.lon=max.lon,
               ocopath=ocopath,txtpath=txtpath,txtname=txtname,headers=headers,
               orig.trajpath=orig.trajpath,new.trajpath=new.trajpath,
	             new.intpath=new.intpath,new.ncdfpath=new.ncdfpath,
               ct.version=ct.version,ctflux.path=ctflux.path,
               ctmole.path=ctmole.path,selTF=selTF, dmassTF=dmassTF,
               columnTF=columnTF,uncertTF=uncertTF,forwardTF=forwardTF,
               biascorrTF=biascorrTF, hourlyTF=hourlyTF, storeTF=storeTF,
							 ak.weight=ak.weight, pw.weight=pw.weight,
               lat.res.anthro=lat.res.anthro, lon.res.anthro=lon.res.anthro,
               lat.res.bio=lat.res.bio, lon.res.bio=lon.res.bio,
               lat.res.ocean=lat.res.ocean, lon.res.ocean=lon.res.ocean)

recp.info$ident<-orig.outname
trajlist<-list(recp.info,agl.info)

# store namelist and trajlist as RData file
#write.list(z=namelist,file=paste("namelist_",timestr,"_",site,".txt",sep=""))
write(x=namelist, file=paste("namelist_",timestr,"_",site,".txt",sep=""),sep="\n")

#------------------------------------------------ STEP 7 --------------------------------------------------------------- #
# start X-STILT runs, call oco2.get.xco2()
run.xstilt <- oco2.get.xco2(namelist=namelist,trajlist=trajlist,odiac.co2=odiac.co2)

# end of namelist
