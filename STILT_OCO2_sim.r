# Main script to simulate dCO2 given existing traject
# need existing traj, need ak, pw profiles from oco2
# Written by DW, 01/30/2017

# add simulation using trajec with error components, uncertTF=T, DW, 10/26/2017
# add simulation using trajec with bias corrections, biascorrTF=T, DW, 01/09/2018
# add GDAS 0.5deg simulation, DW, 01/24/2018

##### LATEST VERSION--
sourcepath<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/stiltR/"
source(paste(sourcepath,"sourceall.r",sep=""))
library(ncdf4)

#### 1-4 Riyadh, 5-7 Cairo, 8 for PRD
index<-1
met<-c("1km","gdas","gdas0p5")[2]
tt<-7
riyadh.timestr<-c("2014122710","2014122910","2015012810","2015081710","2015121610","2016011510","2016021610","2016072510")[tt]
cairo.timestr<-c("20150228","20150318","20160224")[tt]
site<-c("Riyadh","Cairo","PRD")[index]
timestr<-c(riyadh.timestr,cairo.timestr,"20150115")[index]  # locate day first, overpassing hour can be determined later from OCO-2 files
npar<-"3700P";if(site=="PRD")napr<-"2900P"
filestr<-paste(substr(timestr,1,4),"x",substr(timestr,5,6),"x",substr(timestr,7,8),sep="")

selTF<-F 			# true for only use 2-days trajs
debugTF<-F		# extended trajec beyond ZERO longitude
dmassTF<-T		# dmassTF when generating trajecfoot()
columnTF<-T   # whether a column receptor or fixed receptor
uncertTF<-F		# sim using trajec with transport error components
forwardTF<-F	# forward or backward traj
biascorrTF<-F # whether use bias-corrected trajec

nummodel<-10	# copy for trajwind()
if(debugTF)nummodel<-99    # use the fixed Fortran code (distribute particles crossing zero longitude)

str<-NULL
if(debugTF==TRUE)str<-"_debug"
if(uncertTF==TRUE)str<-"_uncert"
if(forwardTF==TRUE)str<-"_forward"
if(biascorrTF==TRUE)str<-"_bc"
if(columnTF==FALSE)aglstr<-"fixed_agl"
if(columnTF==TRUE)aglstr<-"multi_agl"

# for reruning ODAIC using TIMES, no need to run other components
# choose from "sel"
all<-c("anthro","bio","ocean","edp","apriori")
all.str<-c("dxco2.anthro","dxco2.bio","dxco2.ocean","xco2.edp","xco2.apriori")
sel<-all;headers<-all.str

odiac.ver<-3		# odiac version, 2 for v2016; 3 for v2017
odiac.vname<-c("2015a","2016","2017")[odiac.ver]
hourlyTF<-F	# true for using hourly ODIAC emissions

#######
# read in original trajs
orig.trajpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/Traj/",site,"/",met,"/",aglstr,"/orig_traj",str,"/",sep="")
if(length(str)==0)orig.trajpath<-paste(orig.trajpath,timestr,"/",sep="")
orig.outname<-list.files(pattern=filestr, path=orig.trajpath)
recp.info<-ident.to.info(ident=orig.outname)[[1]]
agl.info<-ident.to.info(ident=orig.outname)[[2]]
time.str<-unique(recp.info$time.str)

# path and filename for storing OCO-2 info
ocopath<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_input/OCO-2/OCO2_lite_b7rb/"

# path for storing ak.pw weighted traj
new.trajpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/Traj/",site,"/",met,"/",aglstr,"/akpw_traj",str,"/",sep="")
new.intpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/NetCDF/",site,"/",met,"/",aglstr,"/akpw_intfoot",str,"/",sep="")
new.ncdfpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/NetCDF/",site,"/",met,"/",aglstr,"/akpw_foot_emiss",str,"/",sep="")

# Meteological paths,"wrfout_d04_prd_201501.arl" for interpolating ground hgt, not running traj
metpath<-"/uufs/chpc.utah.edu/common/home/u0947337/"    	# path for the ARL format of WRF and GDAS
outpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/trajwind/",site,"/",sep="")
rundir<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/STILT_Exe/"  # Where to run STILT, where Copy lies

####### GRABBING ODIAC first
if(site=="PRD"& selTF==T)odiac.domain<-"15-40N_50-125E"
if(site=="PRD"& selTF==F)odiac.domain<-"10-50N_20-130E"
if(site=="Cairo" & debugTF==F)odiac.domain<-"0-50N_0-60E"
if(site=="Cairo" & debugTF==T)odiac.domain<-"15-60N_-40-40E"
#if(site=="Riyadh"& index==1)odiac.domain<-"5-30N_0-50E"
#if(site=="Riyadh"& index>1)odiac.domain<-"15-40N_0-50E"
if(site=="Riyadh")odiac.domain<-"0-50N_0-60E"

# create txt files for results and read ODIAC emissions to save time
vv<-"v5"
if(hourlyTF){
	filename<-paste("STILT_OCO2_XCO2_",site,"_",time.str,"_uneven_v",odiac.vname,"_",npar,"_hrly_",met,str,vv,".txt",sep="")
	odiac.path<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_input/anthro_inventories/ODIAC/v",odiac.vname,"/hourly/",site,"/",sep="")
	odiac.file<-paste("odiac",odiac.vname,"_1kmx1km_",substr(time.str,1,6),"_",unique(recp.info$recp.hour),"UTC_",odiac.domain,"_hrly.nc",sep="")
}else{
	filename<-paste("STILT_OCO2_XCO2_",site,"_",time.str,"_uneven_v",odiac.vname,"_",npar,"_",met,str,vv,".txt",sep="")
	odiac.path<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_input/anthro_inventories/ODIAC/v",odiac.vname,"/",sep="")
	odiac.file<-paste("odiac",odiac.vname,"_1kmx1km_",substr(time.str,1,6),"_",odiac.domain,".nc",sep="")
}
odiac.dat<-nc_open(paste(odiac.path,odiac.file,sep=""))

# grab emission grid, [LAT, LON, HR.back], lat, lon have already converted to lower left for this version
cat("take a while to readin ODIAC CO2 emissions...\n")
odiac.co2<-ncvar_get(odiac.dat, "odiac_co2_emiss")	# [lat,lon] or [lat,lon,hr]
odiac.lat<-ncvar_get(odiac.dat,"lat");odiac.lon<-ncvar_get(odiac.dat,"lon")
#hours back are negative
if(hourlyTF){odiac.hr<-ncvar_get(odiac.dat,"hr");dimnames(odiac.co2)<-list(odiac.lat,odiac.lon,odiac.hr)}else{dimnames(odiac.co2)<-list(odiac.lat,odiac.lon)}

### X-STILT parameters for generating footprints
# for fossil fuel emissions, from monthly mean ODIAC, res = 1/120deg, about 1km
# if using hourly ODAIC emission, generate hourly footprint,
# if not, generate time-integrated footprint for monthly mean ODIAC emissions

max.hour<-72;if(selTF)max.hour<-48
if(hourlyTF){foottimes<-seq(0,max.hour,1)}else{foottimes<-c(0,max.hour)}

zbot<-0;ztop<-0			# Set to 0 if we want the surface influence volume
lon.res<-1/120;lat.res<-1/120     # horizontal resolution of y grid

# lower left and x,y pix number of footprint to fit ODIAC
if(site=="Riyadh"){lon.ll<- 0;lat.ll<- 0;numpix.x<- 7200;numpix.y<- 6000}
if(site=="Cairo" & debugTF==F){lon.ll<- 0;lat.ll<- 0;numpix.x<- 7200;numpix.y<- 6000}
if(site=="Cairo" & debugTF==T){lon.ll<- -40;lat.ll<- 15;numpix.x<- 9600;numpix.y<- 5400}	# extended 3-day trajec for Cairo
if(site=="PRD" & selTF==T){lon.ll<-50;lat.ll<-15;numpix.x<-9000;numpix.y<-3000}	# 15-40N, 50-125E for all 2-day trajec
if(site=="PRD" & selTF==F){lon.ll<-20;lat.ll<-10;numpix.x<-13200;numpix.y<-4800}	# 10-50N, 20-130E for all 3-day trajec

# update DW, 09/15/2017
foottimes<-seq(0,72,1)
zbot<-0;ztop<-0					# Set to 0 if we want the surface influence volume
lon.res<-1;lat.res<-1     # horizontal resolution of y grid

# lower left and x,y pix number of footprint to fit CT-NRT
#if (site=="Riyadh"& index==1) {lon.ll<-0;lat.ll<-5;numpix.x<-50;numpix.y<-25}	# 5-30N, 0-50E for all 3-day trajec
#if (site=="Riyadh"& index>1) {lon.ll<-0;lat.ll<-15;numpix.x<-50;numpix.y<-25}	# 5-30N, 0-50E for all 3-day trajec
if (site=="Riyadh"){lon.ll<-0;lat.ll<-0;numpix.x<-60;numpix.y<-50}
if (site=="Cairo"& debugTF==F){lon.ll<-0;lat.ll<-0;numpix.x<-60;numpix.y<-50}
if (site=="Cairo"& debugTF==T){lon.ll<- -40;lat.ll<- 15;numpix.x<- 80;numpix.y<- 25}
if (site=="PRD"){lon.ll<-50;lat.ll<-15;numpix.x<-75;numpix.y<-25}	# 15-40N, 50-125E for all 2-day trajec

# if pass new.trajdat onto Trajecfoot(), no need to read RData again.

### carbontracker
ct.version<-"2016-1"
if(substr(ident,1,10)>="2016x01x01")ct.version<-"2017"
ctpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/group_data/CT-NRT_v",ct.version,"/fluxes/optimized/",sep="")

# add headers in txt files
headers<-c("soundingID", "lat", "lon",headers)
if(length(sel)>=3)headers<-c(headers,"obs.total.xco2","sim.total.xco2", "xco2.err")
write(headers, file = filename, ncolumns = length(headers), append = TRUE, sep = ",")

# loop over existing traj
cat(paste("Start OCO-2 simulations for",sel,"...\n"))
endfile<-length(orig.outname)
if(uncertTF)endfile<-length(orig.outname)/2

#### loop over all soundings
for (i in 1:endfile){

  cat("#-------- Working on receptor/sounding #",i,"---", signif(i/length(orig.outname)*100,3),"%-------- #\n")

	# initialize total XCO2,	#dxco2.anthro + dxco2.bio + dxco2.ocean + xco2.edp + xco2.prior
	total.sim.xco2<-0
	if(uncertTF){
		ident<-orig.outname[(2*i-1):(2*i)]
		tmp.recp.info<-recp.info[(2*i-1):(2*i),]
		tmp.agl.info<-agl.info[(2*i-1):(2*i),]
		tmp.lat<-unique(tmp.recp.info$recp.lat)
		tmp.lon<-unique(tmp.recp.info$recp.lon)
	}else{
		ident<-orig.outname[i]
		tmp.recp.info<-recp.info[i,]
		tmp.agl.info<-agl.info[i,]
		tmp.lat<-tmp.recp.info$recp.lat
		tmp.lon<-tmp.recp.info$recp.lon
	}

	#### call main subroutine to get each contributions...
	result<-oco2.get.xco2(ident=ident,sel=sel,site=site,timestr=timestr,ocopath=ocopath,recp.lat=tmp.lat,recp.lon=tmp.lon,
												recp.info=tmp.recp.info,agl.info=agl.info,orig.trajpath=orig.trajpath,new.trajpath=new.trajpath,
												ak.weight=T,pw.weight=T,min.lon=min.lon,max.lon=max.lon,min.lat=min.lat,max.lat=max.lat,foottimes=c(0,72),
												ct.version="2017",odiac.co2=odiac.co2,dmassTF=T,storeTF=F,new.intpath=new.intpath,new.ncdfpath=new.ncdfpath)

	write(result, file=filename, ncolumns = length(headers),append=TRUE, sep=",")	# TRUE, keep adding lines without overwriting
	gc()	# free memory
	cat("\n")
} # end looping over all soundings

# end of main script
