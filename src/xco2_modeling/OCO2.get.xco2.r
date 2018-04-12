### Topmost subroutine to simulate column-averaged CO2 dry mole fractions using X-STILT
#
#---------------------------------------------- STEPS involved -------------------------------------------------------------------- #
# 1) interpolate modeled ground heights based on meteorological fields;
#	2) read in .RData files for trajec,
#	   interpolate satellite profiles (averaging kernel, pressure weighting and prior profiles), return "combine.profile",
#    call "weight.trajecfootv2()" to weight trajec-level footprints using interpolated satellite profiles;
# 3) call "get.foot()" and "odiac.anthro()" to generate 2D column footprint and get anthropogenic XCO2 enhancements;
# 4) call "ctnrt.biov2()" to get biospheric XCO2 changes;
# 5) call "ctnrt.oceanv2()" to get oceanic XCO2 changes;
# 6) call "ctnrt.background()" to get XCO2 boundary conditions based on CarbonTracker global fields;
# 7) call "oco2.apriori()" to get prior portion

####---------------------------------------- Input variables include ---------------------------------------------------------------------------- #
# ident: name for STILT trajectory files, including ".RData", e.g., "2014x12x27x"
# sel: XCO2 contributions needed to be modeled, can be subset of c("anthro","bio","ocean","edp","apriori"),
#      where "anthro","bio","ocean" stand for CO2 sources/sinks, "edp" for CO2 boundary condition using CarbonTracker-NearRealTime
#      and "apriori" for XCO2 contribution from OCO-2 prior profiles, i.e., (I-AK)*PW*CO2.ap;
# site: target city, e.g., "Riyadh";
# timestr: time of satellite overpass, e.g., "";
# ocopath: path that stores OCO-2 Lite files
# recp.lat & recp.lon: receptor lat/lon, same as selected sounding lat/lon;
# recp.info & agl.info: output from ident.to.info(), including receptor time, lat, lon (recp.info) and release levels (agl.info)
# orig.trajpath & new.trajpath: paths that stores initial trajectories & weighted trajectories (e.g., by averaging kernel or pressure weighting)
# ak.weight & pw.weight: logical flags for whether weighting against averaging kernel or pressure weighting, can be T or F
# min.lon,max.lon,min.lat,max.lat: LOWER LEFT corners of spatial domains for generating footprints,
# 																 these variables will determine the numpix.x, numpix.y, lon.ll, lat.ll;
# foottimes: vector forms for footprint duration, e.g., c(0,72) meaning time-integrated footprint up to 3 days back, same as foottimes in Trajecfoot()
# ct.version: version for CarbonTracker, required when simulating bio and oceanic contributions
# odiac.co2: matrix for anthropogenic CO2 emissions from ODIAC, always readin emissions before call this function
#						 fixed dimension: [Lat, Lon,/Hours], fixed units: micromole/m2/s
# dmassTF: flags for whether to turn on mass conservation, check Trajecfoot() for details
# storeTF: flags for whether to store time-integrated footprint and CO2 contributions
# new.intpath, new.ncdfpath: paths for storing time-integrated footprint and CO2 contributions fields
#
#### Written by Dien Wu, 10/20/2017

#-------------------------------------------- Updates --------------------------------------------------------------------------- #
# use ODIACv2016 and CT-NRTv2016-1 as default input dataset
# ADD biospheric CO2 from CT-NRT, 02/10/2017
# ADD ODIACv2016, 02/15/2017
# ADD hourly ODIAC emissions, 03/09/2017
# Get ak, pw profiles interpolation code separated by weight trajecfoot code
# Add grabbing upper background profiles from CT-NRT
# ADD oceanic fluxes from CT-NRT, 07/14/2017
# update two subroutines for bio and ocean, 09/15/2017
# add ak.weight==F, thus, no need to weight through AK and apriori is no longer needed, DW, 02/06/2018
#######################################################################################################################

oco2.get.xco2<-function(ident,sel,site,timestr,ocopath,recp.lat,recp.lon,recp.info,agl.info,orig.trajpath,new.trajpath,
												ak.weight=T,pw.weight=T,min.lon,max.lon,min.lat,max.lat,foottimes=c(0,72),ct.version="2017",
												odiac.co2,dmassTF=T,storeTF=F,new.intpath,new.ncdfpath){

	#------------------------------------------------------------------------------------------------------------------------------ #
	# 1. get OCO2 info from traj info, requires get.oco2info() from "OCO2.get.oco2info.r"
	ocofile<-list.files(pattern=as.character(timestr),path=ocopath)
	oco2.info<-get.oco2info(ocopath=ocopath, ocofile=ocofile, recp.lat=recp.lat, recp.lon=recp.lon)
	oco2.profiles<-oco2.info[[1]]

	# initialize result
	id<-as.character(oco2.info[[2]]$find.id)	# OCO-2 ID
	result<-c(id, recp.lat, recp.lon)	# store ID, recptor/sounding lat/lon

	# WRF or GDAS path and files for interpolating the modeled grdhgt from metfile
	cat("oco2.get.xco2():Interpolating Ground Height from met fields...\n")
	recp.grdhgt<-get.grdhgt(met=met,site=site,recp.info=recp.info,nummodel=nummodel)
	cat(paste("Interpolated Modeled Ground Height:",recp.grdhgt,"meters...\n"))

  #------------------------------------------------------------------------------------------------------------------------------ #
  # 2. read in traj and apply ak and pw profiles onto original trajs, #orig.trajdat<-getr.old(orig.outname[i],path=orig.trajpath)
	ident<-substr(ident,1,nchar(ident)-6)
  orig.trajdat<-getr(ident,path=orig.trajpath)
	# selTF for only use 2-days trajs
	#if(selTF){sel.trajdat<-orig.trajdat[orig.trajdat[,"time"]>= nhrs*60,];orig.trajdat<-sel.trajdat}

	### start weighting trajec-level footprint...call weight.trajecfootv2()
	# new.trajpath: path for storing weighted trajectory
	cat("oco2.get.xco2():Start weighting trajec-level footprint based on satellite profiles...\n")
	adjust.profilev2<-weight.trajecfootv2(ident=ident, trajdat=orig.trajdat, recp.info=recp.info, agl.info=agl.info,
                                  			oco.ak.norm=oco2.profiles$ak.norm, oco.pw=oco2.profiles$pw, oco.pres=oco2.profiles$pres, oco.apriori=oco2.profiles$apriori,
                                  			recp.grdhgt=recp.grdhgt, new.trajpath=new.trajpath, ak.weight=ak.weight, pw.weight=pw.weight)

  ## returns the ak pw profiles at for all levels & the weighted footprint
	combine.profile<-adjust.profilev2[[1]]	# all profile info, ak, pw, apriori
	new.trajdat<-adjust.profilev2[[2]]	# weighted traj still in matrix form

	#### initialize total modeled XCO2
	total.sim.xco2<-0

	#------------------------------------------------ 3. ANTHROPOGENIC  ------------------------------------------------------------------------------ #
	if(!is.na(match("anthro",sel))){

		# 3.1. GENERATE 2D footprint based on the above weighted traj, call get.foot() from "OCO2.get.foot.r"
		# if there are only 2 components in "foottimes",we need the time-integrated 2D footprint (xfoot), not time-varying 3D footprint(foot)
		lon.res<-
		lat.res<-
		numpix.x<-(max.lon-min.lon)/lon.res+1
		numpix.y<-(max.lat-min.lat)/lat.res+1

		cat("oco2.get.xco2():Generating footprint that matches anthropogenic emissions...\n")
		foot.anthro<-get.foot(ident=ident,foot.overwrite=TRUE,part=new.trajdat,trajpath=new.trajpath,footpath=new.intpath,
													foottimes=foottimes,zlim=c(0,0),dmassTF=dmassTF,fluxweighting=NULL,coarse=1,
													numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=min.lon,lat.ll=min.lat,lon.res=lon.res,lat.res=lat.res,storeTF=storeTF)

		# 3.2. CALL the subroutine for matching ODIAC emission with xfoot, and sum the matrix to obtain "dxco2.anthro"
		cat("oco2.get.xco2():Calculating the dCO2 from anthro using ODIAC...\n")
		dxco2.anthro<-odiac.anthro(foot=foot.anthro,odiac.co2=odiac.co2,ident=ident,storeTF=storeTF,ncdfpath=new.ncdfpath) 	# store the foot*emission in newncdfpath
		cat(paste("Column dCO2 from ODIAC emission: ", signif(dxco2.anthro, 4), "ppm\n\n"))

		# store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+dxco2.anthro
		result<-c(result,dxco2.anthro)
	}	# end "anthro" contribution

	#------------------------------------------------ 4. BIOSPERHIC ------------------------------------------------------------------------------ #
	if(!is.na(match("bio",sel))){

		# 4.1. GENERATE 2D footprint based on the above weighted traj, call get.foot() from "OCO2.get.foot.r"
		# no need to store 1deg footprint
		# if there are more than 2 components in "foottimes", meaning we need the time-varying 3D footprint, not integrated footprint
		# update DW, 04/11/2018
		lon.res<-1
		lat.res<-1
		numpix.x<-(max.lon-min.lon)/lon.res+1
		numpix.y<-(max.lat-min.lat)/lat.res+1

		cat("oco2.get.xco2():Generating footprint that matches CT-NRT bio...\n")
		foot.bio<-get.foot(ident=ident,foot.overwrite=TRUE,part=new.trajdat,trajpath=new.trajpath,footpath=new.intpath,
											 foottimes=foottimes,zlim=c(0,0),dmassTF=dmassTF,fluxweighting=NULL,coarse=1,
											 numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res,storeTF=storeTF)

		# 4.2. CALL the subroutine for matching biosperhic fluxes with footprint, and sum the matrix to obtain "dxco2.bio"
		cat("oco2.get.xco2():Calculating the dCO2 from bio using CT-NRT...\n")
		dxco2.bio<-ctnrt.biov2(ident=ident, foot=foot.bio, ct.version=ct.version,storeTF=storeTF, ncdfpath=new.ncdfpath) 	# store the foot*emission in newncdfpath
		cat(paste("Column dCO2 from biosperhic exchange: ", signif(dxco2.bio, 4), "ppm\n\n"))
		# if foot hours differ from CT hours, the subroutine will return NA values

		# 4.3. store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+dxco2.bio
		result<-c(result,dxco2.bio)
	}	# end "bio" sel

	#------------------------------------------------ 5. OCEANIC ------------------------------------------------------------------------------ #
	if(!is.na(match("ocean",sel))){

		# 5.1. GENERATE 2D footprint based on the above newly weighted traj, using "get.intfoot.r"
		# no need to store 1deg footprint
		# if there are more than 2 components in "foottimes", meaning we need the time-varying 3D footprint, not integrated footprint
		cat("oco2.get.xco2():Generating footprint that matches CT-NRT ocean...\n")

		# if we already have biospheric footprint (same as oceanic footprint)
		if(!is.na(match("bio",sel)*match("ocean",sel))){
			foot.ocean<-foot.bio
		}else{
			# update DW, 04/11/2018
			lon.res<-1
			lat.res<-1
			numpix.x<-(max.lon-min.lon)/lon.res+1
			numpix.y<-(max.lat-min.lat)/lat.res+1
			foot.ocean<-get.foot(ident=ident, foot.overwrite=TRUE,part=new.trajdat, trajpath=new.trajpath,footpath=new.intpath,
													 foottimes=foottimes,zlim=c(0,0),fluxweighting=NULL,coarse=1,
													 numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=min.lon,lat.ll=min.lat,lon.res=lon.res,lat.res=lat.res, storeTF=storeTF)
		}

		# 5.2. CALL the subroutine for matching biosperhic fluxes with footprint, and sum the matrix to obtain "dxco2.bio"
		cat("oco2.get.xco2():Calculating the dCO2 from oceanic fluxes using CT-NRT...\n")
		dxco2.ocean<-ctnrt.oceanv2(ident=ident, foot=foot.ocean, storeTF=storeTF, res=lon.res, ncdfpath=new.ncdfpath) 	# store the foot*emission in newncdfpath
		cat(paste("Column dCO2 from oceanic fluxes: ", signif(dxco2.ocean, 4), "ppm\n\n"))
		# if foot hours differ from CT hours, the subroutine will return NA values

		# 5.3. store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+dxco2.ocean
		result<-c(result,dxco2.ocean)

	}	# end "ocean" sel

	#------------------------------------------------ 6. ENDPOINTS for boundary conditions ------------------------------------------------------------------------------ #
	if(!is.na(match("edp",sel))){

		# 6.1. CALL the subroutine for calculating background CO2 (AK PW weighted CO2 concentration derived from CT-NRT at endpoints)
		# sum all to obtain "XCO2.edp", remember to use weighted traj "new.trajdat"
		edp.info<-ctnrt.background(ident=ident, trajdat=new.trajdat, combine.profile=combine.profile)
		xco2.edp<-edp.info[[2]]
		combine.profile<-edp.info[[1]]
		cat(paste("Column CO2 background using CT-NRT: ", signif(xco2.edp,5), "ppm\n\n"))

		# 6.2. store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+xco2.edp
		result<-c(result,xco2.edp)
	}	# end "edp" sel

	#------------------------------------------------ 7. APRIORI contributions from OCO-2 ------------------------------------------------------------------------------ #
	if(!is.na(match("apriori",sel))){

		# 7.1 CALL the subroutine for calculating the contribution of CO2 from a priori profile (I - AK), sum all to obtain "XCO2.apriori"
		xco2.prior<-oco2.apriori(combine.profile=combine.profile)
		cat(paste("Contribution from CO2 apriori profiles:", signif(xco2.prior,5),"ppm\n\n"))

		# 7.2. store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+xco2.prior
		result<-c(result,xco2.prior)
	}	# end "apriori" sel

	######
  #-------------------------------------------------- 8. calculate the TOTAL modeled XCO2 ---------------------------------------------------------------------------- #
	if(length(sel)>=3){
		# 8.1. GRAB observed XCO2 for each sounding
		sel.xco2.oco<-oco2.info[[2]]$xco2.oco
		cat(paste("TOTAL retrieved COLUMN CO2:", signif(sel.xco2.oco,5), "ppm\n"))

		# 8.2. SUM all contributions to get the final total simulated retrieved XCO2.ak, which can be compared with retrieved XCO2 from OCO-2
		cat(paste("TOTAL simulated COLUMN CO2:", signif(total.sim.xco2,5), "ppm\n"))
		diff.xco2<-total.sim.xco2-sel.xco2.oco	# sim - obs
		# end calculating the total and model-data mismatch

		# 8.3. add results in a txt file, i.e., OCO-2 observed XCO2, diff between observed and simulated XCO2
		result<-c(result,sel.xco2.oco,total.sim.xco2,diff.xco2)
	}

	return(result)
}
