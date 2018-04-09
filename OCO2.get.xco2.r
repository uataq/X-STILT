# Topmost subroutine to simulate column-averaged CO2 dry mole fractions using X-STILT
# steps involved: 1) interpolate modeled ground heights based on meteorological fields
#									2)

# Input variables include:
# ident: name for STILT

# use ODIACv2016 and CT-NRTv2016-1 as default input dataset
# Written by Dien Wu, 10/20/2017

# add ak.weight==F, thus, no need to weight through AK and apriori is no longer needed, DW, 02/06/2018
oco2.get.xco2<-function(ident,sel,site,new.trajdat,combine.profile,new.trajpath, new.intpath, new.ncdfpath,
												debugTF,dmassTF=T,storeTF,selTF=F,nhrs=NULL,hourlyTF=F,odiac.co2,odiac.domain,ak.weight){

	#------------------------------------------------------------------------------------------------------------------------------ #
  # 1. get OCO2 info from traj info,
	#oco2.info<-get.oco2info(ocopath=ocopath, ocofile=ocofile, recp.lat=24.0233, recp.lon=47.3054)
	ocofile<-list.files(pattern=as.character(time.str),path=ocopath)
  oco2.info<-get.oco2info(ocopath=ocopath, ocofile=ocofile, recp.lat=tmp.lat, recp.lon=tmp.lon)
  oco2.profiles<-oco2.info[[1]]

	# initialize result
	id<-as.character(oco2.info[[2]]$find.id)
	tmp<-c(id, tmp.recp.info$recp.lat, tmp.recp.info$recp.lon)

  # WRF or GDAS path and files for interpolating the modeled grdhgt from metfile
	cat("Interpolating Ground Height from met fields...\n")
	recp.grdhgt<-get.grdhgt(met=met,site=site,recp.info=tmp.recp.info,nummodel=nummodel)
	print(recp.grdhgt)	# in meter

  #------------------------------------------------------------------------------------------------------------------------------ #
  # 2. read in traj and apply ak and pw profiles onto original trajs, #orig.trajdat<-getr.old(orig.outname[i],path=orig.trajpath)
	ident<-substr(ident,1,nchar(ident)-6)
  orig.trajdat<-getr(ident,path=orig.trajpath)

	# selTF for only use 2-days trajs
	if(selTF){sel.trajdat<-orig.trajdat[orig.trajdat[,"time"]>= nhrs*60,];orig.trajdat<-sel.trajdat}
	adjust.profilev2<-weight.trajecfootv2(ident=ident, trajdat=orig.trajdat, recp.info=tmp.recp.info, agl.info=tmp.agl.info,
                                  			oco.ak.norm=oco2.profiles$ak.norm, oco.pw=oco2.profiles$pw, oco.pres=oco2.profiles$pres, oco.apriori=oco2.profiles$apriori,
                                  			recp.grdhgt=recp.grdhgt, new.trajpath=new.trajpath, ak.weight=TRUE, pw.weight=TRUE)

  # returns the ak pw profiles at for all levels & the weighted footprint
	combine.profile<-adjust.profilev2[[1]]	# all profile info, ak, pw, apriori
	new.trajdat<-adjust.profilev2[[2]]	# weighted traj still in matrix form

	#### initialize total modeled XCO2 and result
	total.sim.xco2<-0
	result<-NULL

	#------------------------------------------------ 3. ANTHROPOGENIC  ------------------------------------------------------------------------------ #
	if(!is.na(match("anthro",sel))){

		# 1. GENERATE 2D footprint based on the above newly weighted traj, using "get.intfoot.r"
		# for fossil fuel emissions, from monthly mean ODIAC, res = 1/120deg, about 1km
		# if using hourly ODAIC emission, generate hourly footprint,
		# if not, generate time-integrated footprint for monthly mean ODIAC emissions
		max.hour<-72;if(selTF)max.hour<-48
		if(hourlyTF){foottimes<-seq(0,max.hour,1)}else{foottimes<-c(0,max.hour)}

		zbot<-0;ztop<-0			# Set to 0 if we want the surface influence volume
		lon.res<-1/120;lat.res<-1/120     # horizontal resolution of y grid

		# lower left and x,y pix number of footprint to fit ODIAC
		#if(site=="Riyadh"&index==1){lon.ll<-0;lat.ll<-5;numpix.x<-6000;numpix.y<-3000}	# 5-30N, 0-50E for all 3-day trajec
		#if(site=="Riyadh"&index>1){lon.ll<-0;lat.ll<-15;numpix.x<-6000;numpix.y<-3000}	# 5-30N, 0-50E for all 3-day trajec
		if(site=="Riyadh"){lon.ll<- 0;lat.ll<- 0;numpix.x<- 7200;numpix.y<- 6000}
		if(site=="Cairo" & debugTF==F){lon.ll<- 0;lat.ll<- 0;numpix.x<- 7200;numpix.y<- 6000}
		if(site=="Cairo" & debugTF==T){lon.ll<- -40;lat.ll<- 15;numpix.x<- 9600;numpix.y<- 5400}	# extended 3-day trajec for Cairo
		if(site=="PRD" & selTF==T){lon.ll<-50;lat.ll<-15;numpix.x<-9000;numpix.y<-3000}	# 15-40N, 50-125E for all 2-day trajec
		if(site=="PRD" & selTF==F){lon.ll<-20;lat.ll<-10;numpix.x<-13200;numpix.y<-4800}	# 10-50N, 20-130E for all 3-day trajec

		# if pass new.trajdat onto Trajecfoot(), no need to read RData again.
		# if there are only 2 components in "foottimes",we need the time-integrated 2D footprint (xfoot), not time-varying 3D footprint(foot)
		cat("oco2.get.xco2():Generating footprint that matches ODIAC...\n")

		foot.anthro<-get.foot(ident=ident, foot.overwrite=TRUE,part=new.trajdat, trajpath=new.trajpath, footpath=new.intpath, foottimes=foottimes,zlim=c(zbot,ztop),
	                        dmassTF=dmassTF,fluxweighting=NULL,coarse=1,numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res, storeTF=storeTF)

		# 2. CALL the subroutine for matching ODIAC emission with xfoot, and sum the matrix to obtain "dxco2.anthro"
		# odiac.ver = 1 for ODIACv2015a; odiac.ver = 2 for ODIACv2016
		# if using hourly ODIAC, need to readin emission grid to save time...
		cat("oco2.get.xco2():Calculating the dCO2 from anthro using ODIAC...\n")
		dxco2.anthro<-odiac.anthro(ident=ident, odiac.co2=odiac.co2, odiac.domain=odiac.domain, site=site, foot=foot.anthro, odiac.ver=2, storeTF=storeTF, res=lon.res, ncdfpath=new.ncdfpath) 	# store the foot*emission in newncdfpath
		cat(paste("Column dCO2 from ODIAC emission: ", signif(dxco2.anthro, 4), "ppm\n\n"))

		# store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+dxco2.anthro
		result<-c(result,dxco2.anthro)
	}	# end "anthro" contribution

	#------------------------------------------------ 4. BIOSPERHIC ------------------------------------------------------------------------------ #
	if(!is.na(match("bio",sel))){

		# 3. GENERATE 2D footprint based on the above newly weighted traj, using "get.intfoot.r"
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
		# no need to store 1deg footprint
		# if there are more than 2 components in "foottimes", meaning we need the time-varying 3D footprint, not integrated footprint
		cat("oco2.get.xco2():Generating footprint that matches CT-NRT bio...\n")
		foot.bio<-get.foot(ident=ident, foot.overwrite=TRUE,part=new.trajdat, trajpath=new.trajpath, footpath=new.intpath, foottimes=foottimes,zlim=c(zbot,ztop),
	                     dmassTF=dmassTF,fluxweighting=NULL,coarse=1,numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res, storeTF=storeTF)

		# 4. CALL the subroutine for matching biosperhic fluxes with footprint, and sum the matrix to obtain "dxco2.bio"
		cat("oco2.get.xco2():Calculating the dCO2 from bio using CT-NRT...\n")
		dxco2.bio<-ctnrt.biov2(ident=ident, foot=foot.bio, storeTF=storeTF, res=lon.res, ncdfpath=new.ncdfpath) 	# store the foot*emission in newncdfpath
		cat(paste("Column dCO2 from biosperhic exchange: ", signif(dxco2.bio, 4), "ppm\n\n"))
		# if foot hours differ from CT hours, the subroutine will return NA values

		# store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+dxco2.bio
		result<-c(result,dxco2.bio)

	}	# end "bio" sel

	#------------------------------------------------ 5. OCEANIC ------------------------------------------------------------------------------ #
	if(!is.na(match("ocean",sel))){

		# 5. GENERATE 2D footprint based on the above newly weighted traj, using "get.intfoot.r"
		# for oceanic CO2 fluxes, from 3 hourly CT-NRT_v2016, res = 1 deg
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
		# no need to store 1deg footprint
		# if there are more than 2 components in "foottimes", meaning we need the time-varying 3D footprint, not integrated footprint
		cat("oco2.get.xco2():Generating footprint that matches CT-NRT ocean...\n")
		if(!is.na(match("bio",sel)*match("ocean",sel))){	# if we already have biospheric footprint (same as oceanic footprint)
			foot.ocean<-foot.bio
		}else{
			foot.ocean<-get.foot(ident=ident, foot.overwrite=TRUE,part=new.trajdat, trajpath=new.trajpath, footpath=new.intpath, foottimes=foottimes,zlim=c(zbot,ztop),
																fluxweighting=NULL,coarse=1,numpix.x=numpix.x,numpix.y=numpix.y,lon.ll=lon.ll,lat.ll=lat.ll,lon.res=lon.res,lat.res=lat.res, storeTF=storeTF)
		}

		# 6. CALL the subroutine for matching biosperhic fluxes with footprint, and sum the matrix to obtain "dxco2.bio"
		cat("oco2.get.xco2():Calculating the dCO2 from oceanic fluxes using CT-NRT...\n")
		dxco2.ocean<-ctnrt.oceanv2(ident=ident, foot=foot.ocean, storeTF=storeTF, res=lon.res, ncdfpath=new.ncdfpath) 	# store the foot*emission in newncdfpath
		cat(paste("Column dCO2 from oceanic fluxes: ", signif(dxco2.ocean, 4), "ppm\n\n"))
		# if foot hours differ from CT hours, the subroutine will return NA values

		# store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+dxco2.ocean
		result<-c(result,dxco2.ocean)

	}	# end "ocean" sel

	######
	#------------------------------------------------ 6. ENDPOINTS for boundary conditions ------------------------------------------------------------------------------ #
	if(!is.na(match("edp",sel))){

		# 7. CALL the subroutine for calculating background CO2 (AK PW weighted CO2 concentration derived from CT-NRT at endpoints)
		# sum all to obtain "XCO2.edp", remember to use weighted traj "new.trajdat"
		edp.info<-ctnrt.background(ident=ident, trajdat=new.trajdat, combine.profile=combine.profile)
		xco2.edp<-edp.info[[2]]
		combine.profile<-edp.info[[1]]
		cat(paste("Column CO2 background using CT-NRT: ", signif(xco2.edp,5), "ppm\n\n"))

		# store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+xco2.edp
		result<-c(result,xco2.edp)

	}	# end "edp" sel

	######
	#------------------------------------------------ 7. APRIORI contributions from OCO-2 ------------------------------------------------------------------------------ #
	if(!is.na(match("apriori",sel))){

		# 8. CALL the subroutine for calculating the contribution of CO2 from a priori profile (I - AK), sum all to obtain "XCO2.apriori"
		xco2.prior<-oco2.apriori(combine.profile=combine.profile)
		cat(paste("Contribution from CO2 apriori profiles:", signif(xco2.prior,5),"ppm\n\n"))

		# store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+xco2.prior
		result<-c(result,xco2.prior)

	}	# end "apriori" sel

	######
  #-------------------------------------------------- 8. calculate the TOTAL modeled XCO2 ---------------------------------------------------------------------------- #
	if(length(sel)>=3){

		# 9. GRAB observed XCO2 for each sounding
		sel.xco2.oco<-oco2.info[[2]]$xco2.oco
		cat(paste("TOTAL retrieved COLUMN CO2:", signif(sel.xco2.oco,5), "ppm\n"))

		# SUM all contributions to get the final total simulated retrieved XCO2.ak, which can be compared with retrieved XCO2 from OCO-2
		cat(paste("TOTAL simulated COLUMN CO2:", signif(total.sim.xco2,5), "ppm\n"))
		diff.xco2<-total.sim.xco2-sel.xco2.oco	# sim - obs
		# end calculating the total and model-data mismatch

		# add results in a txt file, i.e., OCO-2 observed XCO2, diff between observed and simulated XCO2
		result<-c(result,sel.xco2.oco,total.sim.xco2,diff.xco2)

	}

return(result)
}
