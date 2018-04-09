# Topmost subroutine to model XCO2 in X-STILT
# use ODIACv2016 and CT-NRTv2016-1 as default input dataset
# DW, 10/20/2017

# add ak.weight==F, thus, no need to weight through AK and apriori is no longer needed, DW, 02/06/2018
oco2.get.xco2<-function(ident=ident,sel=sel,site=site, new.trajdat=new.trajdat, combine.profile=combine.profile,
												new.trajpath=new.trajpath, new.intpath=new.intpath, new.ncdfpath=new.ncdfpath,
												debugTF=debugTF, dmassTF=T, storeTF=storeTF,selTF=F, hourlyTF=F, odiac.co2=odiac.co2,odiac.domain=odiac.domain,ak.weight=akTF){

	total.sim.xco2<-0
	result<-NULL

	#------------------------------------------------ ANTHROPOGENIC  ------------------------------------------------------------------------------ #
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

	#------------------------------------------------ BIOSPERHIC ------------------------------------------------------------------------------ #
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

				# for plotting footprint
				# flip footprint
				plotTF<-FALSE
				if(plotTF){
					library(reshape);library(ggplot2)
					melt.foot<-melt(foot.bio)#;colnames(melt.foot)<-c("lat","lon","foot")
					m1<-ggplot.map(minlat=10,maxlat=50,minlon=10,maxlon=50)
					f1<-m1+geom_point(data=melt.foot[melt.foot$value!=0 &melt.foot$X3<60,],aes(x=X2,y=X1,colour=value))
					f2<-ggplot()+geom_point(data=melt.foot[melt.foot$value!=0 &melt.foot$X3<48,],aes(x=X2,y=X1,colour=value))
					#melt.foot<-melt.foot[melt.foot$value>signal,]

					xfoot<-apply(foot.bio,c(1,2),sum)
					#f.foot.bio<-aperm(foot.bio, c(2,1,3))
					foot.lat<-as.numeric(dimnames(xfoot)[[1]]);foot.lon<-as.numeric(dimnames(xfoot)[[2]])
					library(fields);library(RColorBrewer)
					cols<-colorRampPalette(brewer.pal(9,"Blues"))(100)
					image.plot(foot.lon, foot.lat, t(xfoot),col=cols)
					image.plot(foot.lon, foot.lat, t(foot.bio[,,10]),col=cols)
				}

				# 4. CALL the subroutine for matching biosperhic fluxes with footprint, and sum the matrix to obtain "dxco2.bio"
				cat("oco2.get.xco2():Calculating the dCO2 from bio using CT-NRT...\n")
				dxco2.bio<-ctnrt.biov2(ident=ident, foot=foot.bio, storeTF=storeTF, res=lon.res, ncdfpath=new.ncdfpath) 	# store the foot*emission in newncdfpath
				cat(paste("Column dCO2 from biosperhic exchange: ", signif(dxco2.bio, 4), "ppm\n\n"))
				# if foot hours differ from CT hours, the subroutine will return NA values

				# store result and sum up total XCO2
				total.sim.xco2<-total.sim.xco2+dxco2.bio
				result<-c(result,dxco2.bio)
	}	# end "bio" sel

	#------------------------------------------------ OCEANIC ------------------------------------------------------------------------------ #
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

	#------------------------------------------------ ENDPOINTS ------------------------------------------------------------------------------ #
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

	#------------------------------------------------ APRIORI ------------------------------------------------------------------------------ #
	if(!is.na(match("apriori",sel))){

			# 8. CALL the subroutine for calculating the contribution of CO2 from a priori profile (I - AK), sum all to obtain "XCO2.apriori"
			xco2.prior<-oco2.apriori(combine.profile=combine.profile)
			cat(paste("Contribution from CO2 apriori profiles:", signif(xco2.prior,5),"ppm\n\n"))

			# store result and sum up total XCO2
			total.sim.xco2<-total.sim.xco2+xco2.prior
			result<-c(result,xco2.prior)

	}	# end "apriori" sel

	# check all profiles
	if(F){

			#x11(width=12, height=6)
			png(filename=paste("./all_profiles_",recp.lat[i],"N.png",sep=""),width=1200, height=720)
			par(mfrow=c(1,2))
			par(mar=c(5,5,5,2))
			cex<-2
			#plot(combine.profile$combine.ak.norm, combine.profile$combine.pres, ylim=c(1000,0), pch=19, xlim=c(0,1.2), main="Normalized AK profiles\nfor combined (BLACK) and STILT (RED) levels")
			plot(oco2.profiles$ak.norm, oco2.profiles$pres, ylim=c(1000,0), pch=20, xlim=c(0,1.2), main=expression(atop("Normalized averaging kernel profiles","at OCO-2 vs. Model levels")),xlab="Averaging Kernel (AK)",ylab="Pressure [hPa]",cex=cex,cex.axis=cex, cex.lab=cex, cex.main=cex)
			lines(oco2.profiles$ak.norm, oco2.profiles$pres,lty=2)
			points(combine.profile$ak.norm[combine.profile$stiltTF], combine.profile$pres[combine.profile$stiltTF], col="red", pch=1,cex=cex)
			points(combine.profile$ak.norm[combine.profile$stiltTF==F], combine.profile$pres[combine.profile$stiltTF==F], col="blue", pch=1,cex=cex)
			#abline(h=combine.profile[31,"pres"],lty=3,col="red",cex=cex);abline(h=combine.profile[37,"pres"],lty=3,col="red",cex=cex);abline(h=combine.profile[nrow(combine.profile),"pres"],lty=3,col="red",cex=cex)
			abline(h=combine.profile[1,"pres"],lty=3,col="red",cex=cex);abline(h=combine.profile[37,"pres"],lty=3,col="red",cex=cex)
			abline(h=combine.profile[38,"pres"],lty=2,col="blue",cex=cex);abline(h=combine.profile[nrow(combine.profile),"pres"],lty=2,col="blue",cex=cex)
			legend(0.02, 800,c("OCO-2","STILT","CT-NRT"),col=c("black","red","blue"),text.col=c("black","red","blue"),pch=c(20, 1,1),lty=c(2,NA,NA),bty="n",cex=cex-0.3)

			if(F){
				plot(oco2.profiles$pw, oco2.profiles$pres, ylim=c(1000,0), pch=20, xlim=c(0,0.06), main=expression(atop("Pressure weighting functions","at OCO-2 vs. Model levels")),xlab="Pressure Weighting Function (PWF)",ylab="Pressure [hPa]",cex=cex,cex.axis=cex, cex.lab=cex, cex.main=cex)
				lines(oco2.profiles$pw, oco2.profiles$pres,lty=2)
				points(combine.profile$pw[combine.profile$stiltTF], combine.profile$pres[combine.profile$stiltTF], col="red", pch=1,cex=cex)
				points(combine.profile$pw[combine.profile$stiltTF==F], combine.profile$pres[combine.profile$stiltTF==F], col="blue", pch=1,cex=cex)
				abline(h=combine.profile[1,"pres"],lty=3,col="red",cex=cex);abline(h=combine.profile[37,"pres"],lty=3,col="red",cex=cex)
				abline(h=combine.profile[38,"pres"],lty=2,col="blue",cex=cex);abline(h=combine.profile[nrow(combine.profile),"pres"],lty=2,col="blue",cex=cex)
				legend(0.002, 50,c("OCO-2","STILT","CT-NRT"),col=c("black","red","blue"),text.col=c("black","red","blue"),pch=c(20, 1,1),lty=c(2,NA,NA),bty="n",cex=cex-0.3)
			}

			plot(combine.profile$oco2.prior, combine.profile$pres, ylim=c(1000,0), pch=20, xlim=c(385, 405), main=expression(atop("STILT-CT background profiles", "vs. OCO-2 a priori profiles")),xlab="CO2 [ppm]",ylab="Pressure [hPa]",cex=cex,cex.axis=cex, cex.lab=cex, cex.main=cex)
			lines(oco2.profiles$oco2.prior, oco2.profiles$pres,lty=2)
			points(combine.profile$ctnrt.edp[combine.profile$stiltTF==F], combine.profile$pres[combine.profile$stiltTF==F], col="blue", pch=1,cex=cex)
			points(combine.profile$ctnrt.edp[combine.profile$stiltTF], combine.profile$pres[combine.profile$stiltTF], col="red")
			abline(h=combine.profile[1,"pres"],lty=3,col="red",cex=cex);abline(h=combine.profile[37,"pres"],lty=3,col="red",cex=cex)
			abline(h=combine.profile[38,"pres"],lty=2,col="blue",cex=cex);abline(h=combine.profile[nrow(combine.profile),"pres"],lty=2,col="blue",cex=cex)
			legend("bottomleft",c("OCO-2 prior","Background CO2 from traj endpoints","CT CO2 at local receptor"),col=c("black","red","blue"),text.col=c("black","red","blue"),pch=c(20, 1,1),lty=c(2,NA,NA),bty="n",cex=cex-0.3)
			dev.off()
	}

  #-------------------------------------------------- TOTAL ---------------------------------------------------------------------------- #
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
