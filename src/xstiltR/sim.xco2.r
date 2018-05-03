### Topmost subroutine to simulate column-averaged CO2 dry mole fractions
#   using X-STILT; required subroutines include: get.oco2info(),
#
#----------------------- STEPS involved -------------------------------------- #
# 1) interpolate modeled ground heights based on meteorological fields;

#	2) read in .RData files for trajec, interpolate satellite profiles (averaging
#    kernel AK, pressure weighting PW, and prior profiles),
#    return "combine.prof"; call "weight.trajecfootv2()" to
#    weight trajec-level footprints using interpolated satellite profiles;

# 3) call "get.foot()" and "odiac.anthro()" to generate 2D column footprint
#    and get anthropogenic XCO2 enhancements;

# 4) call "ctnrt.biov2()" to get biospheric XCO2 changes;
# 5) call "ctnrt.oceanv2()" to get oceanic XCO2 changes;
# 6) call "ctnrt.background()" to get XCO2 boundary conditions
#    based on CarbonTracker global fields;
# 7) call "oco2.apriori()" to get prior portion

####------------- Input variables in the namelist include -------------------- #
# sel: XCO2 contributions needed to be modeled, can be subset of
#      c("anthro","bio","ocean","edp","apriori"), where
#      "anthro","bio","ocean" stand for CO2 sources/sinks;
#      "edp" for CO2 boundary condition using CarbonTracker-NearRealTime; and
#      "apriori" for XCO2 contribution from OCO-2 prior profiles.
# site: target city, e.g., "Riyadh";
# timestr: time of satellite overpass, e.g., "";
# ocopath: path that stores OCO-2 Lite files
# recp.lat & recp.lon: receptor lat/lon, same as selected sounding lat/lon;
# recp.info & agl.info: output from ident.to.info(), including receptor time,
#                       lat, lon (recp.info) and release levels (agl.info)
# orig.trajpath & new.trajpath: paths that stores initial trajectories &
#                               weighted trajectories (e.g., by AK or/and PW)

# ak.weight & pw.weight: logical flags for whether weighting against averaging
#                        kernel or pressure weighting, can be T or F
# minlon,maxlon,minlat,maxlat: spatial domains for generating footprints,
# 														 these variables will determine the numpix.x,
#                              numpix.y, lon.ll, lat.ll;

# ct.version: version for CarbonTracker,
#             required when simulating bio and oceanic contributions
# odiac.co2: matrix for anthropogenic CO2 emissions from ODIAC,
#            always readin emissions before call this function
#						 fixed dimension: [Lat, Lon,/Hours], fixed units: micromole/m2/s

# dmassTF: whether to turn on mass conservation, check Trajecfoot() for details
# storeTF: whether to store time-integrated footprint and CO2 contributions
# foot.path, xco2.path: paths for storing time-integrated footprint
#                       and XCO2 contributions fields

#### Written by Dien Wu, 10/20/2017

#------------------------- Updates ------------------------------------------ #
# use ODIACv2016 and CT-NRTv2016-1 as default input dataset
# ADD biospheric CO2 from CT-NRT, 02/10/2017
# ADD ODIACv2016, 02/15/2017
# ADD hourly ODIAC emissions, 03/09/2017
# Get ak, pw profiles interpolation code separated by weight trajecfoot code
# Add grabbing upper background profiles from CT-NRT
# ADD oceanic fluxes from CT-NRT, 07/14/2017
# update two subroutines for bio and ocean, 09/15/2017
# add ak.weight --> no need to weight through AK and apriori=0, DW, 02/06/2018
# use namelist instead of individual variables, DW, 05/02/2018
# combine oceanic and biospheric fluxes, DW, 05/03/2018
###############################################################################

sim.xco2 <- function(namelist, recp.info, odiac.co2){

  #--------------------------------------------------------------------------- #
  # 1. get OCO-2 info given recp.info
	# requires get.oco2info() from "OCO2.get.oco2info.r"
  YYYYMMDD <- substr(namelist$timestr,1,8)
  oco2.file <- list.files(pattern=YYYYMMDD, path=namelist$ocopath)

	cat(paste("sim.xco2(): grabbing satellite profiles given recp lat/lon...\n"))
  oco2.info <- get.oco2info(ocopath=namelist$ocopath, ocofile=oco2.file,
		                        recp.lat=recp.info$recp.lat,
														recp.lon=recp.info$recp.lon)

  wgt.profile <- oco2.info[[1]]    # lists of vertical profiles for each recp
  oco2.dat <- oco2.info[[2]]      # include obs xco2, other info
  id <- as.numeric(as.character(oco2.dat$find.id))	# OCO-2 ID

  # get LOWER LEFT and lat/lon domain for footprint
	lon.ll <- namelist$minlon
	lat.ll <- namelist$minlat
  dlat <- namelist$maxlat - lat.ll
  dlon <- namelist$maxlon - lon.ll

  # write headers in txt file
  filenm <- file.path(namelist$txtpath, namelist$txtname)

  # headers for txtfile
  headers <- c("SoundingID", "recp.lat", "recp.lon", "sim.grdhgt",
                paste("xco2.", namelist$sel, sep=""),
               "obs.xco2", "sim.xco2", "diff.xco2")
  ncol <- length(headers)
  write(headers, file=filenm, sep=",", ncol=ncol)

  # loop over each receptor...
  for (i in 1:nrow(recp.info)){

    cat(paste("#------", signif(i/nrow(recp.info)*100,3), "% Done ------- #\n"))

		## get receptor info, in vector form
		tmp.lat <- recp.info$recp.lat[i]  # receptor lat
		tmp.lon <- recp.info$recp.lon[i]  # receptor lon
		tmp.ident <- recp.info$recp.ident[i]   # Rdata filename of STILT trajec
		tmp.grdhgt <- recp.info$recp.grdhgt[i]   # get modeled ground heights
    tmp.recp.info <- recp.info[i,]
    tmp.agl.info <- agl.info[i,]
    tmp.wgt.prof <- wgt.profile[[i]]  # get weighting profiles

		## initialize result, store sounding ID, recptor lat/lon, ground hgt
	  tmp.result <- c(as.character(id[i]), tmp.lat, tmp.lon, tmp.grdhgt)
		total.sim.xco2<-0   		# initialize total modeled XCO2

		#------------------------------------------------------------------------- #
		# 2. read in traj and apply ak and pw profiles onto original trajs
		cat("sim.xco2(): Reading unweighted trajec ...\n")
		tmp.ident <- substr(tmp.ident, 1, nchar(tmp.ident)-6)
		orig.trajdat <- getr(tmp.ident, path=namelist$orig.trajpath)

		### start weighting trajec-level footprint...call weight.trajecfootv2()
		# new.trajpath: path for storing weighted trajectory
		cat("sim.xco2(): Weighting trajec-level footprint...\n")
		scale.prof <- weight.trajecfootv2(ident = tmp.ident, trajdat = orig.trajdat,
																		  recp.info = tmp.recp.info,
																			agl.info = tmp.agl.info,
																		  oco.ak.norm = tmp.wgt.prof$ak.norm,
																		  oco.pw = tmp.wgt.prof$pw,
																		  oco.pres = tmp.wgt.prof$pres,
																		  oco.apriori = tmp.wgt.prof$apriori,
																		  recp.grdhgt = tmp.grdhgt,
																		  new.trajpath = namelist$wgt.trajpath,
																		  ak.weight = namelist$ak.weight,
																			pw.weight = namelist$pw.weight)

		## returns the ak pw profiles at for all levels & the weighted footprint
		combine.prof <- scale.prof[[1]]	# all profile info, ak, pw, apriori
		wgt.trajdat  <- scale.prof[[2]]	# weighted traj still in matrix form

		#------------------------ 3. ANTHROPOGENIC  ------------------------------ #
		if(!is.na(match("anthro", namelist$sel))){

			cat("\nsim.xco2(): Dealing with enhancements due to FFCO2 emissions...\n")

			# 3.1. GENERATE 2D footprint based on the above weighted traj,
			# call get.foot() from "OCO2.get.foot.r"
			lat.res <- namelist$lat.res.anthro; numpix.y <- dlat/lat.res
			lon.res <- namelist$lon.res.anthro; numpix.x <- dlon/lon.res

			# compute footprint time stamps, if only 2 elements in "foottimes",
			# we need to return the time-integrated 2D footprint (xfoot),
      foottimes <- c(0,abs(namelist$max.hour))  # monthly ODIAC

			# if >2 elements in "foottimes", return time-varying 3D footprint (foot)
			# when coupling with hourly ODIAC
			if(namelist$hourlyTF)foottimes <- seq(0,abs(namelist$max.hour),1)
			foot.anthro <-get.foot(ident=tmp.ident,foot.overwrite=T, lat.res=lat.res,
				                     trajdat=wgt.trajdat,trajpath=namelist$wgt.trajpath,
				                     footpath=namelist$foot.path, foottimes=foottimes,
									           dmassTF=namelist$dmassTF, numpix.x=numpix.x,
														 storeTF=namelist$storeTF, numpix.y=numpix.y,
														 lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res)

			# 3.2. CALL the subroutine for matching ODIAC emission with xfoot,
			#      and sum the matrix to obtain "dxco2.anthro"
			cat("sim.xco2(): Calculating the dCO2 from anthro using ODIAC...\n")
			xco2.anthro <- odiac.anthro(foot=foot.anthro, odiac.co2=odiac.co2,
				                          ident=tmp.ident, storeTF=namelist$storeTF,
																	ncdfpath=namelist$xco2.path)
			cat(paste("sim.xco2(): XCO2.anthro: ", signif(xco2.anthro, 4), "ppm\n"))

			# store result and sum up total XCO2
			total.sim.xco2 <- total.sim.xco2 + xco2.anthro
			tmp.result <- c(tmp.result, xco2.anthro)
		}	# end "anthro" contribution


		#-------------------- 4. BIOSPERHIC + OCEANIC ---------------------------- #
		if(!is.na(match("bio", namelist$sel))|!is.na(match("ocean", namelist$sel))){

			cat("\nsim.xco2(): Dealing with XCO2 changes due to CT-NRT bio/ocean...\n")

			# 4.1. GENERATE 2D footprint that matches CT
			lat.res <- namelist$lat.res.ct; numpix.y <- dlat/lat.res
			lon.res <- namelist$lon.res.ct; numpix.x <- dlon/lon.res
      foottimes <- seq(0, abs(namelist$max.hour), 1)
			foot.ct <- get.foot(ident=tmp.ident, foot.overwrite=T,
				                   trajdat=wgt.trajdat, trajpath=namelist$wgt.trajpath,
													 foottimes=foottimes, footpath=namelist$foot.path,
									         dmassTF=namelist$dmassTF, storeTF=namelist$storeTF,
													 numpix.x=numpix.x, numpix.y=numpix.y, lon.ll=lon.ll,
                           lat.ll=lat.ll, lon.res=lon.res, lat.res=lat.res)

			# 4.2. CALL the subroutine for matching biosperhic fluxes with footprint,
			#      and sum the matrix to obtain "xco2.bio" or "xco2.ocean"
			cat("sim.xco2(): Calculating the dCO2 from bio using CT-NRT...\n")
			xco2.ct <- ctnrt.bio.ocean(ident=tmp.ident, recp.info=tmp.recp.info,
				                         ct.version=namelist$ct.version,
															   ctpath=namelist$ctflux.path,
				                         storeTF=namelist$storeTF, foot=foot.ct,
															   ncdfpath=namelist$xco2.path)

      xco2.bio <- xco2.ct[1]; xco2.ocean <- xco2.ct[2]
			cat(paste("sim.xco2(): XCO2.bio: ", signif(xco2.bio, 4), "ppm\n"))
      cat(paste("sim.xco2(): XCO2.ocean: ", signif(xco2.ocean, 4), "ppm\n\n"))
			# if foot hours differ from CT hours, the subroutine will return NA values

			# 4.3. store result and sum up total XCO2
			total.sim.xco2 <- total.sim.xco2 + xco2.bio + xco2.ocean
			tmp.result <- c(tmp.result, xco2.bio, xco2.ocean)
		}	# end "bio" and "ocean" sections


		#----------------- 5. ENDPOINTS for boundary conditions ------------------ #
		if(!is.na(match("edp", namelist$sel))){

			# 5.1. CALL the subroutine for calculating background CO2
      # (AK PW weighted CO2 concentration derived from CT-NRT at endpoints)
			# sum all to get "XCO2.edp", remember to use weighted traj "wgt.trajdat"
			edp.info <- ctnrt.traj.edp(ident=tmp.ident, trajdat=wgt.trajdat,
                                 recp.info=tmp.recp.info,
                                 combine.profile=combine.prof,
                                 ct.version=namelist$ct.version,
															   ctpath=namelist$ctmole.path)

    	xco2.edp<-edp.info[[2]]
			combine.prof<-edp.info[[1]]
			cat(paste("XCO2 boundary condition: ", signif(xco2.edp,5), "ppm\n\n"))

			# 6.2. store result and sum up total XCO2
			total.sim.xco2 <- total.sim.xco2 + xco2.edp
			tmp.result <- c(tmp.result, xco2.edp)
		}	# end "edp" sel

		#-------------- 6. APRIORI contributions from OCO-2 ---------------------- #
		if(!is.na(match("apriori", namelist$sel))){

			# 6.1 CALL the subroutine for calculating the contribution of CO2 from
      # a priori profile (I - AK), sum all to obtain "XCO2.apriori"
			xco2.prior <- oco2.apriori(combine.profile=combine.prof)
			cat(paste("XCO2 prior portion:", signif(xco2.prior,5), "ppm\n\n"))

			# 6.2. store result and sum up total XCO2
			total.sim.xco2 <- total.sim.xco2 + xco2.prior
			tmp.result <- c(tmp.result, xco2.prior)
		}	# end "apriori" sel

		######
		#---------------- 7. calculate the TOTAL modeled XCO2 -------------------- #
		if(length(namelist$sel) >= 3){

			# 7.1. GRAB observed XCO2 for each sounding
			tmp.obs.xco2 <- oco2.dat$xco2[i]
			cat(paste("TOTAL retrieved XCO2:", signif(tmp.obs.xco2, 5), "ppm\n"))

			# 7.2. SUM all contributions to get the final total simulated
      # retrieved XCO2.ak, which can be compared with retrieved XCO2 from OCO-2
			cat(paste("TOTAL modeled XCO2:", signif(total.sim.xco2, 5), "ppm\n\n"))
			diff.xco2 <- total.sim.xco2 - tmp.obs.xco2	           # sim - obs
			# end calculating the total and model-data mismatch

			# 7.3. add results in a txt file,
      # i.e., OCO-2 observed XCO2, diff between observed and simulated XCO2
			tmp.result <- c(tmp.result, tmp.obs.xco2, total.sim.xco2, diff.xco2)
		}

    # write in txtfile...
    write(tmp.result, file=filenm, append=T, sep=",", ncol=ncol)
	}

  cat(paste("Model result in", filenm,"\n\n"))
  cat("End of X-STILT XCO2 simulations!\n")
	#return(result)
}


# end of subroutine



### old code
if(F){

	#--------------------- 5. OCEANIC ---------------------------------------- #
	if(!is.na(match("ocean", namelist$sel))){

		# 5.1. GENERATE 2D footprint based on the above newly weighted traj
		# no need to store 1 deg footprint

		# if we already have biospheric footprint (same as oceanic footprint)
		if(!is.na(match("bio", namelist$sel) * match("ocean", namelist$sel))){
			foot.ocean <- foot.bio
		}else{
      lat.res <- namelist$lat.res.ocean; numpix.y <- dlat/lat.res
			lon.res <- namelist$lon.res.ocean; numpix.x <- dlon/lon.res
      foottimes <- seq(0, abs(namelist$max.hour), 1)

      cat("sim.xco2(): Generating footprint that matches CT-NRT ocean...\n")
      foot.ocean <-get.foot(ident=tmp.ident,foot.overwrite=T, lat.res=lat.res,
				                    trajdat=wgt.trajdat, storeTF=namelist$storeTF,
				                    footpath=namelist$foot.path, foottimes=foottimes,
									          dmassTF=namelist$dmassTF, numpix.x=numpix.x,
                            trajpath=namelist$wgt.trajpath,numpix.y=numpix.y,
													  lon.ll=lon.ll, lat.ll=lat.ll, lon.res=lon.res)
		}

		# 5.2. CALL the subroutine for matching biosperhic fluxes with footprint,
    #      and sum the matrix to obtain "xco2.bio"
		cat("sim.xco2(): Calculating the dCO2 from CT-NRT oceanic fluxes...\n")
		xco2.ocean <- ctnrt.oceanv2(ident=tmp.ident, foot=foot.ocean,
                                ct.version=namelist$ct.version,
                                ctpath=namelist$ctflux.path,
                                storeTF=namelist$storeTF, res=lon.res,
                                ncdfpath=namelist$xco2.path)
		cat(paste("XCO2 from oceanic fluxes: ", signif(xco2.ocean, 4), "ppm\n\n"))
		# if foot hours differ from CT hours, the subroutine will return NA values

		# 5.3. store result and sum up total XCO2
		total.sim.xco2<-total.sim.xco2+dxco2.ocean
		tmp.result<-c(tmp.result,dxco2.ocean)

	}	# end "ocean" sel

}
