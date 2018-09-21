# Main script to obtain the dCO2 for each trajectories and variances of dCO2
# need existing traj, need ak, pw profiles from oco2
# by DW, 01/11/2017

# merge two traject (2 error stat), 03/15/2017
# add biospheric and background portions, 04/11/2017
# use un-weighted trajec for estimating CO2.true (co2.ff, co2.bio and co2.edp)
#   profiles at each level, DW, 04/20/1027
# add more tracks, DW, 11/25/2017

# generalize as a subroutine, DW, 07/19/2018
# path1, path2 for paths that store original trajec vs. trajec with error component
# add vertical integration of errors, DW
# set an upper limit of sd in ppm, e.g., max.sd.trans <- 200

cal.trans.err <- function(namelist, max.sd.trans = 200) {

	# checking...
	check <- function(x, combine.prof) {
		weight <- data.frame(x, level = rep(seq(1, nrow(x)/100), each = 100))
		sum <- sum(tapply(weight$co2.sum, weight$level, mean) *
			combine.prof[combine.prof$stiltTF, 'ak.pwf'])
		return(sum)
	}

	# grab trajec info
	timestr <- namelist$timestr
	pattern <- '_X_wgttraj.rds'
	files1 <- file.path(namelist$traj.path1,
		list.files(namelist$traj.path1, pattern, recursive = T))
	files2 <- file.path(namelist$traj.path2,
		list.files(namelist$traj.path2, pattern, recursive = T))

	outname1 <- gsub(pattern, '', basename(files1))
	recp.info1 <- data.frame(matrix(unlist(strsplit(outname1, '_')), byrow = T,
	  ncol = 3), stringsAsFactors = F)
	colnames(recp.info1) <- c('timestr', 'recp.lon', 'recp.lat')

	recp.lat <- recp.info1$recp.lat
	recp.lon <- recp.info1$recp.lon

  # get emissions
  emiss <- raster(namelist$emiss.file)
	dpar <- namelist$dpar
  xco2.trans <- NULL

	### loop over all unique lat lon and merge two traj
	for (i in 1:length(files1)){

		cat('#--------- ', i / length(recp.lat[i]) * 100, '% ------#\n')

	  # 2. directly read from exiting weighted trajs ----------------------------
		cat('Reading in trajec...\n')
		r1 <- readRDS(files1[i]); p1 <- r1$particle
    r2 <- readRDS(files2[i]);	p2 <- r2$particle
    combine.prof <- r1$wgt.prof

    # combine particle first, correct indx
	  #p2 <-

    # correct trajec format, for STILTv1 output, DW, 07/23/2018
		if (namelist$stilt.ver == 1) {
			p1 <- as.data.frame(p1) %>% dplyr::select(time = time, lati = lat,
				long = lon, zagl = agl, zsfc = zi, foot = foot, indx = index,
				dmas = dmass)
			p2 <- as.data.frame(p2) %>% dplyr::select(time = time, lati = lat,
				long = lon, zagl = agl, zsfc = zi, foot = foot, indx = index,
				dmas = dmass)
		}


		#### 3. USE non-weighted trajec (orig vs. err) ----------------------------
		# to calculate the dCO2.anthro
		# return data frame with index number and corresponding modeled co2
		co2.ff1 <- ff.trajfoot(trajdat = p1, emiss)
		co2.ff2 <- ff.trajfoot(trajdat = p2, emiss)

		# checking...
		xco2.ff1 <- check(x = co2.ff1, combine.prof)
		xco2.ff2 <- check(x = co2.ff2, combine.prof)
		cat(paste('Checking...co2.ff.orig =', signif(xco2.ff1, 4), '[ppm-CO2];',
		  'co2.ff.err =', signif(xco2.ff2, 4), '[ppm-CO2]\n'))


		#### 4. USE non-weighted trajec (orig vs. err) ----------------------------
		# to calculate the dCO2.bio
		co2.bio1 <- bio.trajfoot(trajdat = p1, timestr, namelist$ctflux.path)
	  co2.bio2 <- bio.trajfoot(trajdat = p2, timestr, namelist$ctflux.paht)

		# checking...
		xco2.bio1 <- check(x = co2.bio1, combine.prof)
		xco2.bio2 <- check(x = co2.bio2, combine.prof)
		cat(paste('Checking...bio.orig =', signif(xco2.bio1, 4),
		  '[ppm-CO2]; bio.err =', signif(xco2.bio2, 4), '[ppm-CO2]\n'))


		### 5. USE non-weighted trajec --------------------------------------------
		# to grab background CO2 for STILT levels
		cat('For endpoints Contribution...\n')

		# !!! remember to weight background concentration with Averaging Kernel
		co2.edp1 <- endpts.trajfoot(trajdat = p1, namelist$ctmole.path, combine.prof)
		co2.edp2 <- endpts.trajfoot(trajdat = p2, namelist$ctmole.path, combine.prof)

		### 6. NOW, sum up all 3 contributions ------------------------------------
		# to calculate CO2.true profiles (w/wout errors) for each particle
		merge.co2 <- co2.ff1 %>% full_join(co2.ff2,  by = 'indx') %>%
		  full_join(co2.bio1, by = 'indx') %>% full_join(co2.bio2, by = 'indx') %>%
		  full_join(co2.edp1, by = 'indx') %>% full_join(co2.edp2, by = 'indx') %>%
			mutate(
				tot1 = ff.sum.x + bio.sum.x + edp.sum.x,
				tot2 = ff.sum.y + bio.sum.y + edp.sum.y) %>%
			na.omit()

		# assign level index to co2 for each traj
		merge.co2$level<-	rep(seq(1, nrow(merge.co2)/dpar), each = dpar)

		### 7. calculate the variances at each level ------------------------------
		# v3 for non-removal, using dVAR*AK^2*PW^2
		norm.info <- get.xco2.uncert.v3(logTF = F, rmTF = F, recp.lat = recp.lat[i],
			recp.lon = recp.lon[i], site = site, combine.prof = combine.prof,
			merge.co2 = merge.co2, traj.info = sel.orig.info)

		ln.info <- get.xco2.uncert.v3(logTF = T, rmTF = F, recp.lat = recp.lat[i],
			recp.lon = recp.lon[i], site = site, combine.prof = combine.prof,
			merge.co2 = merge.co2, traj.info = sel.orig.info)

    ### 8. calculate CO2 variances before and after randomizations ------------
		stat1 <- cal.varv2(x = merge.co2, func = 'normal', colname = 'tot1', pct = 0.99)
		stat2 <- cal.varv2(x = merge.co2, func = 'normal', colname = 'tot2', pct = 0.99)
		colnames(stat1) <- c('mean1','sd1','var1','level')
		colnames(stat2) <- c('mean2','sd2','var2','level')

	  # merge statistics and calculate diff in variances, i.e., transort errors
		co2.stat <- full_join(stat.orig, stat.err, by = 'level') %>% na.omit() %>%
		  mutate(dVAR = var2 - var1, dSD = sqrt(dVAR))


		### 9. additional step ----------------------------------------------------
		# to remove negative trans error and output in txt file
		# scale trans errors based on weighted linear regression lines
    lr.stat.info <- scale.dvar(co2.stat = co2.stat)
    co2.stat.wgt <- lr.stat.info[[1]] %>% mutate(sd.trans = sqrt(scaled.dvar),
			sd.trans = replace(sd.trans, is.na(sd.trans), 0))
		#lr <- lr.stat.info[[2]]

		# store vertical profile of trans errors
		if(storeTF) saveRDS(co2.stat.wgt, file = file.path(namelist$store.path,
			paste0(outname1, '_', met, '_info')))


    ### 10. Vertically integrate transport errors for each recp ----------------
    # up to now, we get a transport error per release level
		# we need to vertically integrate errors with correlation length scale
		Lh <- 356  # empirical mean length scale

		# w for total weighting function (AK*PWF)
		prod.trans <- cal.var.cov(sd = sd.trans[sd.trans < max.sd.trans], L = Lh,
			w = combine.prof[combine.prof$stiltTF & sd.trans < max.sd.trans, 'ak.pwf'],
			x = co2.stat.wgt[sd.trans < max.sd.trans,], type = "ver")

		# multiple grid.wgt with prod.trans
		# total SD for trans error for each receptor, aggregated over vertial levels
		tot.trans.sd <- sqrt(sum(prod.trans, na.rm=T))
		print(tot.trans.sd)		# now in ppm
		xco2.trans <- c(xco2.trans, tot.trans.sd)
	} # end of lat loop i

	# write in txtfile
	result <- cbind(timestr, recp.lon, recp.lat, xco2.trans)
  write.table(result, file = namelist$txtfile, sep = ',', quote = F, row.names = F)
}
