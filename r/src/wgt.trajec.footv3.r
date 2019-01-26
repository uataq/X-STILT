#### subroutine to weight the column of footprint for each particle
# create new profile based on OCO-2 pressure weighting and averaging kernel
#     and then apply to STILT footprint, by weighting the footprints of
#     particles, based on different releasing heights
# OCO-2 only provides the PW, AK and a priori at 20 levels, use linear
#     interpolation to 'approx' values at given STILT releasing levels
# written by Dien Wu

# updates--
# add 'get.weight.func()' for obtaining interpolated weighting function
# trajec with level column is not passed to weight.trajecfoot(), DW, 05/22/2017
# version 3 for matching Ben's STILTv2, DW, 05/29/2018
# use 'xhgt' instead of 'level', DW, 06/01/2018
# add PW and Ak weighting for trajec with error and fix a bug, DW, 10/21/2018

wgt.trajec.footv3 <- function(output, oco2.info, ak.wgt = T, pwf.wgt = T){

  # read trajectory before weighting
  trajdat <- output$particle  # now a data.frame
  trajdat <- trajdat[order(abs(trajdat$time)), ]    # order by time

	# add PW and Ak weighting for trajec with error as well, DW, 10/21/2018
	errTF <- 'particle_error' %in% names(output)
	if (errTF) {
		trajdat.err <- output$particle_error 
		trajdat.err <- trajdat.err[order(abs(trajdat.err$time)), ]
	}

	# HERE, ak.wgt and pwf.wgt is passed on for weighting trajec
	combine.prof <- get.wgt.funcv3(output = output, oco2.info = oco2.info,
		                             ak.wgt = ak.wgt, pwf.wgt = pwf.wgt)

	### STARTing weighting trajec based on profiles
	if (ak.wgt == F & pwf.wgt == F) {

		# if ak.wgt == F && pwf.wgt == F, return trajec with original footprint,
		# no longer need any following weighting...
		# !!! but still need to return weighting functions and other info
		cat('weight.trajecfootv3(): NO weighting turned on...\n')
		wgt.prof <- trajdat  # still the same trajec
    if (errTF) wgt.prof <- list(trajdat, trajdat.err)

	} else {

		#### --------- START WEIGHTING FOOTPRINT COLUMN FROM .RData FILE --------- #
		# group particles, sort traj files by 'indx', 05/22/2017
		# add one more column for release level to which particles belong
    # since 'xhgt' has been stored using Ben's code, indicating initial release
    # hgts, use 'xhgt' instead of 'level', DW, 06/01/2018
    uni.xhgt <- unique(trajdat$xhgt)
    nlevel <- length(uni.xhgt)

		# initialize weighted foot column with normal footprint
		trajdat$newfoot <- NA; if (errTF) trajdat.err$newfoot <- NA

		# weighting newfoot by multipling AK and PW profiles from 'combine.prof',
		# along with number of STILT levels
		stilt.prof <- combine.prof[combine.prof$stiltTF == TRUE, ]

		# DW, 04/20/2017, add pwf.wgt flag too
		# only weight footprint in trajec if one of the two flags/or both are TRUE
		if (ak.wgt == T & pwf.wgt == T) {
		  cat('weight trajec by both AK & PW profiles...\n')
			wgt.prof <- stilt.prof$ak.pwf

		} else if (ak.wgt == F & pwf.wgt == T) {
      cat('weight trajec only by PW profiles\n')
			wgt.prof <- stilt.prof$pwf

		} else if (ak.wgt == T & pwf.wgt == F) {
			cat('weight trajec only by AK profiles\n')
			wgt.prof <- stilt.prof$ak.norm
		}

		# start weighting for unique release levels
		for (h in 1:length(uni.xhgt)) {
			hgt.indx <- which(trajdat$xhgt == uni.xhgt[h])

      # need to multiple by number of levels, as trajecfoot()/calc_footprint()
      # calculates the spatial footprint based on average footprint in a column
      # thus, resultant 'newfoot' should have similar order of mag as of 'foot'
			trajdat$newfoot[hgt.indx] <- trajdat$foot[hgt.indx] * wgt.prof[h] * nlevel

			# also weight the trajec with error, DW, 10/21/2018
			if (errTF) 
			  trajdat.err$newfoot[hgt.indx] <- 
				  trajdat.err$newfoot[hgt.indx] * wgt.prof[h] * nlevel
		} # end loop h

	} # end if all flags, ak.wgt & pwf.wgt

	# for testing, store two sets of trajdat
	# one weighting over AK.norm * PW, newfoot are much smaller than original foot
	newtraj <- trajdat[, -which(colnames(trajdat) == 'foot')]
	colnames(newtraj)[colnames(newtraj) == 'newfoot'] <- 'foot'

  if (errTF) {
		newtraj.err <- trajdat.err[, -which(colnames(trajdat.err) == 'foot')]
		colnames(newtraj.err)[colnames(newtraj.err) == 'newfoot'] <- 'foot'
	}

	# add interpolated profiles in RData files as well, DW, 04/19/2017
	# put 'newtraj' back to 'output'
	wgt.output <- output
	wgt.output$particle <- newtraj # overwrite with weighted trajec
	if (errTF) wgt.output$particle.err <- newtraj.err

  # add interpolated AK, PW profiles
	wgt.output$wgt.prof <- combine.prof  
	wgt.output$file <- gsub('X_traj.rds', 'X_wgttraj.rds', output$file)
  saveRDS(wgt.output, wgt.output$file)

	# return both weighting profiles and weighted trajec
	return(wgt.output)

}  # end of subroutine
