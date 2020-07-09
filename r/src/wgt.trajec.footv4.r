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
# add PW and Ak weighting for trajec with error and fix a bug, DW, 10/21/2018
# minor bug and typo fixed for weighting trajec with error, DW, 02/21/2019 
# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# stop generating a new "wgttraj.rds", merge two columns of footprint in one file
#    to reduce the storage, DW, 07/03/2020 
# interpolate vertical profiles onto each particle instead of each release level, 
#    due to changes in HYSPLIT compiler, DW, 07/07/2020

wgt.trajec.footv4 <- function(output, oco.info, ak.wgt = T, pwf.wgt = T){

	if (ak.wgt == F & pwf.wgt == F) {

		# if ak.wgt == F && pwf.wgt == F, return trajec with original footprint,
		# no longer need any following weighting...
		# !!! but still need to return weighting functions and other info
		cat('wgt.trajec.footv3(): NO need to perform vertical weighting...\n')
		return(output)
	} 


	# read trajectory before weighting
	trajdat <- output$particle  # now a data.frame
	trajdat <- trajdat[order(abs(trajdat$time)), ]    # order by time

	# double check to see if 'foot_before_weight' exists, DW, 07/03/2020
	# if TRUE, it means that footprint has been weighted by AK and PW, 
	# so, use 'foot_before_weight' as initial footprint ('foot') to redo the weighting
	if ( 'foot_before_weight' %in% colnames(trajdat) ) {
		trajdat$foot <- NULL 
		trajdat <- trajdat %>% dplyr::rename(foot = foot_before_weight)
	}
	
	# add PW and Ak weighting for trajec with error as well, DW, 10/21/2018
	errTF <- 'particle_error' %in% names(output)
	if (errTF) {
		trajdat.err <- output$particle_error 
		trajdat.err <- trajdat.err[order(abs(trajdat.err$time)), ]

		if ( 'foot_before_weight' %in% colnames(trajdat.err) ) {
			trajdat.err$foot <- NULL 
			trajdat.err <- trajdat.err %>% dplyr::rename(foot = foot_before_weight)
		}
	}	# end if errTF

	# ------------------------------------------------------------------------ #
	# HERE, ak.wgt and pwf.wgt is passed on for weighting trajec
	combine.prof <- get.wgt.funcv6(output, oco.info, ak.wgt, pwf.wgt)

	### STARTing weighting trajec based on profile
	# weighting newfoot by multipling AK and PW profiles from 'combine.prof',
	# along with number of STILT levels
	stilt.prof <- combine.prof[combine.prof$stiltTF == TRUE, ] %>% 
				  dplyr::select(indx, ak.norm, ap, pwf, ak.pwf)

	# perform weighting for each particle, so merge combine.prof with trajdat 
	trajdat <- trajdat %>% left_join(stilt.prof, by = 'indx') 
	max.indx <- max(trajdat$indx)

	# DW, 04/20/2017, add pwf.wgt flag too
	# only weight footprint in trajec if one of the two flags/or both are TRUE
	if (ak.wgt == T & pwf.wgt == T) {
		cat('weight trajec by both AK & PW profiles...\n')
		trajdat$wgt <- trajdat$ak.pwf

	} else if (ak.wgt == F & pwf.wgt == T) {
		cat('weight trajec only by PW profiles\n')
		trajdat$wgt <- trajdat$pwf

	} else if (ak.wgt == T & pwf.wgt == F) {
		cat('weight trajec only by AK profiles\n')
		trajdat$wgt <- trajdat$ak.norm
	}

	trajdat <- trajdat %>% mutate(foot.wgt = foot * wgt * max.indx) %>% 
						   dplyr::select(-c(ak.norm, ap, ak.pwf, wgt, pwf))
	if (errTF) trajdat.err <- trajdat.err %>% mutate(foot.wgt = foot * wgt * max.indx) %>%
							  dplyr::select(-c(ak.norm, ap, ak.pwf, wgt, pwf))

	# stop generating a new "wgttraj.rds", merge two columns of footprint in one file
	#    to reduce the storage, DW, 07/03/2020 
	# and put 'newtraj' back to 'output'
	newtraj <- trajdat %>% rename(foot_before_weight = foot, foot = foot.wgt)
	output$particle <- newtraj			 # overwrite with weighted trajec

  	if (errTF) {
		newtraj.err <- trajdat.err %>% rename(foot_before_weight = foot, foot = foot.wgt)
		output$particle_error <- newtraj.err
	}

	# add interpolated AK and PW profiles, DW, 04/19/2017
	saveRDS(output, output$file) 	# overwrite the "X_traj.rds" file

	# return both weighting profiles and weighted trajec
	return(output)
}  # end of subroutine
