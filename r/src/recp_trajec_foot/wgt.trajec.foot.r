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

# oco.fn can be NA for ideal simulation by setting ak.wgt as FALSE (AK = 1), DW, 12/06/2020

wgt.trajec.foot = function(output, oco.fn = NA, ak.wgt = T, pwf.wgt = T){

	if (ak.wgt == F & pwf.wgt == F) {
		# if ak.wgt == F && pwf.wgt == F, return trajec with original footprint,
		# no longer need any following weighting...
		# !!! but still need to return weighting functions and other info
		cat('wgt.trajec.foot.oco (): NO need to perform ANY vertical weighting...\n')
		return(output)
	} 

	# read trajectory before weighting
	trajdat = output$particle %>% arrange(abs(time)) 	# now a data.frame

	# ------------------------------------------------------------------------ #
	# double check to see if 'foot_before_weight' exists, DW, 07/03/2020
	# if TRUE, it means that footprint has been weighted by AK and PW, 
	# so, use 'foot_before_weight' as initial footprint ('foot') to redo the weighting
	if ( 'foot_before_weight' %in% colnames(trajdat) ) 
		trajdat = trajdat %>% dplyr::select(-starts_with('ak.norm'), 
											-starts_with('pwf'), 
				  							-starts_with('ap'), 
											-starts_with('ak.pwf'), 
											-starts_with('wgt'), 
											-c('foot', 'stiltTF')) %>% 
				  			  dplyr::rename(foot = foot_before_weight)

	# add PW and Ak weighting for trajec with error as well, DW, 10/21/2018
	errTF = 'particle_error' %in% names(output)
	if (errTF) {
		trajdat.err = output$particle_error %>% arrange(abs(time)) 

		if ( 'foot_before_weight' %in% colnames(trajdat.err) ) 
			trajdat.err = trajdat.err %>% dplyr::select(-starts_with('ak.norm'), 
														-starts_with('pwf'), 
														-starts_with('ap'), 
														-starts_with('ak.pwf'), 
														-starts_with('wgt'), 
														-c('foot', 'stiltTF')) %>% 
				  			  			  dplyr::rename(foot = foot_before_weight)

			#trajdat.err[, c('foot', 'ak.norm', 'pwf', 'ap', 'ak.pwf', 'stiltTF', 'wgt')] = NULL
			#trajdat.err = trajdat.err %>% dplyr::rename(foot = foot_before_weight)
	}	# end if errTF


	# ------------------------------------------------------------------------ #
	# if ak.wgt == FALSE, remove dependence of any satellite, DW, 09/15/2020 
	# use modeled surface pressure and height for pressure weighting
	if (!ak.wgt) xstilt.prof = get.wgt.mod.func(output)
	
	# if ak.wgt == TRUE, keep OCO dependence, DW, 09/15/2020 
	# use retrieved surface pressure and height for pressure weighting
	# use retrieved Ak for AK weighting
	if (ak.wgt) xstilt.prof = get.wgt.oco.func(output = output, oco.fn = oco.fn) %>% 
				    		  dplyr::select(c('indx', 'ak.norm', 'pwf', 'ap', 
							   				  'ak.pwf', 'stiltTF')) %>% 
							  filter(stiltTF == TRUE)


	### STARTing weighting trajec based on profile
	# weighting foot by multipling AK and PW profiles from 'xstiltprof'
	# need to multiple PWF for model levels with number of particles 
	# since calc_footprint() calculates the average spatial foot based on trajec-level foot
	
	# perform weighting for each particle or level, so merge xstilt.prof with trajdat 
	npar = max(trajdat$indx)
	trajdat.wgt = trajdat %>% left_join(xstilt.prof, by = 'indx') 
	if (ak.wgt == T & pwf.wgt == T) trajdat.wgt$wgt = trajdat.wgt$ak.pwf * npar
	if (ak.wgt == F & pwf.wgt == T) trajdat.wgt$wgt = trajdat.wgt$pwf * npar 
	if (ak.wgt == T & pwf.wgt == F) trajdat.wgt$wgt = trajdat.wgt$ak.norm 
	trajdat.wgt = trajdat.wgt %>% mutate(foot_wgt = foot * wgt) %>% 
			      rename(foot_before_weight = foot, foot = foot_wgt)
	output$particle = trajdat.wgt 


	if (errTF) {
		trajdat.err.wgt = trajdat.err %>% left_join(xstilt.prof, by = 'indx') 
		if (ak.wgt == T & pwf.wgt == T) trajdat.err.wgt$wgt = trajdat.err.wgt$ak.pwf * npar
		if (ak.wgt == F & pwf.wgt == T) trajdat.err.wgt$wgt = trajdat.err.wgt$pwf * npar 
		if (ak.wgt == T & pwf.wgt == F) trajdat.err.wgt$wgt = trajdat.err.wgt$ak.norm 
		trajdat.err.wgt = trajdat.err.wgt %>% mutate(foot_wgt = foot * wgt) %>% 
						  rename(foot_before_weight = foot, foot = foot_wgt)
		output$particle_error = trajdat.err.wgt
	}


	# return both weighting profiles and weighted trajec
	saveRDS(output, output$file) 	# overwrite the "X_traj.rds" file
	return(output)

}  # end of subroutine




	
# check the sensitivity of having AK and diff in surface pres 
if (F) {

	x1 = get.wgt.mod.func(output)
	x2 = get.wgt.oco.func(output, oco.fn) %>% filter(stiltTF == TRUE)

	xx = ggplot() + xlim(c(0.15, 1.1)) + theme_bw() + xlab('weighting functions') + 
		geom_point(data = x1, aes(ak.norm, indx, col = 'AK_norm if ak.wgt == F'), size = 0.3) + 
		geom_point(data = x2, aes(ak.norm, indx, col = 'AK_norm if ak.wgt == T'), size = 0.3) + 
		geom_point(data = x1, aes(pwf * 3000, indx, col = 'PWF if ak.wgt == F (using modeled Psfc)'), size = 0.3) + 
		geom_point(data = x2, aes(pwf * 3000, indx, col = 'PWF if ak.wgt == T (using retrieved Psfc)'), size = 0.3) +
		scale_color_manual(values = c('red', 'orange', 'purple', 'blue'))

}



