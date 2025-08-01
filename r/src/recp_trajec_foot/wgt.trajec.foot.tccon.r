

# perform trajec-level footprint weighting using TCCON weighting profiles
# DW, 04/21/2023
wgt.trajec.foot.tccon = function(output, tccon.fn = NA, 
								 tccon.species = c('CH4', 'CO', 'CO2', 
								 				   'H2O', 'N2O')[3]){

	# read trajectory before weighting
	p = output$particle %>% arrange(abs(time)) 	# now a data.frame

	# ------------------------------------------------------------------------ #
	# double check to see if 'foot_before_weight' exists, DW, 07/03/2020
	# if TRUE, it means that footprint has been weighted by AK and PW, 
	# use 'foot_before_weight' as initial footprint ('foot') to redo weighting
	if ( 'foot_before_weight' %in% colnames(p) ) 
		p = p %>% dplyr::select(-starts_with('ak.norm'), -starts_with('ap'), 
							 -c('foot', 'stiltTF', 'pwf', 'ak.pwf', 'wgt')) %>%
				  dplyr::rename(foot = foot_before_weight)

	# add PW and Ak weighting for trajec with error as well, DW, 10/21/2018
	errTF = 'particle_error' %in% names(output)
	if (errTF) {
		p.err = output$particle_error %>% arrange(abs(time)) 
		if ( 'foot_before_weight' %in% colnames(p.err) ) 
			p.err = p.err %>% 
					dplyr::select(-starts_with('ak.norm'), -starts_with('ap'), 
								  -starts_with('wgt'), 
								  -c('foot', 'stiltTF', 'pwf', 'ak.pwf')) %>% 
				  	rename(foot = foot_before_weight)
	}	# end if errTF


	# ------------------------------------------------------------------------ #
	# use retrieved surface pressure and height for pressure weighting
	# use retrieved Ak for AK weighting
    #source('r/dependencies.r')

	combine.prof = get.wgt.tccon.func(output, tccon.fn, tccon.species)
	colnames(combine.prof)[ colnames(combine.prof) == 'ap_gas_wet'] = 
		paste0('ap_', tolower(tccon.species), '_wet')

	output$combine.prof = combine.prof
	
	xstilt.prof = combine.prof %>% filter(stiltTF == TRUE) %>% 
					dplyr::select(-c('pres'))

	### STARTing weighting trajec based on profile
	# weighting foot by multipling AK and PW profiles from 'xstiltprof'
	# need to multiple PWF for model levels with number of particles 
	# since calc_footprint() calculates the average spatial foot based on trajec-level foot
	
	# perform weighting for each particle or level
	npar = max(p$indx)
	p.wgt = p %>% left_join(xstilt.prof, by = 'indx') %>% 
			mutate(wgt = ak.pwf * sf.wet * npar, foot_wgt = foot * wgt) %>%
			rename(foot_before_weight = foot, foot = foot_wgt)

	if (F) {	# SANITY CHECK
		w1 = p.wgt %>% filter(time == -1)
		p1 = ggplot(data = w1) + 
			 geom_point(aes(ak.norm, xhgt, color = 'ak of xco2'), size = 0.1) + 
			 geom_point(aes(sf.wet, xhgt, color = 'sf for wetting foot'), size = 0.1) + 
			 geom_point(aes(pwf * npar, xhgt, color = 'interpolated int operator'), size = 0.1)
	}

	colnames(p.wgt)
	output$particle = p.wgt 

	if (errTF) {
		p.err.wgt = p.err %>% left_join(xstilt.prof, by = 'indx') %>%
					mutate(wgt = ak.pwf * sf.wet * npar, 
						   foot_wgt = foot * wgt) %>% 
					rename(foot_before_weight = foot, foot = foot_wgt)
		output$particle_error = p.err.wgt
	}


	# return both weighting profiles and weighted trajec
	#saveRDS(output, output$file) 	# overwrite the "X_traj.rds" file
	return(output)

}  # end of subroutine

