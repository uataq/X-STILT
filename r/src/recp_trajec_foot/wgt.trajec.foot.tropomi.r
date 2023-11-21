#### subroutine to weight the column of footprint for each particle
# create new profile based on TROPOMI CO averaging kernel
#     and then apply to STILT footprint, by weighting the footprints of
#     particles, based on different releasing heights
# written by Dien Wu
if (F) {
	tropomi.fn = obs_fn
    tropomi.species = obs_species 
    ak.wgt = ak_wgt
	pwf.wgt = pwf_wgt
}

wgt.trajec.foot.tropomi = function(output, tropomi.fn, 
								   tropomi.species = c('CO', 'CH4', 'NO2')[1], 
								   ak.wgt = T, pwf.wgt = T, overwriteTF = T) {

	# read trajectory before weighting
	trajdat  = output$particle %>% arrange(abs(time))  # now a data.frame
	receptor = output$receptor 
	info = as.data.frame(receptor[c('lati', 'long', 'zsfc', 'psfc')]) %>% 
		   dplyr::rename(mod.zsfc = zsfc, mod.psfc = psfc)

	# before weighting trajec-level footprint by AK & PWF
	# get specific humidity and temp profiles that have been extracted via 
	# before_trajec_xstilt()
	qt.prof = output$qt_prof
	npar = max(trajdat$indx)

	# double check to see if 'foot_before_weight' exists, DW, 07/03/2020
	# if TRUE, it means that footprint has been weighted by AK and PW, 
	# use 'foot_before_weight' as footprint ('foot') to redo weighting
	if ( 'foot_before_weight' %in% colnames(trajdat) ) 
		trajdat = trajdat %>% dplyr::select(-starts_with('ak.norm'), 
											-starts_with('pwf'), 
				  							-starts_with('ap'), 
											-starts_with('ak.pwf'), 
											-starts_with('wgt'), 
											-c('foot')) %>% 
				  			  dplyr::rename(foot = foot_before_weight)


	# ------------------------------------------------------------------------ #
	# STARTing weighting trajec-level footprint based on vertical profile
	# ------------------------------------------------------------------------ #
	# weighting foot with AK and PW profiles and # of STILT levels/particles
	if ( tropomi.species == 'CO' ) {
		
		cat('\n\nwgt.trajec.foot.tropomi(): weighting using TROPOMI CO profiles...\n')

		# grab press weighting function, pressures, normalized AK
		co.fn  = tropomi.fn[grep('CO', tropomi.fn)]
		out.co = get.wgt.tropomi.func(output, co.fn, 'CO') 

		# combine weighting functions with particles
		xstilt.prof.co = out.co$combine.prof %>% filter(stiltTF == TRUE) %>%
				  		 dplyr::select(indx, ak.norm.co = ak.norm, 
						   			   pwf.co = pwf, ak.pwf.co = ak.pwf)

		# perform weighting for each particle 
		trajdat.co = trajdat %>% left_join(xstilt.prof.co, by = 'indx') 
		if (ak.wgt == T & pwf.wgt == T) 
			trajdat.co$wgt.co = trajdat.co$ak.pwf.co * npar
		if (ak.wgt == F & pwf.wgt == T) 
			trajdat.co$wgt.co = trajdat.co$pwf.co * npar 
		if (ak.wgt == T & pwf.wgt == F) 
			trajdat.co$wgt.co = trajdat.co$ak.norm.co 
		if (ak.wgt == F & pwf.wgt == F) trajdat.co$wgt.co = 1
		trajdat = trajdat.co %>% dplyr::rename(foot_before_weight = foot) %>% 
								 mutate(foot = foot_before_weight * wgt.co)
	}	# end if TROPOMI CO weighting


	# ------------------------------------------------------------------------ #
	if ( tropomi.species == 'NO2' ) {
		
		cat('\n\nwgt.trajec.foot.tropomi(): weighting using TROPOMI NO2 profiles...\n')

		# grab press weighting function, pressures, normalized AK
		no2.fn  = tropomi.fn[grep('NO2', tropomi.fn)]
		out.no2 = get.wgt.tropomi.func(output, no2.fn, 'NO2')
		
		# combine NO2 weighting profiles with particles
		xstilt.prof.no2 = out.no2$combine.prof %>% filter(stiltTF == TRUE) %>%
				  		  dplyr::select(indx, ak.norm, pwf, ak.pwf) %>% 
						  rename(ak.norm.no2 = ak.norm, pwf.no2 = pwf, 
						  		 ak.pwf.no2 = ak.pwf)

		# perform weighting for each particle
		trajdat.no2 = trajdat %>% left_join(xstilt.prof.no2, by = 'indx') 
		if (ak.wgt == T & pwf.wgt == T) 
			trajdat.no2$wgt.no2 = trajdat.no2$ak.pwf.no2 * npar
		if (ak.wgt == F & pwf.wgt == T) 
			trajdat.no2$wgt.no2 = trajdat.no2$pwf.no2 * npar 
		if (ak.wgt == T & pwf.wgt == F) 
			trajdat.no2$wgt.no2 = trajdat.no2$ak.norm.no2 
		if (ak.wgt == F & pwf.wgt == F) trajdat.no2$wgt.no2 = 1

		# change initial foot to foot_before_weight, and calculate weighted foot
		trajdat = trajdat.no2 %>% rename(foot_before_weight = foot) %>% 
								  mutate(foot = foot_before_weight * wgt.no2)
	} 	# end if TROPOMI NO2 weighting


	# ------------------------------------------------------------------------ #
	# weighting foot with AK and PW profiles and # of STILT levels/particles
	if ( tropomi.species == 'CH4' ) {
		
		cat('\n\nwgt.trajec.foot.tropomi(): weighting using TROPOMI CH4 profiles...\n')

		# grab press weighting function, pressures, normalized AK
		ch4.fn  = tropomi.fn[grep('CH4', tropomi.fn)]
		out.ch4 = get.wgt.tropomi.func(output, ch4.fn, 'CH4') 

		# combine weighting functions with particles
		xstilt.prof.ch4 = out.ch4$combine.prof %>% filter(stiltTF == TRUE) %>%
				  		  dplyr::select(indx, ak.norm, pwf, ak.pwf) %>% 
						  rename(ak.norm.ch4 = ak.norm, pwf.ch4 = pwf, 
						  		 ak.pwf.ch4 = ak.pwf)
		
		# perform weighting for each particle
		trajdat.ch4 = trajdat %>% left_join(xstilt.prof.ch4, by = 'indx') 
		if (ak.wgt == T & pwf.wgt == T) 
			trajdat.ch4$wgt.ch4 = trajdat.ch4$ak.pwf.ch4 * npar
		if (ak.wgt == F & pwf.wgt == T) 
			trajdat.ch4$wgt.ch4 = trajdat.ch4$pwf.ch4 * npar 
		if (ak.wgt == T & pwf.wgt == F) 
			trajdat.ch4$wgt.ch4 = trajdat.ch4$ak.norm.ch4 
		if (ak.wgt == F & pwf.wgt == F) trajdat.ch4$wgt.ch4 = 1
		trajdat = trajdat.ch4 %>% rename(foot_before_weight = foot) %>% 
								  mutate(foot = foot_before_weight * wgt.ch4)
	}	# end if TROPOMI CH4 weighting

    output$particle = trajdat     # put weighted particle back to output
	if (overwriteTF) saveRDS(output, output$file)

	return(output)
}  # end of subroutine

