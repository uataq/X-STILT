#### subroutine to weight the column of footprint for each particle
# create new profile based on TROPOMI CO averaging kernel
#     and then apply to STILT footprint, by weighting the footprints of
#     particles, based on different releasing heights
# written by Dien Wu

wgt.trajec.foot.tropomi = function(output, tropomi.fn, tropomi.species, 
								   ak.wgt = T, pwf.wgt = T, overwriteTF = F) {

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
	# so, use 'foot_before_weight' as initial footprint ('foot') to redo the weighting
	if ( 'foot_before_weight' %in% colnames(trajdat) ) 
		trajdat = trajdat %>% dplyr::select(-starts_with('ak.norm'), 
											-starts_with('pwf'), 
				  							-starts_with('ap'), 
											-starts_with('ak.pwf'), 
											-starts_with('wgt'), 
											-c('foot', 'stiltTF')) %>% 
				  			  dplyr::rename(foot = foot_before_weight)


	# ------------------------------------------------------------------------ #
	# STARTing weighting trajec-level footprint based on vertical profile
	# ------------------------------------------------------------------------ #
	# weighting foot by multipling AK and PW profiles and # of STILT levels/particles
	if ( tropomi.species == 'CO' ) {
		
		cat('\n\nwgt.trajec.foot.tropomi(): weighting using TROPOMI CO profiles...\n')

		# grab press weighting function, pressures, normalized AK
		co.fn  = tropomi.fn[grep('CO', tropomi.fn)]
		out.co = get.wgt.tropomi.func(output, co.fn, 'CO') 

		co.info = as.data.frame(out.co$tropomi.info)
		colnames(co.info)[grepl('tropomi.', colnames(co.info))] = 
			gsub('tropomi', 'tropomi.co', colnames(co.info)[grepl('tropomi.', colnames(co.info))])
		colnames(co.info)[colnames(co.info) == 'xdry.tot'] = 'tropomi.co.xdry.tot'
		co.info$tropomi.co.ak.sfc = out.co$combine.prof$ak.norm[1]

		# combine weighting functions with particles
		xstilt.prof.co = out.co$combine.prof %>% filter(stiltTF == TRUE) %>%
				  		 dplyr::select(indx, ak.norm.co = ak.norm, 
						   			   pwf.co = pwf, ak.pwf.co = ak.pwf)

		# perform weighting for each particle, so merge combine.prof with trajdat 
		trajdat.co = trajdat %>% left_join(xstilt.prof.co, by = 'indx') 
		if (ak.wgt == T & pwf.wgt == T) trajdat.co$wgt.co = trajdat.co$ak.pwf.co * npar
		if (ak.wgt == F & pwf.wgt == T) trajdat.co$wgt.co = trajdat.co$pwf.co * npar 
		if (ak.wgt == T & pwf.wgt == F) trajdat.co$wgt.co = trajdat.co$ak.norm.co 
		if (ak.wgt == F & pwf.wgt == F) trajdat.co$wgt.co = 1
		trajdat = trajdat.co %>% rename(foot_before_wgt = foot) %>% 
								 mutate(foot = foot_before_weight * wgt.no2)

		info = cbind(info, co.info)
	}	# end if TROPOMI CO weighting


	# ------------------------------------------------------------------------ #
	if ( tropomi.species == 'NO2' ) {
		
		cat('\n\nwgt.trajec.foot.tropomi(): weighting using TROPOMI NO2 profiles...\n')

		# grab press weighting function, pressures, normalized AK
		no2.fn  = tropomi.fn[grep('NO2', tropomi.fn)]
		out.no2 = get.wgt.tropomi.func(output, no2.fn, 'NO2')
		
		no2.info = as.data.frame(out.no2$tropomi.info)
		colnames(no2.info)[grepl('tropomi.', colnames(no2.info))] = 
			gsub('tropomi', 'tropomi.no2', colnames(no2.info)[grepl('tropomi.', colnames(no2.info))])
		colnames(no2.info)[colnames(no2.info) == 'xdry.tot'] = 'tropomi.no2.xdry.tot'
		no2.info$tropomi.no2.ak.sfc = out.no2$combine.prof$ak.norm[1]


		# combine NO2 weighting profiles with particles
		xstilt.prof.no2 = out.no2$combine.prof %>% filter(stiltTF == TRUE) %>%
				  		  dplyr::select(indx, ak.norm, pwf, ak.pwf) %>% 
						  rename(ak.norm.no2 = ak.norm, pwf.no2 = pwf, ak.pwf.no2 = ak.pwf)

		# perform weighting for each particle, so merge combine.prof with trajdat 
		trajdat.no2 = trajdat %>% left_join(xstilt.prof.no2, by = 'indx') 
		if (ak.wgt == T & pwf.wgt == T) trajdat.no2$wgt.no2 = trajdat.no2$ak.pwf.no2 * npar
		if (ak.wgt == F & pwf.wgt == T) trajdat.no2$wgt.no2 = trajdat.no2$pwf.no2 * npar 
		if (ak.wgt == T & pwf.wgt == F) trajdat.no2$wgt.no2 = trajdat.no2$ak.norm.no2 
		if (ak.wgt == F & pwf.wgt == F) trajdat.no2$wgt.no2 = 1

		# change former foot to foot_before_wgt, and calculate the weighted foot
		trajdat = trajdat.no2 %>% rename(foot_before_wgt = foot) %>% 
								  mutate(foot = foot_before_weight * wgt.no2)

		info = cbind(info, no2.info)
	} 	# end if TROPOMI NO2 weighting


	# weighting foot by multipling AK and PW profiles and # of STILT levels/particles
	if ( tropomi.species == 'CH4' ) {
		
		cat('\n\nwgt.trajec.foot.tropomi(): weighting using TROPOMI CH4 profiles...\n')

		# grab press weighting function, pressures, normalized AK
		ch4.fn  = tropomi.fn[grep('CH4', tropomi.fn)]
		out.ch4 = get.wgt.tropomi.func(output, ch4.fn, 'CH4') 

		ch4.info = as.data.frame(out.ch4$tropomi.info)
		colnames(ch4.info)[grepl('tropomi.', colnames(ch4.info))] = 
			gsub('tropomi', 'tropomi.ch4', colnames(ch4.info)[grepl('tropomi.', colnames(ch4.info))])
		colnames(ch4.info)[colnames(ch4.info) == 'xdry.tot'] = 'tropomi.ch4.xdry.tot'
		ch4.info$tropomi.ch4.ak.sfc = out.ch4$combine.prof$ak.norm[1]

		# combine weighting functions with particles
		xstilt.prof.ch4 = out.ch4$combine.prof %>% filter(stiltTF == TRUE) %>%
				  		  dplyr::select(indx, ak.norm, pwf, ak.pwf) %>% 
						  rename(ak.norm.ch4 = ak.norm, pwf.ch4 = pwf, ak.pwf.ch4 = ak.pwf)

		# perform weighting for each particle, so merge combine.prof with trajdat 
		trajdat.ch4 = trajdat %>% left_join(xstilt.prof.ch4, by = 'indx') 
		if (ak.wgt == T & pwf.wgt == T) trajdat.ch4$wgt.ch4 = trajdat.ch4$ak.pwf.ch4 * npar
		if (ak.wgt == F & pwf.wgt == T) trajdat.ch4$wgt.ch4 = trajdat.ch4$pwf.ch4 * npar 
		if (ak.wgt == T & pwf.wgt == F) trajdat.ch4$wgt.ch4 = trajdat.ch4$ak.norm.ch4 
		if (ak.wgt == F & pwf.wgt == F) trajdat.ch4$wgt.ch4 = 1
		trajdat = trajdat.ch4 %>% rename(foot_before_wgt = foot) %>% 
								  mutate(foot = foot_before_weight * wgt.ch4)

		info = cbind(info, ch4.info)
	}	# end if TROPOMI CH4 weighting


	# save TROPOMI info in txt file in each by-id folder
	#fn = file.path(dirname(output$file), paste0('tropomi_info_', receptor$long, 
	#											'_', receptor$lati, '.txt'))
	#if (overwriteTF) write.table(info, file = fn, row.names = F, sep = ',', quote = F)

    output$particle = trajdat     # put weighted particle back to output
	saveRDS(output, output$file)

	return(output)
}  # end of subroutine

