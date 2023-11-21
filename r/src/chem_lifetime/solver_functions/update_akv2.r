
# output - particle
# no2_fn - updated NO2 file
# no2_info_v1 - no2 profile from old calcs
if (F) {
    receptor = tmp_receptor
    p_chem = tmp_chem
    no2_fn = pal_fn 
    no2_info_v1 = tmp_info
}

update_akv2 = function(receptor, p_chem, no2_fn, no2_info_v1) {    
    
    air_vcd_stilt = no2_info_v1$air_vcd_stilt
    air_vcd_tropo = no2_info_v1$air_vcd_tropo

    # get NO2 profile from updated version
    no2_info_v2 = get.tropomi.prof(receptor, 'NO2', tropomi.fn = no2_fn)
    no2_info_v1$tno2_ppb_v2 = no2_info_v2$no2_vcd_tropo / air_vcd_tropo * 1e9
    no2_info_v1$tno2_uncert_ppb_v2 = no2_info_v2$no2_vcd_tropo_uncert / air_vcd_tropo * 1e9

    # correct for diff in ZSFC between sensor and particle ZSFC
    ak_prof = as.data.frame(no2_info_v2[c('ak_tropo', 'lower_pres')])
    t.psfc = no2_info_v2$tropomi_psfc
    t.zsfc = no2_info_v2$tropomi_zsfc	

    if ('time' %in% colnames(p_chem))
        p_chem = p_chem %>% filter(time == max(time))
        
	p.zagl.corr = p_chem$zagl + p_chem$zsfc - t.zsfc
	p.pres = p_chem$pres 
    
	# use satellite surface altitude and pres to calc ZAGL that match TROPOMI 
	# use particle-level pressure and ASL to calculate the coeffient
	nls = stats::nls(p.pres ~ t.psfc * exp(a * (-p.zagl.corr)), 
					 start = list(a = 1E-4))		
	a.tropomi = coef(nls)[[1]]    # a stands for g/RTv_mean using TROPO atmos-X

    # calc modeled tNO2 for version 2
    p_wgt = p_chem %>% rename(ak_tno2_v1 = ak_tno2) %>% 
            mutate(lower_pres = t.psfc * exp(a.tropomi*(-xhgt)), 
                   ak_tno2 = approx(ak_prof$lower_pres, ak_prof$ak_tropo, 
                                    xout = lower_pres, rule = 2)$y) %>% 
            dplyr::select(-c('lower_pres'), -ends_with('_wgt'))%>%
            mutate(
                p_nox_mix_wgt = p_nox_mix * ak_tno2 * air_vcd_stilt / air_vcd_tropo,

                p_nox_nomix_wgt = p_nox_nomix * ak_tno2 * air_vcd_stilt / air_vcd_tropo,

                p_nox_nochem_wgt = p_nox_nochem * ak_tno2 * air_vcd_stilt / air_vcd_tropo,

                p_no2_nochem_wgt = p_nox_nochem * 0.74 * ak_tno2 * air_vcd_stilt / air_vcd_tropo, 

                # p_no2_mix_wgt with or without chem still in ppm 
                # but weighted by AK and PWF
                p_no2_mix_wgt = p_no2_mix * ak_tno2 * air_vcd_stilt / air_vcd_tropo, 

                p_no2_nomix_wgt = p_no2_nomix * ak_tno2 * air_vcd_stilt / air_vcd_tropo
            ) %>% summarise_all(mean) %>% 
            dplyr::select(contains('_wgt'), contains('ak_tno2')) %>% 
            mutate_at(vars(matches('p_no')), function(x) x * 1e3) 

    # plot(ak_prof$lower_pres, ak_prof$ak_tropo)
    # points(p_wgt$pres, p_wgt$ak_tno2, col = 'red', cex = 0.1)
    # points(p_wgt$pres, p_wgt$ak_tno2_v1, col = 'blue', cex = 0.1)

    return(list(p_sim = p_wgt, no2_info = no2_info_v1))
}

