
# output - particle
# no2_fn - updated NO2 file
# p_chem - particle with chemical calcs
# tmp_info - no2 profile from old calcs
update_ak = function(output, no2_fn, p_chem, tmp_info) {    
    
    # get NO2 profile from updated version
    no2_all  = get.wgt.tropomi.func(output, no2_fn, 'NO2')
    no2_prof = no2_all$combine.prof %>% filter(stiltTF) %>% 
               dplyr::select(indx, ak_tno2 = ak.norm)
    no2_info = no2_all$tropomi.info     # version 2.3

    # calc tNO2 for version 2
    air_vcd_stilt = tmp_info$air_vcd_stilt
    air_vcd_tropo = tmp_info$air_vcd_tropo
    tmp_info$tno2_ppb_v2 = no2_info$no2_vcd_tropo / air_vcd_tropo * 1e9
    tmp_info$tno2_uncert_ppb_v2 = no2_info$no2_vcd_tropo_uncert / air_vcd_tropo * 1e9

    p_wgt = p_chem %>% rename(ak_tno2_v1 = ak_tno2) %>%  
            left_join(no2_prof, by = 'indx') %>% 
            filter(time == max(time)) %>%
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

    return(list(p_sim = p_wgt, tmp_info = tmp_info))
}

