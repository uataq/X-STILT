# ---------------------------------------------------------------------------- #
prep_sim_obs = function(pr_edgar, bg, pr_epa = NULL) {
    
    # merging and add NO2-NOx ratio
    colnms = c('corner', 'lons', 'lats', 'lon', 'lat')
    pr = pr_edgar %>% dplyr::select(corner, lat, lon, lats, lons, obs_tno2,
                                    obs_tno2_uncert, #obs_tno2_ff, 
                                    edgar_tno2_nochem = sim_tno2_nochem_ppb, 
                                    edgar_tno2_nomix = sim_tno2_nomix_ppb, 
                                    edgar_tno2_mix2 = sim_tno2_mix2_ppb,
                                    edgar_tno2_nomix2 = sim_tno2_nomix2_ppb,
                                    edgar_xco2 = sim_xco2_ppm, 
                                    contains('obs_xco2_ff')) %>% 
        mutate(edgar_ratio = edgar_tno2_nomix2 / edgar_xco2)

    if ( !is.null(pr_epa) ) 
        pr = pr_epa %>% dplyr::select(corner, lat, lon, lats, lons, 
                                      epa_tno2_nochem = sim_tno2_nochem_ppb, 
                                      epa_tno2_nomix = sim_tno2_nomix_ppb, 
                                      epa_tno2_mix2 = sim_tno2_mix2_ppb,
                                      epa_tno2_nomix2 = sim_tno2_nomix2_ppb,
                                      epa_xco2 = sim_xco2_ppm) %>% 
             right_join(pr, by = colnms) %>%                       
             mutate(epa_ratio = epa_tno2_nomix2 / epa_xco2)

    # choose satellite polygons with some XCO2 enhancements
    pr = pr %>% mutate(enhTF = edgar_tno2_nomix2 > bg | obs_tno2 > bg + 0.1, 
                       polygon = findInterval(seq(1, nrow(pr)), 
                                              which(pr$corner == 1)))
    
    pm = pr %>% dplyr::select(-ends_with('_ratio')) %>% 
         reshape2::melt(id.vars = c(colnms, 'enhTF', 'polygon')) %>% 
         mutate(variable = toupper(variable))

    return(list(pr = pr, pm = pm))
}



# ---------------------------------------------------------------------------- #
prep_sim_obsv3 = function(pr_edgar = NULL, pr_epa = NULL, pr_odiac = NULL) {
    
    # merging and add NO2-NOx ratio
    colnms = c('corner', 'lons', 'lats', 'lon', 'lat')
    pr = data.frame(corner = NA, lons = NA, lats = NA, lon = NA, lat = NA)
    if ( !is.null(pr_edgar) ) 
        pr = pr_edgar %>% dplyr::select(corner, lat, lon, lats, lons, 
                                        obs_tno2_tot = obs_tno2,
                                        obs_tno2_uncert, 
                                        contains('obs_tno2_ff'),
                                        contains('obs_xco2_ff'), 
                                        edgar_tno2_nochem = sim_tno2_nochem, 
                                        edgar_tno2_nomix = sim_tno2_nomix, 
                                        edgar_tno2_mix = sim_tno2_mix,
                                        edgar_tno2_ff = sim_tno2_ff, 
                                        edgar_xco2 = sim_xco2_mix) %>% 
            left_join(pr, by = colnms) %>%
            mutate(edgar_ratio = edgar_tno2_nomix / edgar_xco2)

    if ( !is.null(pr_epa) & is.null(pr_edgar)) {
        pr = pr_epa %>% dplyr::select(corner, lat, lon, lats, lons, 
                                      obs_tno2_tot = obs_tno2,
                                      obs_tno2_uncert, 
                                      contains('obs_tno2_ff'),
                                      contains('obs_xco2_ff'), 
                                      epa_tno2_nochem = sim_tno2_nochem, 
                                      epa_tno2_nomix = sim_tno2_nomix, 
                                      epa_tno2_mix = sim_tno2_mix,
                                      epa_tno2_ff = sim_tno2_ff, 
                                      epa_xco2 = sim_xco2_mix) %>% 
             left_join(pr, by = colnms) %>%                       
             mutate(epa_ratio = epa_tno2_nomix / epa_xco2)

    } else if ( !is.null(pr_epa) & !is.null(pr_edgar) ) {
        pr = pr_epa %>% dplyr::select(corner, lat, lon, lats, lons, 
                                      epa_tno2_nochem = sim_tno2_nochem, 
                                      epa_tno2_nomix = sim_tno2_nomix, 
                                      epa_tno2_mix = sim_tno2_mix,
                                      epa_tno2_ff = sim_tno2_ff, 
                                      epa_xco2 = sim_xco2_mix) %>% 
             left_join(pr, by = colnms) %>%                       
             mutate(epa_ratio = epa_tno2_nomix / epa_xco2)
    }
        

    if ( !is.null(pr_odiac) ) 
        pr = pr_odiac %>% dplyr::select(corner, lat, lon, lats, lons, 
                                        odiac_tno2_nochem = sim_tno2_nochem, 
                                        odiac_tno2_nomix = sim_tno2_nomix, 
                                        odiac_tno2_mix = sim_tno2_mix,
                                        odiac_xco2 = sim_xco2_mix) %>% 
             left_join(pr, by = colnms) %>%                       
             mutate(odiac_ratio = odiac_tno2_nomix / odiac_xco2)

    # choose satellite polygons with some XCO2 enhancements
    pr = pr %>% mutate(polygon = findInterval(seq(1, nrow(pr)), 
                                              which(pr$corner == 1)))
    
    return(pr)
}