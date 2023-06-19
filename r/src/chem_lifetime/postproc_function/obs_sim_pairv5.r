
# ---------------------------------------------
obs_sim_pairv5 = function(site, timestr, out_path, rds_path, rds_patt, met, 
                          overwriteTF = T, transerrTF = F, site_lon = NULL, site_lat = NULL) {
    
    fn = file.path(rds_path, paste0(site, '_', timestr, '_', rds_patt, '_', 
                                    met, '.rds'))
    if ( transerrTF) fn = file.path(rds_path, paste0(site, '_', timestr, '_', 
                                    rds_patt, '_', met, '_transerr.rds'))
    
    # ---------------------------------------------
    if ( !file.exists(fn) | overwriteTF ) {
        rds_fns = list.files(out_path, paste0(rds_patt, '.rds'), 
                             full.names = T, recursive = T)
        if (is.null(rds_fns)) { cat('NO rds files found...\n'); return() }
        if (length(rds_fns) == 0) { cat('NO rds files found...\n'); return() }
        
        rds_info = strsplit.to.df(basename(rds_fns)) %>%
                   dplyr::select(time = V1, lon = V2, lat = V3) %>% 
                   mutate_if(is.character, as.numeric) %>% mutate(fn = rds_fns)
        rds_fns = rds_info$fn; print(length(rds_fns))
        
        # ---------------------------------------------
        # go over each file
        x_df = NULL 
        for ( x in 1 : length(rds_fns) ) {

            if (x %% 20 == 0) 
            cat(paste('# ---', signif(x / length(rds_fns), 3) *100,'% --- #\n'))
            dat = readRDS(rds_fns[x])
            recp_lon = strsplit.to.df(basename(rds_fns[x]))$V2
            recp_lat = strsplit.to.df(basename(rds_fns[x]))$V3

            # calculate the column enhancement by subtracting the initial cond 
            # DW, 04/04/2022 
            tmp_recp = dat$p_recp #%>% na.omit()
            tmp_info = dat$no2_info 
            tmp_chem = dat$p_chem
            air_vcd_stilt = tmp_info$air_vcd_stilt
            air_vcd_tropo = tmp_info$air_vcd_tropo
            
            # convert TNO2 uncertainty to ppb
            tmp_info$tno2_uncert_ppb = 
                tmp_info$no2_vcd_tropo_uncert / air_vcd_tropo * 1e9

            tmp_info = tmp_info[names(tmp_info) %in% 
                                c('tropomi_lat', 'tropomi_lon', 'tropomi_lats', 
                                  'tropomi_lons', 'tno2_ppb', 'tno2_uncert_ppb',
                                  'air_vcd_tropo', 'air_vcd_stilt')]
            tmp_info$corner = seq(1, 4)

            # ---------------------------------------------
            if ( !grepl('emiss', rds_patt) & !grepl('tau', rds_patt) ) 
                tmp_recp = tmp_recp %>% 
                           mutate(sim_tno2ff_mix = (p_no2_mix - no2_ic) * 
                                  ak_tno2 * air_vcd_stilt / air_vcd_tropo, 
                                  
                                  sim_tno2ff_nomix = (p_no2_nomix - no2_ic) * ak_tno2 * air_vcd_stilt / air_vcd_tropo)

            # ---------------------------------------------
            tmp_x = tmp_recp %>% 
                    dplyr::select(mean_foot = foot_wgt, ends_with('_wgt')) %>% 
                    summarise_all(list(~mean(., na.rm = TRUE))) %>% 

                    # convert ppm to ppb for NOx and CO
                    mutate_at(vars(matches('p_no')), function(x) x * 1e3) %>% 
                    mutate_at(vars(matches('p_co_')), function(x) x * 1e3) %>% 

                    # NOx, NO2, CO have units in ppb
                    rename(sim_tnox_nochem = p_nox_nochem_wgt, 
                           sim_tnox_mix = p_nox_mix_wgt, 
                           sim_tnox_nomix = p_nox_nomix_wgt,

                           sim_tno2_nochem = p_no2_nochem_wgt,
                           sim_tno2_mix = p_no2_mix_wgt, 
                           sim_tno2_nomix = p_no2_nomix_wgt, 

                           #sim_xco_mix = p_co_mix_wgt, 
                           #sim_xco_nomix = p_co_nomix_wgt, 
                           sim_xco2_mix = p_co2_mix_wgt, 
                           sim_xco2_nomix = p_co2_nomix_wgt)

            tmp_df = as.data.frame(tmp_info) %>% mutate(indx = x) %>% 
                     bind_cols(tmp_x) %>% 
                     dplyr::select(indx, corner, lon = tropomi_lon, 
                                   lat = tropomi_lat, lons = tropomi_lons, 
                                   lats = tropomi_lats, mean_foot, 
                                   obs_tno2 = tno2_ppb, 
                                   obs_tno2_uncert = tno2_uncert_ppb, 
                                   starts_with('sim_')) %>% arrange(corner) %>% 
                     mutate(recp_lon = as.numeric(recp_lon), 
                            recp_lat = as.numeric(recp_lat))

            x_df = rbind(x_df, tmp_df)
        }   # end for

        saveRDS(x_df, file = fn)
    } else {
        cat('Found rds file, loading results...\n')
        print(fn)
        x_df = readRDS(fn)
    }

    return(x_df)
}

