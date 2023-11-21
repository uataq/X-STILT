
correct_emiss_epa = function(emiss_cp, epa_name, epa_loc, epa_tz, epa_fn, 
                             epa_species = c('CO2', 'NOx', 'SO2')) {

    # cannot use the max value as pp emissions, need to locate the EDGAR pixel that PP falls into 
    epa_loc = load_epa_pp(epa_name)
    
    # grab the initial NO emission and its coordinates from EDGAR
    emiss_pp = raster::extract(emiss_cp, epa_loc)

    # convert hourly total pounds to mean flux on 0.1deg in umol m-2 s-1 
    # 1 pound = 0.453592 kg; 1E3 for kg to g; 46 for g to mol; 
    # 1E6 for mol to umol; 3600 for hourly to second; 
    # 0.1*111 for EDGAR grid spacing; 1E6 for km2 to m2
    # !!! EPA record data in local standard time
    hr_df = read.csv(file = epa_fn, row.names = NULL, stringsAsFactors = F) %>%
            filter(grepl(epa_name, FACILITY_NAME)) %>% 
            mutate(
                NOx = NOX_MASS..lbs. * 0.453592 * 1E3 / 46 * 1E6 / 3600 / (0.1 * cos(epa_loc$lat * pi / 180) * 111 * 0.1 * 111) / 1E6, 

                CO2 = CO2_MASS..tons. * 907.185 * 1E3 / 44 * 1E6 / 3600 / (0.1 * cos(epa_loc$lat * pi / 180) * 111 * 0.1 * 111) / 1E6, 

                # EPA always report local standard time
                # correct the tz of epa_df$date to pp-specific tz, 
                # otherwise it'll be the user-specific timezone (e.g., PDT)
                Date = paste(OP_DATE, formatC(OP_HOUR, width = 2,flag = 0)),
                date_lst = as.POSIXct(Date, epa_tz, format = '%m-%d-%Y %H')
            ) 

    # sum by facility ID rather than unit ID
    if (epa_species == 'NOx') {
        epa_df = hr_df %>% 
                 dplyr::select(FACILITY_NAME, date_lst, FAC_ID, NOx) %>%
                 group_by(FACILITY_NAME, date_lst, FAC_ID) %>% 
                 dplyr::summarise(NOx = sum(NOx), .groups = 'drop') %>% 

                 # calc scaling factor for correcting NOx emission for PP
                 ungroup() %>% mutate(sf_pp_no = NOx / emiss_pp) %>% 
                 dplyr::select(date_lst, sf_pp_no, epa_pp = NOx)
    
        #cat(paste('load_eno(): on average, scale ENO by', 
        #          signif(mean(epa_df$sf_pp_no, na.rm = T), 3), 
        #          'according to EPA\n\n'))

    } else if (epa_species == 'CO2') {

        epa_df = hr_df %>% 
                 dplyr::select(FACILITY_NAME, date_lst, FAC_ID, CO2) %>%
                 group_by(FACILITY_NAME, date_lst, FAC_ID) %>% 
                 dplyr::summarise(CO2 = sum(CO2), .groups = 'drop') %>% 

                 # calc scaling factor between EDGAR and EPA
                 ungroup() %>% mutate(sf_pp_co2 = CO2 / emiss_pp) %>% 
                 dplyr::select(date_lst, sf_pp_co2, epa_pp = CO2)
        
        #cat(paste('load_eco2(): on average, scale ECO2 by', 
        #          signif(mean(epa_df$sf_pp_co2, na.rm = T), 3), 
        #          'according to EPA\n\n'))

    } else cat('correct_emiss_epa(): wrong EPA species...\n')

    return(list(epa_df = epa_df, emiss_pp = emiss_pp))
}



# !!! this function is outdated...
load_epa_hrly = function(epa_fn, epa_name, epa_tz, epa_loc) {

    hr_df = read.csv(file = epa_fn, row.names = NULL, stringsAsFactors = F) %>%
            filter(grepl(epa_name, FACILITY_NAME)) %>% 
            mutate(
                NOx = NOX_MASS..lbs. * 0.453592 * 1E3 / 30 * 1E6 / 3600 / (0.1 * cos(epa_loc$lat * pi / 180) * 111 * 0.1 * 111) / 1E6, 

                CO2 = CO2_MASS..tons. * 907.185 * 1E3 / 44 * 1E6 / 3600 / (0.1 * cos(epa_loc$lat * pi / 180) * 111 * 0.1 * 111) / 1E6, 

                # EPA always report local standard time
                # correct the tz of epa_df$date to pp-specific tz, 
                # otherwise it'll be the user-specific timezone (e.g., PDT)
                Date = paste(OP_DATE, formatC(OP_HOUR, width = 2, flag = 0)),
                date_lst = as.POSIXct(Date, epa_tz, format = '%m-%d-%Y %H')
            ) %>% 

            # sum up all units to facility level
            dplyr::select(FACILITY_NAME, date_lst, FAC_ID, NOx, CO2) %>%
            group_by(FACILITY_NAME, date_lst, FAC_ID) %>% 
            dplyr::summarise(NOx = sum(NOx), CO2 = sum(CO2), 
                            .groups = 'drop') %>% ungroup()

    return(hr_df)
}
