
correct_emiss_epav2 = function(emiss_cp, epa_name, epa_tz, epa_fn, 
                               epa_species = c('CO2', 'NOx', 'SO2')) {

    # cannot use the max value as pp emissions, need to locate the EDGAR pixel that PP falls into 
    epa_loc = load_epa_pp(epa_name)
    
    # grab the initial NO emission and its coordinates from EDGAR
    emiss_pp = raster::extract(emiss_cp, epa_loc)

    # convert hourly total pounds to mean flux on 0.1deg in umol m-2 s-1 
    # 1 pound = 0.453592 kg; 1E3 for kg to g; 46 for g to mol; 
    # 1E6 for mol to umol; 3600 for hourly to second; 
    # 0.1*111 for EDGAR grid spacing; 1E6 for km2 to m2

    # !!! EPA record data in local standard time, so the LST - UTC offset is 
    # constant throughout the year. e.g., for Chicago, dtz = -6 
    dtz = lutz::tz_offset('2020-01-01', epa_tz)$utc_offset_h
    hr_df = read.csv(file = epa_fn, row.names = NULL, stringsAsFactors = F) %>%
            
            # also remove NA NOx first before summing up unit-level emissions
            filter(grepl(epa_name, FACILITY_NAME), !is.na(NOX_MASS..lbs.)) %>% 
            mutate(
                NOx = NOX_MASS..lbs. * 0.453592 * 1E3 / 46 * 1E6 / 3600 / (0.1 * cos(epa_loc$lat * pi / 180) * 111 * 0.1 * 111) / 1E6, 
                
                CO2 = CO2_MASS..tons. * 907.185 * 1E3 / 44 * 1E6 / 3600 / (0.1 * cos(epa_loc$lat * pi / 180) * 111 * 0.1 * 111) / 1E6, 

                # EPA always report local standard time
                timestr_lst = paste0(substr(OP_DATE, 7, 10), 
                                     substr(OP_DATE, 1, 2), 
                                     substr(OP_DATE, 4, 5), 
                                     formatC(OP_HOUR, width = 2, flag = 0)), 
                
                # although LST, we state in "UTC" since we are correcting for the UTC offset
                datestr_utc = as.POSIXct(timestr_lst, 'UTC', 
                                         format = '%Y%m%d%H') - dtz * 3600,     

                # rewrite as character
                timestr_utc = format(datestr_utc, format = '%Y%m%d%H')
            ) %>% 

            # sum up individual units to facility total
            group_by(FACILITY_NAME, FAC_ID, timestr_lst, timestr_utc) %>% 
            summarise_if(is.numeric, sum) %>% ungroup() %>% 
            dplyr::select(FACILITY_NAME, timestr_lst, timestr_utc, NOx, CO2)

    # calc scaling factor for correcting NOx emission for PP
    if (epa_species == 'NOx') {
        epa_df = hr_df %>% mutate(sf_pp_no = NOx / emiss_pp) %>% 
                 dplyr::select(timestr_lst, timestr_utc, sf_pp_no, epa_pp = NOx)
    
    } else if (epa_species == 'CO2') {
        epa_df = hr_df %>% mutate(sf_pp_co2 = CO2 / emiss_pp) %>% 
                dplyr::select(timestr_lst, timestr_utc, sf_pp_co2, epa_pp = CO2)

    } else cat('correct_emiss_epa(): wrong EPA species...\n')

    return(list(epa_df = epa_df, emiss_pp = emiss_pp))
}

