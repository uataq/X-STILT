
download.met.arl = function(timestr, met_file_format, nhrs, met_path, 
                            met = c('wrf27km', 'hrrr', 'gfs0p25')[2]) {
    
    options(timeout = max(1500, getOption('timeout')))   # for download.file

    # quickly check if met fields exist 
    t_start = as.POSIXct(timestr, 'UTC', format = '%Y%m%d%H')
    met_fns = find_met_files(t_start, met_file_format, n_hours = nhrs, met_path)
    
    t_end = t_start + nhrs * 3600
    t_met = unique(format(seq(t_end, t_start, by = 'hour'), format = '%Y%m%d'))
    n_met_min = length(t_met)
    
    if (met == 'hrrr') n_met_min = n_met_min * 4
    met_fns = do.call(c, lapply(t_met, function(x) met_fns[grepl(x, met_fns)] ))

    # downloading met fields if not found
    if ( length(met_fns) < n_met_min ) {
        cat('download.met.arl(): missing metfiles...downloading now...\n')

        arl_path = file.path('ftp://arlftp.arlhq.noaa.gov/archives', met)
        if (met == 'wrf27km') 
            urls = file.path(arl_path, substr(timestr, 1, 4), 
                             paste0('wrfout_d01_', timestr, '.ARL'))
        
        if (met == 'gfs0p25') {
            if (timestr <= '2019050100' & timestr >= '2016051300') 
            arl_path = gsub('gfs0p25', 'gfs0p25.v1', arl_path)
            urls = file.path(arl_path, paste0(t_met, '_', met))
        }

        if (met == 'hrrr') {
            # ftp://arlftp.arlhq.noaa.gov/archives/hrrr.v1 for data before July 2019
            if (timestr < '2019070100') {
                arl_path = gsub('hrrr', 'hrrr.v1', arl_path)
                hrrr_strs = c('00z', '06z', '12z', '18z')
                hrrr_df = expand.grid(t_met, hrrr_strs)
                urls = file.path(arl_path, paste('hysplit', hrrr_df$Var1, 
                                              hrrr_df$Var2, 'hrrra', sep = '.'))

            } else {
                hrrr_strs = c('00-05', '06-11', '12-17', '18-23')
                hrrr_df = expand.grid(t_met, hrrr_strs)
                urls = file.path(arl_path, paste(hrrr_df$Var1, hrrr_df$Var2, 
                                                 met, sep = '_'))
            }
        }   # end if

        for (url in urls) 
            download.file(url, destfile = file.path(met_path, basename(url)))
    } # end if

    print(met_fns)
    return(n_met_min)
}