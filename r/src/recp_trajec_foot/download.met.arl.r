
# fix timestring, always convert it to character format
# add a line to check if any met files exist; if so, skip it, DW, 10/20/2024
# add n_met_min back, in case users prefer their own fields, DW, 10/22/2024
download.met.arl = function(timestr, met_file_format, nhrs, met_path, 
                            met = c('wrf27km', 'hrrr', 'gfs0p25')[2], 
                            met_file_tres = c('1 hour', '6 hours', '1 day')[2], 
                            run_trajec = F, n_met_min = NA, selfTF = FALSE) {

    # quickly check if met fields exist 
    for (tmp_timestr in timestr) {
        
        if ( nchar(tmp_timestr) == 8) ff = '%Y%m%d'
        if ( nchar(tmp_timestr) == 10) ff = '%Y%m%d%H'
        t_start = as.POSIXct(as.character(tmp_timestr), 'GMT', format = ff)
        met_fns = find_met_files(t_start, n_hours = nhrs, met_path, 
                                 met_file_format, met_file_tres)

        # according to ARL, 1hr for WRF27km, 6hr for HRRR, 1d for GFS
        if (!selfTF) {
            t_end = t_start + nhrs * 3600
            t_min = min(c(t_start, t_end))
            t_max = max(c(t_start, t_end))
            t_hrs = seq(t_min, t_max, by = 'hour')
            t_dys = unique(trunc(t_hrs, units = 'days')) + c(-1, 1) * 3600 * 24
            if (met == 'wrf27km') arl_format = '%Y%m%d%H'
            if (met == 'hrrr') 
                arl_format = ifelse(tmp_timestr < '2019070100', 
                                    '%Y%m%d.%Hz', '%Y%m%d_%H') 
        
            if (met == 'gfs0p25') arl_format = '%Y%m%d'
            t_bound = seq(min(t_dys), max(t_dys), by = met_file_tres)
            t_find_time = unique(t_bound[findInterval(t_hrs, t_bound)])
            t_met = strftime(t_find_time, format = arl_format, tz = 'GMT')
            if (is.na(n_met_min)) n_met_min = length(t_met)
            met_fns = do.call(c, lapply(t_met, 
                                    function(x) met_fns[grepl(x, met_fns)] ))
        }   # end if for downloading files from ARL
        
        
        # downloading met fields if not found ----------------------------
        if ( length(met_fns) < n_met_min & run_trajec & !selfTF) {
            cat('download.met.arl(): missing metfiles...downloading now...\n')
            
            arl_path = file.path('ftp://arlftp.arlhq.noaa.gov/archives', met)
            if (met == 'wrf27km') 
                urls = file.path(arl_path, 'avg', substr(t_met, 1, 4), 
                                paste0('wrfout_d01_', t_met, '.ARL'))
            
            if (met == 'gfs0p25') {
                if (tmp_timestr <= '2019050100' & tmp_timestr >= '2016051300') 
                arl_path = gsub('gfs0p25', 'gfs0p25.v1', arl_path)
                urls = file.path(arl_path, paste0(t_met, '_', met))
            }

            if (met == 'hrrr') {
                # ftp://arlftp.arlhq.noaa.gov/archives/hrrr.v1 for data before July 2019
                if ( tmp_timestr < '2019070100' ) {
                    arl_path = gsub('hrrr', 'hrrr.v1', arl_path)
                    urls = file.path(arl_path, paste0('hysplit.', t_met, '.hrrra'))
                } else {
                    hh_end = formatC(as.numeric(substr(t_met, 10, 11)) + 5, 
                                    width = 2, flag = 0) 
                    urls = file.path(arl_path, paste0(t_met, '-', hh_end, '_', met))
                }
            }   # end if

            for (url in urls) {
                dest_fn = file.path(met_path, basename(url))
                if (!file.exists(dest_fn)) 
                    download.file(url, destfile = dest_fn)
            }
            
        } else if ( length(met_fns) < n_met_min & run_trajec ) {
            stop('Since you chose to your own met fields, insufficient met fields found...\n')
        } # end if

        print(met_fns)
    }

    return(n_met_min)
}