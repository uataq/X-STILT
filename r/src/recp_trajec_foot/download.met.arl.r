

download.met.arl = function(timestr, nhrs, run_trajec = F, met_path, 
                            met = c('wrf27km', 'hrrr', 'gfs0p25')[2], 
                            met_file_format = NULL) {
    
    # if only day is given, need to download more met files just in case
    if ( nchar(timestr) == 8) timestr = paste0(timestr, c('00', '23'))
    t_start = as.POSIXct(timestr, 'UTC', format = '%Y%m%d%H')
    t_end = t_start + nhrs * 3600
    t_hrs = seq(min(t_end, t_start), max(t_end, t_start), by = 'hour')
    t_day = unique(as.Date(t_hrs, format = '%Y%m%d%H'))

    # ------------------------------------------------------------------------
    # initialize temporal resolution of met files from NOAA ARL 
    #arl_path = file.path('ftp://arlftp.arlhq.noaa.gov/archives', met)
    arl_path = file.path('https://www.ready.noaa.gov/data/archives', met)
    if (met == 'wrf27km') { 
        t_met = t_hrs
        met_file_format = 'wrfout_d01_%Y%m%d%H.ARL'
        arl_path = file.path(arl_path, 'avg', unique(substr(t_hrs, 1, 4)))
    }

    if (met == 'hrrr') {
        if ( max(t_day) < as.Date('2019-07-01') ) {
            met_file_format = 'hysplit.%Y%m%d.%Hz.hrrra'
            arl_path = gsub('hrrr', 'hrrr.v1', arl_path)
        } else met_file_format = '%Y%m%d_%H'
        t_met = unique(lubridate::floor_date(t_hrs, '6 hours'))
    }

    # no GFS0p25 data before May 13, 2016
    # GFS0p25 version 1 between May 13, 2016 and June 12, 2019
    # GFS0p25 default after June 13, 2019
    if (met == 'gfs0p25') { 
        if ( min(as.Date(t_hrs)) < as.Date('2016-05-13') ) 
            stop('download.met.arl(): NO GFS0p25 from NOAA ARL before May 13, 2016...\n')

        if ( min(as.Date(t_hrs)) >= as.Date('2016-05-13') & 
             max(as.Date(t_hrs)) <= as.Date('2019-06-12')) 
            arl_path = gsub('gfs0p25', 'gfs0p25.v1', arl_path)
        
        t_met = unique(lubridate::floor_date(t_hrs, 'day'))
        met_file_format = '%Y%m%d_gfs0p25'
    }

    if ( is.null(met_file_format) ) 
        stop('download.met.arl(): User needs to provide data format for meteo files...\n')

    # form filename that matches NOAA ARL
    fns = strftime(t_met, tz = 'UTC', format = met_file_format)

    # hrrr.v1 for data before July 2019; hrrr for data after July 2019
    if ( met == 'hrrr' & min(t_day) > as.Date('2019-07-01') ) 
        fns = paste0(fns, '-', substr(t_met + 5 * 3600, 12, 13), '_hrrr')
    cat('download.met.arl(): require the following met files:\n'); print(fns)

    # ------------------------------------------------------------------------
    #t_met = unique(format(seq(t_end, t_start, by = 'hour'), format = '%Y%m%d'))
    n_met_min = length(t_met)
    ava_fns = basename(find_met_files(t_start, met_file_format, nhrs, met_path))
    
    find_fns = do.call(c, lapply(fns, function(x) ava_fns[grepl(x, ava_fns)] ))
    miss_fns = do.call(c, lapply(fns, function(x) {
        if ( length(ava_fns[grepl(x, ava_fns)]) == 0 ) return(x) } ))

    # ------------------------------------------------------------------------
    # downloading met fields if not found
    if ( length(ava_fns) < n_met_min & run_trajec ) {
        cat('download.met.arl(): missing metfiles...downloading now...\n')

        urls = file.path(arl_path, miss_fns)
        for (url in urls) 
            download.file(url, destfile = file.path(met_path, basename(url)))
    } else cat('Found meteo files...\n') # end if

    return(list(n_met_min = n_met_min, met_file_format = met_file_format))
}