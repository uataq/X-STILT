# script to get overpass time string from OCO or TROPOMI (real simulation)
# or check receptor time from self-provided 'recp_fn' file (ideal sim)
# DW, 07/05/2021

# plotTF  = FALSE             # T: plot OCO XCO2 and SIF on maps and lat series
# nfTF = TRUE                 # T: also count overpasses over near-field region
# searchTF = FALSE            # T: loop over all OCO overpasses and find all 
                              #  the overpasses that pass by your site
if (F) {
    nfTF = T 
    searchTF = F
    date_range = c('20140101', '20221231')
    obs_ver = 'V10r'
    qfTF = T 
    plotTF = T
}

get_timestr = function(site, lon_lat, obs_sensor, obs_ver = 'V10r', obs_path, 
                       store_path, recp_fn = NULL, plotTF = F, qfTF = F, 
                       nfTF = T, nf_dlon = 0.5, nf_dlat = 0.5, sif_path = NULL, 
                       searchTF = F, date_range = c('20140101', '20211231')) {

    timestr = NA        # initialize
    if ( is.na(obs_sensor) )  {                 # NA for ideal simulation

        cat('\n\nTHIS IS AN IDEAL X-SIMULATION (obs_sensor is NA), no satellite data will be used.')
        cat('\n\nChecking .csv/.txt file indicated by @param recp_fn...\n')
        
        if (is.null(recp_fn)) 
            stop('NO receptor file found, please check recp_fn\n')
        if (!file.exists(recp_fn)) 
            stop('NO receptor file found, please check recp_fn\n')
        recp_info = read.csv(recp_fn, sep = ',')
        
        # check if txt file contains lati and long with correct column names
        colTF = colnames(recp_info) %in% c('lati', 'long')
        if (FALSE %in% colTF) stop(paste0('Incorrect column names for',recp_fn))
        if ('timestr' %in% colnames(recp_info)) {
            cat('Reading time string from the receptor file as receptor time\n')
            timestr = min(unique(recp_info))
        } 

    } else if ( grepl('OCO', obs_sensor) ) {        # simulations using OCO-2/3
        
        # qfTF only controls whether data points will be used for plotting
        oco_track = get.site.track(site, oco.sensor = obs_sensor, 
                                   oco.ver = obs_ver, oco.path = obs_path, 
                                   searchTF = searchTF, 
                                   date.range = date_range, 
                                   thred.count.per.deg = 100, 
                                   lon.lat = lon_lat, nfTF = nfTF, 
                                   nf.dlon = nf_dlon, nf.dlat = nf_dlat, 
                                   thred.count.per.deg.nf = 50, 
                                   plotTF, store.path = store_path, 
                                   sif.path = sif_path, qfTF)   
                                   
        all.timestr = oco_track$timestr; print(oco_track)

        # *** NOW pick the time string that you'd like to work on...
        cat('Enter the row index from the df above for your target overpass: ')
        tt = as.numeric(readLines('stdin', n = 1))
        timestr = all.timestr[tt]

    } else if ( grepl('TROPOMI', obs_sensor) ) {

        # as TROPOMI has daily coverage, no need to look for overpass info for your site
        # *** NOW choose the timestr that you'd like to work on..
        cat('\n\nSince TROPOMI possesses daily coverage, no need to search for overpass info\n')
        cat('Instead please enter the TROPOMI overpass time in format of YYYYMMDD:\n')
        timestr = as.numeric(readLines('stdin', n = 1))

    } else stop('Incorrect @param obs_sensor, please check...')

    # if timestr is still not found, mannually insert one
    #if ( is.na(timestr) ) {
    #    cat('Please enter a time string to work on in form of YYYYMMDDHH:\n')
    #    timestr = readLines('stdin', n = 1)#as.numeric()
    #}   

    #cat(paste('Working on:', timestr, 'for site/region:', site, '...\n\n'))
    #return(timestr)
}