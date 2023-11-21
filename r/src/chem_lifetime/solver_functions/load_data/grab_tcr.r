# function to read OH and NO2 fields from Chemical reanalysis product from JPL, 
# DW, 09/28/2020 

if (F) {
    tcr.fn = 'mon_emi_nox_anth_2019.nc'
    #xmn = 109; xmx = 110.5; ymn = 39.5; ymx = 41.2
    xmn = -112.5; xmx = -111.3; ymn = 32.5; ymx = 34
    tmn = '2019-01-01 00:00:00'
    tmx = '2019-12-31 00:00:00'
    tcr.varname = 'nox'

}

# allow for extracting variables other than OH
grab_tcr = function(tcr.fn, xmn, xmx, ymn, ymx, tmn, tmx, 
                    tcr.varname = c('oh', 'nox', 'no', 'no2')[2], 
                    tres = c('hourly', 'monthly')[1], levelTF = T, 
                    tformat = '%Y-%m-%d %H:%M:%S') {

    # if trajec time went beyond the latest available year
    # use the latest available year by modifying the trajec year
    if ( as.numeric(substr(tmn, 1, 4)[1]) > 2019 ) {
        tmn = paste0('2019', substr(as.character(tmn), 5, 
                                    nchar(as.character(tmn))) )
        tmx = paste0('2019', substr(as.character(tmx), 5, 
                                    nchar(as.character(tmx))) )
    }   # end if

    tmn = as.POSIXct(tmn, 'UTC', format = tformat)
    tmx = as.POSIXct(tmx, 'UTC', format = tformat)
    
    # load data ---------------------------------------------------------------
    dat = nc_open(tcr.fn)
    lat = ncvar_get(dat, 'lat')                # uneven latitudes
    lon = ncvar_get(dat, 'lon')                # degree east, [0, 360]
    lon[lon > 180] = lon[lon > 180] - 360      # convert to [-180, 180]

    if ( tres == 'hourly' ) {    # hourly data
        mins = ncvar_get(dat, 'time')          # minutes since 2019-01-01 01:00
        yr   = substr(tmn, 1, 4)[1]
        times = as.POSIXct(mins * 60, 'UTC', origin = paste0(yr,'-01-01 01:00'))
        time.indx = which(times >= tmn & times <= tmx)

    } else if (tres == 'monthly') {    
        
        # monthly data, one file for all years from 2005
        # as.POSIXlt()$year indicates # of years since 1900, e.g., 0 for 1900-MM-DD
        # as.POSIXlt()$mon indicates # of months since 01, e.g., 0 for Jan 
        mons = ncvar_get(dat, 'time')   # months since yr-01-01 00:00
        dyr = mons %/% 12; dmon = mons %% 12  

        mon_0 = 1
        yr_0 = as.numeric(substr(tmn, 1, 4))
        yr_n = formatC(yr_0 + dyr, flag = 0, width = 4)
        mon_n = formatC(mon_0 + dmon, flag = 0, width = 2) 

        times = as.Date(paste0(yr_n, mon_n, '01'), 'UTC', format = '%Y%m%d')
        time.indx = seq(findInterval(as.Date(tmn), times), 
                        findInterval(as.Date(tmx), times))
    }   # end if

    lat.indx  = which(lat >= ymn & lat <= ymx)
    lon.indx  = which(lon >= xmn & lon <= xmx)

    # 4D array with an extra dim of altitude for concentration fields
    if ( levelTF ) {  
        pres = ncvar_get(dat, 'lev')  # mb for 27 pressure levels

        # mol/cc for OH, ppt for NO2 concentration
        var = ncvar_get(dat, tcr.varname, 
                        start = c(lon.indx[1], lat.indx[1], 1, time.indx[1]), 
                        count = c(length(lon.indx), length(lat.indx), 
                                  length(pres), length(time.indx)))
        
        # if 4D array drops as 3D or even matrix, 
        # need to select dimname to match dim of array, DW, 10/04/2020
        if ( tres == 'hourly' ) {
            dim.list = list(lon[lon.indx], lat[lat.indx], pres, mins[time.indx])
            tcol = 'mins'
        } else if ( tres == 'monthly' ) {
            dim.list = list(lon[lon.indx], lat[lat.indx], pres, mons[time.indx])
            tcol = 'mons'
        } 

        ndim = as.numeric(lapply(dim.list, length))
        dim.list = dim.list[which(ndim != 1)]  # select dimnames with ndim > 1
        dimnames(var) = dim.list

        # convert to data frame, e.g., OH concentration in mol cm-3
        colnm = c('lon', 'lat', 'pres', tcol)[which(ndim != 1)]
        
    } else {    # no altitude dim for emission fields

        var = ncvar_get(dat, tcr.varname, 
                        start = c(lon.indx[1], lat.indx[1], time.indx[1]), 
                        count = c(length(lon.indx), length(lat.indx), 
                                  length(time.indx)))

        # convert from kgN/m2/s to umol/m2/s - N or NOx
        if ( grepl('emi', tcr.fn) ) var = var * 1E3 / 14 * 1E6

        # if only spatial matrix 
        dim.list = list(lon[lon.indx], lat[lat.indx], 
                        as.character(times[time.indx]))

        ndim     = as.numeric(lapply(dim.list, length))
        dim.list = dim.list[which(ndim != 1)]  # select dimnames with ndim > 1
        dimnames(var) = dim.list

        # convert to data frame, e.g., OH concentration in mol cm-3
        colnm = c('lon', 'lat', 'time')[which(ndim != 1)]
    }  # end if

    df = melt(var); colnames(df) = c(colnm, tcr.varname)
    nc_close(dat)
    
    if (!'lon' %in% colnm) df$lon = lon[lon.indx] 
    if (!'lat' %in% colnm) df$lat = lat[lat.indx] 
    if (tres == 'hourly' & !'mins' %in% colnm) df$mins = mins[time.indx] 
    if ('mins' %in% colnm) df$time = as.POSIXct(df$mins * 60, 'UTC', 
                                            origin = paste0(yr, '-01-01 01:00'))
    if (!'time' %in% colnames(df)) df$time = times[time.indx]
    
    # force increasing trend of lats and lons
    df = df %>% arrange(lat, lon, time) 
    df
}
