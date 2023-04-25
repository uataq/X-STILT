
# tccon species: 'xh2o', 'xhdo', 'xco', 'xn2o', 'xhf', 'xch4', 'xlco2', 'xwco2', 'xco2', 

# tccon aks: 'ak_xh2o', 'ak_xhdo', 'ak_xco', 'ak_xn2o', 'ak_xo2', 'ak_xhf', 'ak_xch4', 'ak_xlco2', 'ak_xwco2', 'ak_xco2', 'ak_pressure'

grab_tccon = function(tccon_fn, timestr = NULL) {

    library(ncdf4)
    dat = nc_open(tccon_fn)
    lat = ncvar_get(dat, 'lat')
    lon = ncvar_get(dat, 'long')
    #alt = ncvar_get(dat, 'zobs')    # geometric altitude in km
    alt = ncvar_get(dat, 'zmin')     # pressure altitdude 
    sza = ncvar_get(dat, 'solzen')
    saa = ncvar_get(dat, 'azim')     # solar azimuth angles from true north

    # load time in UTC
    # sec_prior = ncvar_get(dat, 'prior_time')   # UTC time for prior GEOS5, seconds since 1970-01-01 00:00:00
    # time_prior = as.POSIXct(sec_prior, tz = 'UTC', origin = '1970-01-01 00:00:00')
    
    yr = ncvar_get(dat, 'year')
    doy = ncvar_get(dat, 'day')     # days of the year 1 - 366
    hr_dec = ncvar_get(dat, 'hour')     # decimal_hour
    date = as.Date(doy - 1, origin = paste0(yr, '-01-01'))
    datestr = as.POSIXct(hr_dec * 3600, 'UTC', date)
    all_timestr = format(datestr, tz = 'UTC', format = '%Y%m%d%H')

    # wind obs
    ws = ncvar_get(dat, 'wspd')     # wind speed in m s-1
    wd = ncvar_get(dat, 'wdir')     # wind dir from 

    # column average dry mole fraction, 0.2095 * column CO / column O2
    xco_ppb = ncvar_get(dat, 'xco')     # ppb
    xco_err = ncvar_get(dat, 'xco_error')

    xch4_ppm = ncvar_get(dat, 'xch4')   # WMO CH4 X2004, ppm 
    xch4_err = ncvar_get(dat, 'xch4_error')

    # 0.2095*column_co2/column_o2 lco2 is the strong CO2 band centered at 4852.87 cm-1 and does not contribute to the xco2 calculation. These data are EXPERIMENTAL. If you plan to use them, please work with the site PI.
    #xco2_strong_ppm = ncvar_get(dat, 'xlco2_experimental')  # WMO CO2 X2007
    #xco2_weak_ppm = ncvar_get(dat, 'xwco2_experimental')
    
    xco2_2007_ppm = ncvar_get(dat, 'xco2')                       
    xco2_2019_ppm = ncvar_get(dat, 'xco2_x2019')            # WMO CO2 X2019
    xco2_2019_err = ncvar_get(dat, 'xco2_error_x2019')

    nc_close(dat)

    tccon_df = data.frame(lon, lat, alt, sza, saa, yr, datestr, 
                          timestr = all_timestr, ws, wd, xco_ppb, xco_err, 
                          xch4_ppm, xch4_err, xco2_2019_ppm, xco2_2019_err)

    if ( !is.null(timestr) ) tccon_df = tccon_df[tccon_df$timestr %in% timestr,]

    return(tccon_df)
}


