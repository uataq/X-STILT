# function to read OH and NO2 fields from WACCM, DW, 09/28/2020 

grab_waccm <- function(waccm.fn, xmn, xmx, ymn, ymx, tmn, tmx) {

    dat <- nc_open(waccm.fn)
    waccm.timestr <- substr(waccm.fn, nchar(waccm.fn) - 18, nchar(waccm.fn) - 9)

    lat <- ncvar_get(dat, 'lat')
    lon <- ncvar_get(dat, 'lon')        # degree east, [0, 360]
    lon[lon > 180] <- lon[lon > 180] - 360      # convert to [-180, 180]
    lev  <- ncvar_get(dat, 'lev')        # hybrid level at midpoints (1000*(A+B)) in hPa, 88
    ilev <- ncvar_get(dat, 'ilev')      # hybrid level at interfaces (1000*(A+B)) in hPa, 89
    nhrs <- ncvar_get(dat, 'time') * 24     # hours since 00:00:00
    
    P0 <- ncvar_get(dat, 'P0') / 100    # reference press in Pa, convert to hPa now
    PS <- ncvar_get(dat, 'PS') / 100    # sruface press in Pa, convert to hPa now
    hyai <- ncvar_get(dat, 'hyai')      # hybrid A coefficient at layer interfaces
    hybi <- ncvar_get(dat, 'hybi')      # hybrid B coefficient at layer interfaces
    #hrs  <- seq(0, 23, 6)               # start hour of the 6-hour interval 

    lat.indx <- which(lat >= ymn & lat <= ymx)
    lon.indx <- which(lon >= xmn & lon <= xmx)

    # grab OH in mol / mol, [lon, lat, lev, time]
    oh <- ncvar_get(dat, 'OH', start = c(lon.indx[1], lat.indx[1], 1, 1), 
                    count = c(length(lon.indx), length(lat.indx), length(lev), length(nhrs)))
    dimnames(oh) <- list(lon[lon.indx], lat[lat.indx], seq(1, length(lev)), nhrs)
    oh_df <- reshape2::melt(oh); colnames(oh_df) <- c('lon', 'lat', 'l', 'hrs', 'oh')

    # calculate the mean coefficient for middle pressure
    # atmosphere_hybrid_sigma_pressure_coordinate 
    # p(n,k,j,i) = a(k)*p0 + b(k)*ps(n,j,i), from TOA to the surface at interfaces
    ab_df <- data.frame(a = zoo::rollmean(hyai, 2), b = zoo::rollmean(hybi, 2)) %>% 
             mutate(l = seq(1, length(lev)))
    dimnames(PS) <- list(lon, lat, nhrs)
    ps_df <- reshape2::melt(PS); colnames(ps_df) <- c('lon', 'lat', 'hrs', 'ps')

    # convert to data frame 
    df <- oh_df %>% left_join(ps_df, by = c('lon', 'lat', 'hrs')) %>% 
          left_join(ab_df, by = 'l') %>% 
          mutate(pmid = a * P0 + b * ps, 
                 time = as.POSIXct(hrs * 3600, 'UTC', origin = paste0(waccm.timestr, '00:00:00'))) %>% 
          filter(time >= tmn, time <= tmx)
    
    df
    
}
