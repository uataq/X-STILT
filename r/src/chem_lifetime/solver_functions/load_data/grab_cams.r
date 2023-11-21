# function to read OH and NO2 fields from Chemical reanalysis product from JPL, 
# DW, 09/28/2020 

grab_cams <- function(cams.fn, xmn, xmx, ymn, ymx, tmn, tmx, 
                      cams.speci = c('oh', 'co', 'no', 'no2', 'o3', 'hno3')[1], 
                      xstilt_wd) {

    # load data
    dat <- nc_open(cams.fn)
    lat <- ncvar_get(dat, 'latitude')
    lon <- ncvar_get(dat, 'longitude')       # degree east, [0, 360]
    lon[lon > 180] <- lon[lon > 180] - 360   # convert to [-180, 180]
    nhrs <- ncvar_get(dat, 'time')           # hours since 1900-01-01 00:00:00.0
    time <- as.POSIXct(nhrs * 3600, 'UTC', origin = '1900-01-01 00:00:00')
    
    # L60 model levels, convert levels to pressure
    lev <- ncvar_get(dat, 'level')     

    # load table for coefficients of 60 sigma/hybrid levels for CAMS/ECWMF
    ab_df <- read.csv(file.path(xstilt_wd, 'r/src/chem_lifetime/CAMS_L60.csv'), 
                      stringsAsFactor = F) %>% filter(n %in% lev) %>%
             mutate(pf_hPa = as.numeric(pf_hPa))
    pmid <- ab_df$pf_hPa    # midpoint pressure, 'full-level'

    # locate dim index to subset array 
    lat.indx  <- which(lat >= ymn & lat <= ymx)
    lon.indx  <- which(lon >= xmn & lon <= xmx)
    time.indx <- which(time >= tmn & time <= tmx)
    #print(time[time.indx]); print(lat[lat.indx]); print(lon[lon.indx])

    # kg / kg [lon, lat, level, time]
    oh <- ncvar_get(dat, cams.speci, 
                    start = c(lon.indx[1], lat.indx[1], 1, time.indx[1]), 
                    count = c(length(lon.indx), length(lat.indx), 
                              length(pmid), length(time.indx)))
    
    #t <- ncvar_get(dat, 't', start = c(lon.indx[1], lat.indx[1], 1, time.indx[1]), 
    #                         count = c(length(lon.indx), length(lat.indx), 
    #                                   length(pmid), length(time.indx)))

    # if 4D array drops as 3D or even matrix, 
    # need to select dimname to match dim of array, DW, 10/04/2020
    dim.list <- list(lon[lon.indx], lat[lat.indx], pmid, nhrs[time.indx])
    ndim <- as.numeric(lapply(dim.list, length))
    dim.list <- dim.list[which(ndim != 1)]  # select dimnames with ndim > 1
    dimnames(oh) <- dim.list
  
    # convert to data frame, OH concentration in mol cm-3
    df <- melt(oh) 
    colnm <- c('lon', 'lat', 'pmid', 'nhrs')[which(ndim != 1)]
    colnames(df) <- c(colnm, 'oh')
    
    # convert to data frame
    df <- melt(oh); colnames(df) <- c(colnm, cams.speci)

    if (!'lon' %in% colnm) df$lon <- lon[lon.indx] 
    if (!'lat' %in% colnm) df$lat <- lat[lat.indx] 
    if (!'nhrs' %in% colnm) df$nhrs <- nhrs[time.indx] 
    
    # force increasing trend of lats and lons
    df <- df %>% mutate(time = as.POSIXct(nhrs * 3600, 'UTC', origin = '1900-01-01 00:00:00')) %>% 
                 arrange(lat, lon, time) %>% 
                 left_join(ab_df[, c('ph_hPa', 'pf_hPa')], by = c('pmid' = 'pf_hPa')) %>% 
                 rename(lower.pres = ph_hPa)

    df
}
