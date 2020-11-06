# function to read OH and NO2 fields from Chemical reanalysis product from JPL, 
# DW, 09/28/2020 

grab_tcr <- function(tcr.fn, xmn, xmx, ymn, ymx, tmn, tmx) {

    # load data
    dat <- nc_open(tcr.fn)
    lat <- ncvar_get(dat, 'lat')                # uneven latitudes
    lon <- ncvar_get(dat, 'lon')                # degree east, [0, 360]
    lon[lon > 180] <- lon[lon > 180] - 360      # convert to [-180, 180]
    
    pres <- ncvar_get(dat, 'lev')               # mb for 27 pressure levels
    mins <- ncvar_get(dat, 'time')              # minutes since 2019-01-01 01:00
    yr   <- substr(tmn, 1, 4)[1]
    time <- as.POSIXct(mins * 60, 'UTC', origin = paste0(yr, '-01-01 01:00'))
    
    lat.indx  <- which(lat >= ymn & lat <= ymx)
    lon.indx  <- which(lon >= xmn & lon <= xmx)
    time.indx <- which(time >= tmn & time <= tmx)

    # mol/cc for OH, and ppt for NO2
    oh <- ncvar_get(dat, 'oh', start = c(lon.indx[1], lat.indx[1], 1, time.indx[1]), 
                               count = c(length(lon.indx), length(lat.indx), 
                                         length(pres), length(time.indx)))
    
    # if 4D array drops as 3D or even matrix, 
    # need to select dimname to match dim of array, DW, 10/04/2020
    dim.list <- list(lon[lon.indx], lat[lat.indx], pres, mins[time.indx])
    ndim     <- as.numeric(lapply(dim.list, length))
    dim.list <- dim.list[which(ndim != 1)]  # select dimnames with ndim > 1
    dimnames(oh) <- dim.list
  
    # convert to data frame, OH concentration in mol cm-3
    df <- melt(oh); colnm <- c('lon', 'lat', 'pres', 'mins')[which(ndim != 1)]
    colnames(df) <- c(colnm, 'oh')
    
    if (!'lon' %in% colnm) df$lon <- lon[lon.indx] 
    if (!'lat' %in% colnm) df$lat <- lat[lat.indx] 
    if (!'mins' %in% colnm) df$mins <- mins[time.indx] 
    
    df$time <- as.POSIXct(df$mins * 60, 'UTC', origin = paste0(yr, '-01-01 01:00'))
    df <- df %>% arrange(lat, lon, time)  # force increasing trend of lats and lons

    df
}
