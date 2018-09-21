# subroutine to write 3D array of biospheric or oceanic contributions to nc files
# originated from 'write_footprint.r' by Ben Fasoli
# DW, 09/15/2018 

write_contri <- function(contribution, time_out, lon, lat, 
    store.file, str = c('bio', 'ocn')) {

    # write 3D array in nc file
    xdim <- ncdim_def('lon', 'degrees_east', lon)
    ydim <- ncdim_def('lat', 'degrees_north', lat)
    tdim <- ncdim_def('time', 'seconds since 1970-01-01 00:00:00Z',
                      as.numeric(time_out))
    fvar <- ncvar_def(paste0('xco2.', str), 'ppm', 
                      list(xdim, ydim, tdim), -1)
                      
    system(paste('rm -f ', store.file))
    nc <- nc_create(store.file, list(fvar), force_v4 = T)
    ncatt_put(nc, 'lon', 'standard_name', 'longitude')
    ncatt_put(nc, 'lon', 'long_name', 'longitude at cell center')
    ncatt_put(nc, 'lat', 'standard_name', 'latitude')
    ncatt_put(nc, 'lat', 'long_name', 'latitude at cell center')

    # Insert contribution data
    ncvar_put(nc, fvar, contribution)

    ncatt_put(nc, 'time', 'standard_name', 'time')
    ncatt_put(nc, 'time', 'long_name', 'utc time')
    ncatt_put(nc, 'time', 'calendar', 'standard')
    
    ncatt_put(nc, paste0('xco2.', str), 'standard_name', 'contribution')
    ncatt_put(nc, paste0('xco2.', str), 'long_name', 
              paste0('XCO2 contribution due to ', str, ' fluxes'))
    
    ncatt_put(nc, 0, 'crs', '+proj=longlat')
    ncatt_put(nc, 0, 'crs_format', 'PROJ.4')
    ncatt_put(nc, 0, 'documentation', 'github.com/uataq/stilt')
    ncatt_put(nc, 0, 'time_created', format(Sys.time(), tz = 'UTC'))
    return(store.file)
}
