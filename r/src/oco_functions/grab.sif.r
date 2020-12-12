
# ---------------------------------------------------------------------------- #
# subroutine to read and select SIF data, DW 08/07/2018
# ---------------------------------------------------------------------------- #
# sif.ver same as oco.ver, e.g., V7r, V8r, V10r, VEarlyr
# ---------------------------------------------------------------------------- #
grab.sif <- function(sif.path, timestr, lon.lat, sif.ver) {

    # get Sif file name
    sif.file <- list.files(pattern = paste0('LtSIF_', substr(timestr, 3, 8)),
                           path = sif.path, full.names = T)
    oco.sensor <- substr(basename(sif.file), 1, 4)

    if (length(sif.file) == 0) {
      warnings('NO SIF file found for this OCO-2 overpass...'); return()

    } else {  # if file found...

      library(ncdf4)
      sif.dat <- nc_open(sif.file)
      
      if (oco.sensor == 'oco2' & sif.ver != 'V10r') {
        lat     <- ncvar_get(sif.dat, 'latitude')
        lon     <- ncvar_get(sif.dat, 'longitude')
        sif757  <- ncvar_get(sif.dat, 'SIF_757nm') # unit in W/m2/sr/µm
        sif771  <- ncvar_get(sif.dat, 'SIF_771nm')
        igbp    <- ncvar_get(sif.dat, 'IGBP_index')

        sif <- data.frame(timestr = as.numeric(timestr), lat = as.numeric(lat),
                          lon = as.numeric(lon), sif757 = as.numeric(sif757),
                          sif771 = as.numeric(sif771), igbp = as.numeric(igbp))

      } else {
        
        #"seconds since 1993-01-01 00:00:00 UTC"
        sec     <- ncvar_get(sif.dat, 'Geolocation/time_tai93')
        lat     <- ncvar_get(sif.dat, 'Latitude')
        lon     <- ncvar_get(sif.dat, 'Longitude')
        sif757  <- ncvar_get(sif.dat, 'Science/SIF_757nm') # unit in W/m2/sr/µm
        sif771  <- ncvar_get(sif.dat, 'Science/SIF_771nm')
        igbp    <- ncvar_get(sif.dat, 'Science/IGBP_index')

        sif <- data.frame(timestr = as.numeric(timestr), lat = as.numeric(lat),
                          lon = as.numeric(lon), sif757 = as.numeric(sif757),
                          sif771 = as.numeric(sif771), igbp = as.numeric(igbp))
      } 
     
      # select SIF in given region and scale SIF_771nm with sacling factor of
      # 1.35 (used in Luus et al., 2017) to calculate an averaged SIF
      sel.sif <- sif %>% filter(lon >= lon.lat$minlon & lon <= lon.lat$maxlon &
                                lat >= lon.lat$minlat & lat <= lon.lat$maxlat) %>%
                         mutate(avg.sif = (sif757 + sif771 * 1.35)/2)

      # assign months and seasons, only for NH for now, DW
      sel.sif <- sel.sif %>% 
                 mutate(mon = substr(timestr, 5, 6),
                        season = ifelse(mon %in% c('12', '01', '02'), 'WINTER',
                                 ifelse(mon %in% c('03', '04', '05'), 'SPRING',
                                 ifelse(mon %in% c('06', '07', '08'), 'SUMMER', 
                                                                      'FALL'))))
      nc_close(sif.dat)
      return(sel.sif)
    } # end if

} # end if subroutine

