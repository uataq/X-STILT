#' separate script that includes all site information, including timestr, lat, lon \
#' @author: Dien Wu, 06/15/2018

#' @param lon.lat generate from get.lon.lat()
#' @param dlon @param dlat for find.overpass(), e.g., city.lat +/- dlat
#' @param site name of site, e.g., "Riyadh"
#' @param thred.count.per.deg sounding thredshold, # of soundings per 1 degree lat

# ------------------------------ Update ------------------------------------- #
# add season selction, rmTF, remove overpasses during growing seasons, DW, 01/25/2019
# update function for loading OCO-3 data, DW, 06/28/2020 

get.site.track <- function(site, oco.sensor, oco.ver, oco.path, searchTF = F,
                           date.range = c('20140101', '20201231'), 
                           thred.count.per.deg = 100, lon.lat, urbanTF, 
                           dlon.urban = NULL, dlat.urban = NULL,
                           thred.count.per.deg.urban = NULL, rmTF = F){

    library(dplyr)

    # instead of inputting all info mannually, use find.overpass() to find
    # overpasses given a city lat.lon
    # once you have coordinate info -> get OCO-2 overpasses,
    # either from txt file or scanning through all OCO-2 files
    # first need timestr for SIF files, look for txt.file from find.overpass()
    
    # txt.path, path for storing output from 'get.site.track()'
    txt.path <- file.path(dirname(oco.path), 'overpass_city') 
    dir.create(txt.path, showWarnings = F, recursive = T)
    txt.file <- file.path(txt.path, paste0(oco.sensor, '_overpass_', site, '_', 
                                           oco.ver, '.txt'))

    # if not call find.overpass(); if exists, read from txt.file
    if (!file.exists(txt.file) | searchTF == T) {
      cat('NO overpass txt file found or need overwrite...searching now...\n')

      # find overpasses over all OCO-2 time period
      oco.track <- find.overpass(date.range, lon.lat, oco.ver, oco.path, 
                                 urbanTF, dlon.urban, dlat.urban)
      write.table(oco.track, file = txt.file, sep = ',', row.names = F, quote = F)

    } else {   # txt.file found--
      oco.track <- read.table(txt.file, header = T, sep = ',', stringsAsFactors = F)

    } # end if !file.exists()

    # select time range and remove tracks with too few soundings
    thred.count <- thred.count.per.deg * abs(diff(c(lon.lat$minlat, lon.lat$maxlat)))
    cat('Only return overpass dates that have >', thred.count, 'sounding...\n')
    oco.track <- oco.track %>% filter(timestr >= date.range[1] &
                                      timestr <= date.range[2] & 
                                      tot.count >= thred.count)

    # at least one sounding near the city
    if (urbanTF) {
      thred.count.urban <- thred.count.per.deg.urban * dlat.urban * 2
      oco.track <- oco.track %>% filter(tot.urban.count > thred.count.urban)
    } # end if urbanTF

    # track selection, whether to remove summtertime tracks
    if (rmTF) {
      if (lon.lat$citylat > 0) {   
        # northern Hemi
        oco.track <- oco.track %>% filter(substr(timestr, 5, 6) < '05' | 
                                          substr(timestr, 5, 6) > '08')
      } else {  

        # southern Hemi
        oco.track <- oco.track %>% filter(substr(timestr, 5, 6) < '12' & 
                                          substr(timestr, 5, 6) > '03')
      } # end if Hemisphere
    }  # end if rmTF

    return(oco.track)
}

# end of script
