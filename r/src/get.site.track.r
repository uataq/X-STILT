#' separate script that includes all site information, including timestr, lat, lon \
#' @author: Dien Wu, 06/15/2018

#' @param:
#' urbanTF, dlon, dlat: for find.overpass(), e.g., city.lat +/- dlat
#' site: name of site, e.g., "Riyadh"
#' thred.count.per.deg: sounding thredshold, # of soundings per 1 degree lat
#' txtpath: path to store an output txtfile
#' lon.lat: generate from get.lon.lat()

#' @updates:
#' add season selction, rmTF, remove overpasses during growing seasons, DW, 01/25/2019

get.site.track <- function(site, oco2.ver, oco2.path, searchTF = F,
                           date.range = c('20140101', '20181231'), 
                           thred.count.per.deg = 100, lon.lat, urbanTF, 
                           dlon.urban = NULL, dlat.urban = NULL,
                           thred.count.per.deg.urban = NULL, txtpath, 
                           rmTF = F){

  library(dplyr)

  # instead of inputting all info mannually, use find.overpass() to find
  # overpasses given a city lat.lon
  # once have coordinate info, get OCO-2 overpasses,
  # either from txt file or scanning through all OCO-2 files
  # first need timestr for SIF files, look for txtfile from find.overpass()
  txtfile <- paste0('oco2_overpass_', site, '_', oco2.ver, '.txt')

  # if not call find.overpass(); if exists, read from txtfile
  if (!file.exists(file.path(txtpath, txtfile)) | searchTF == T) {
    cat('NO overpass txt file found or need overwrite...searching now...\n')

    # find overpasses over all OCO-2 time period
    oco2.track <- find.overpass(date.range, lon.lat, oco2.ver, oco2.path, 
                                urbanTF, dlon.urban, dlat.urban)

    write.table(oco2.track, file = file.path(txtpath, txtfile), sep = ',',
                row.names = F, quote = F)

  } else {   # txtfile found--
    oco2.track <- read.table(file.path(txtpath, txtfile), header = T, sep = ',',
                             stringsAsFactors = F)
  } # end if !file.exists()

  # select time range and remove tracks with too few soundings
  thred.count <- thred.count.per.deg * abs(diff(c(lon.lat$minlat, lon.lat$maxlat)))
  cat('Only return overpass dates that have >', thred.count, 'sounding...\n')
  oco2.track <- oco2.track %>% filter(timestr >= date.range[1] &
                                      timestr <= date.range[2] & 
                                      tot.count >= thred.count)

  # at least one sounding near the city
  if (urbanTF) {
    thred.count.urban <- thred.count.per.deg.urban * dlat.urban * 2
    oco2.track <- oco2.track %>% filter(tot.urban.count > thred.count.urban)
  } # end if urbanTF

  # track selection, whether to remove summtertime tracks
  if (rmTF) {
    if (lon.lat$citylat > 0) {   
      # northern Hemi
      oco2.track <- oco2.track %>% filter(substr(timestr, 5, 6) < '05' | 
                                          substr(timestr, 5, 6) > '08')
    } else {  

      # southern Hemi
      oco2.track <- oco2.track %>% filter(substr(timestr, 5, 6) < '12' & 
                                          substr(timestr, 5, 6) > '03')
    } # end if Hemisphere
  }  # end if rmTF

  return(oco2.track)
}

# end of script
