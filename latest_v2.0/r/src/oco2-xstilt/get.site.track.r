# separate script that includes all site information, including timestr, lat, lon \
# by DW, 06/15/2018

# 'urbanTF, dlon, dlat' for find.overpass(), e.g., city.lat +/- dlat

get.site.track <- function(site, oco2.ver, oco2.path, searchTF = F,
  date.range = c('20140101', '20181231'), thred.count.per.deg = 100,
  lon.lat, urbanTF, dlon.urban = NULL, dlat.urban = NULL,
  thred.count.per.deg.urban = NULL, txtpath){

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
    oco2.track <- find.overpass(date = c('20140901', '20181231'),
      lon.lat = lon.lat, oco2.ver = oco2.ver, oco2.path = oco2.path,
      urbanTF, dlon.urban, dlat.urban)

    write.table(oco2.track, file = file.path(txtpath, txtfile), sep = ",",
      row.names = F, quote = F)

  } else {   # txtfile found--
    oco2.track <- read.table(file.path(txtpath, txtfile), header = T, sep = ',',
      stringsAsFactors = F)
  } # end if !file.exists()

  # select time range and remove tracks with too few soundings
  thred.count <- thred.count.per.deg * abs(diff(c(lon.lat$minlat, lon.lat$maxlat)))
  cat('Only return overpass dates that have >', thred.count, 'sounding...\n')

  oco2.track <- oco2.track %>% filter(timestr >= date.range[1] &
    timestr <= date.range[2] & tot.count >= thred.count)

  # at least one sounding near the city
  if (urbanTF) {
    thred.count.urban <- thred.count.per.deg.urban * dlat.urban * 2
    oco2.track <- oco2.track %>%
      filter(tot.urban.count > thred.count.urban)
  }

  return(oco2.track)
}

# end of script
