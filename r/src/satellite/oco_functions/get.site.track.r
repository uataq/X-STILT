
# ---------------------------------------------------------------------------- #
#' separate script that includes all site information, including timestr, lat, lon \
#' @author: Dien Wu, 06/15/2018
# ---------------------------------------------------------------------------- #
#' @param lon.lat generate from get.lon.lat()
#' @param dlon @param dlat for find.overpass(), e.g., city.lat +/- dlat
#' @param site name of site, e.g., "Riyadh"
#' @param thred.count.per.deg sounding thredshold, # of soundings per 1 degree lat
# update:
# add season selction, rmTF, remove overpasses during growing seasons, DW, 01/25/2019
# update function for loading OCO-3 data, DW, 06/28/2020 
# combine ggmap.obs.info() if plotTF == TRUE, DW, 10/26/2020
# ---------------------------------------------------------------------------- #

if (F) {
  oco.sensor = obs_sensor
  oco.ver = obs_ver
  oco.path = obs_path
  date.range = c('20140101', '20211231')
  thred.count.per.deg = 100
  thred.count.per.deg.urban = 50
  urban_dlon = urban_dlat = 0.3
  store.path = store_path
  sif.path = sif_path
}


get.site.track = function(site, oco.sensor, oco.ver, oco.path, searchTF = FALSE,
                          date.range = NULL, thred.count.per.deg = 100, lon.lat, 
                          urbanTF = FALSE, urban_dlon = NULL, urban_dlat = NULL,
                          thred.count.per.deg.urban = NULL, rmTF = FALSE, 
                          plotTF = FALSE, store.path = NULL, sif.path = NULL, 
                          qfTF = T){

    library(dplyr)

    if (is.null(date.range)) 
      date.range = c('20140101', format( round(Sys.time(), 'year'), format = '%Y%m%d'))

    # --------------------------------------------------------------------------
    # instead of inputting all info mannually, use find.overpass() to find
    # overpasses given a city lat.lon
    # once you have coordinate info -> get OCO-2 overpasses,
    # either from txt file or scanning through all OCO-2 files
    # first need timestr for SIF files, look for txt.fn from find.overpass()
    
    # txt.path, path for storing output from 'get.site.track()'
    txt.path = file.path(dirname(oco.path), 'overpass_city') 
    dir.create(txt.path, showWarnings = F, recursive = T)
    txt.fn = file.path(txt.path, paste0(oco.sensor, '_overpass_', site, '_', 
                                        oco.ver, '.txt'))

    # if not call find.overpass(); if exists, read from txt.fn
    if (!file.exists(txt.fn) | searchTF == T) {
      cat('NO overpass txt file found or need overwrite...searching now...\n')

      # find overpasses over all OCO-2 time period
      oco.track = find.overpass(date.range, lon.lat, oco.ver, oco.path, 
                                urbanTF, urban_dlon, urban_dlat)
      
      # write output in a txt file
      write.table(oco.track, file = txt.fn, sep = ',', row.names = F, quote = F)
    } # end if


    # --------------------------------------------------------------------------
    oco.track = read.table(txt.fn, header = T, sep = ',', stringsAsFactors = F)

    # select time range and remove tracks with too few soundings
    thred.count = thred.count.per.deg * abs(diff(c(lon.lat$minlat, lon.lat$maxlat)))
    oco.track = oco.track %>% filter(timestr >= date.range[1] &
                                     timestr <= date.range[2] & 
                                     tot.count >= thred.count)

    # at least one sounding near the city
    if (urbanTF) {
      thred.count.urban = thred.count.per.deg.urban * urban_dlat * 2
      oco.track = oco.track %>% filter(tot.urban.count > thred.count.urban)
    } # end if urbanTF

    # track selection, whether to remove summtertime tracks
    if (rmTF) {

      if (lon.lat$citylat > 0) {    # northern Hemi
        oco.track = oco.track %>% filter(substr(timestr, 5, 6) < '05' | 
                                         substr(timestr, 5, 6) > '08')
      } else oco.track = oco.track %>% filter(substr(timestr, 5, 6) < '12' & 
                                                substr(timestr, 5, 6) > '03')

    }  # end if rmTF


    # whether to plot them on maps, plotTF = T/F,
    # this helps you choose which overpass to simulate, see 'tt' below
    if (nrow(oco.track) > 0) {
      ggmap.obs.info(plotTF, site, store.path, all.timestr = oco.track$timestr, 
                     oco.sensor, oco.ver, oco.path, sif.path, lon.lat, 
                     urban_dlat, urban_dlon, qfTF = qfTF)
                    
    } else cat('get.site.track(): NO OCO track available...\n')

    return(oco.track)
} 
# end of script







# --------------------------------------------------------------------------
# add an option of searching # of TROPOMI soundings and get overpass hour 
# on the same day of OCO-2/3 overpasses, DW, 10/26/2020
if (F) {

  cat('get.site.track(): searching for TROPOMI overpasses...this takes a while\n')
  add.track = NULL             # will return TROPOMI overpass info as well
  
  for (t in 1 : nrow(oco.track)) {

    oco.timestr = substr(oco.track$timestr[t], 1, 8)
    tmp.path = ifelse(oco.timestr >= 20190806, tropomi.hir.path[1], tropomi.path[1]) 
    tropomi.info = find.tropomi(tropomi.path = tmp.path, 
                                  timestr = oco.timestr, 
                                  lon.lat = lon.lat)

    if (!is.null(tropomi.info)) 
      add.track = rbind(add.track, cbind(oco.track[t, ], tropomi.info))
  } # end for t

  oco.track = add.track  
} # end if for TROPOMI searching
