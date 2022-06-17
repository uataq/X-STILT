
# ---------------------------------------------------------------------------- #
#' separate script that includes all site information, including timestr, lat, lon \
#' @author: Dien Wu, 06/15/2018
# ---------------------------------------------------------------------------- #
#' @param lon.lat generate from get.lon.lat()
#' @param dlon @param dlat for find.overpass(), e.g., city.lat +/- dlat
#' @param site name of site, e.g., "Riyadh"
#' @param thred.count.per.deg sounding thredshold, # of soundings per 1 degree lat
# update:
# update function for loading OCO-3 data, DW, 06/28/2020 
# combine ggmap.obs.info() if plotTF == TRUE, DW, 10/26/2020
# get rid of any hard-coding about "urban", DW, 06/16/2022
# ---------------------------------------------------------------------------- #

if (F) {
  oco.sensor = obs_sensor
  oco.ver = obs_ver
  oco.path = obs_path
  date.range = c('20140101', '20211231')
  thred.count.per.deg = 100
  thred.count.per.deg.nf = 50
  nf.dlon = nf.dlat = 0.3
  store.path = store_path
  sif.path = sif_path
  lon.lat = lon_lat
}


# ---------------------------------------------------------------------------- #
get.site.track = function(site, oco.sensor, oco.ver, oco.path, searchTF = F,
                          date.range = NULL, thred.count.per.deg = 100, 
                          lon.lat, nfTF = F, nf.dlon = NULL, nf.dlat = NULL, 
                          thred.count.per.deg.nf = NULL, plotTF = F, 
                          store.path = NULL, sif.path = NULL, qfTF = T) {

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
      oco.track = find.oco.overpass(date.range, lon.lat, oco.ver, oco.path, 
                                    nfTF, nf.dlon, nf.dlat)
      
      # write output in a txt file
      write.table(oco.track, file = txt.fn, sep = ',', row.names = F, quote = F)
    } # end if


    # --------------------------------------------------------------------------
    oco.track = read.table(txt.fn, header = T, sep = ',', stringsAsFactors = F)

    # select time range and remove tracks with too few soundings
    thred.count = thred.count.per.deg * abs(diff(c(lon.lat$minlat, 
                                                   lon.lat$maxlat)))
    oco.track = oco.track %>% filter(timestr >= date.range[1] &
                                     timestr <= date.range[2] & 
                                     tot.count >= thred.count)

    # at least one sounding near the city
    if (nfTF) {
      thred.count.nf = thred.count.per.deg.nf * nf.dlat * 2
      oco.track = oco.track %>% filter(tot.nf.count > thred.count.nf)
    } # end if nfTF

    if (qfTF) {
      oco.track = oco.track %>% filter(qf.count >= 5)
      if (nfTF) oco.track = oco.track %>% filter(qf.nf.count >= 5)
    }

    # whether to plot them on maps, plotTF = T/F,
    # this helps you choose which overpass to simulate, see 'tt' below
    if (nrow(oco.track) > 0) {
      ggmap.obs.info(plotTF, site, store.path, all.timestr = oco.track$timestr, 
                     oco.sensor, oco.ver, oco.path, sif.path, lon.lat, 
                     nf.dlat, nf.dlon, qfTF = qfTF)
                    
    } else cat('get.site.track(): NO OCO track available...\n')

    return(oco.track)
} 
# end of script
