# subroutine to determine receptor locations based on obs soundings (centered lat/lon)
# DW, 07/02/2021 

# in this version, we always use OCO soundings as receptors
get.recp.info = function(timestr, oco.ver, oco.path, lon.lat, selTF, peak.lat,
                         num.bg, num.peak, recp.num = NULL, find.lat = NULL, agl, 
                         run_trajec, output_wd = NULL, data.filter = c('QF', 0), 
                         tropomi_info){

  # ------------------- Step 1. READ IN OCO-2 LITE FILES ------------------- #
  oco.dat = grab.oco(oco.path, timestr, lon.lat, oco.ver) %>% 
            dplyr::select(-c('lons', 'lats', 'vertices')) %>% unique()
  if (nrow(oco.dat) == 0) cat('No sounding found over this overpass for this region\n')
  
  # filter by quality flag too, more receptors when XCO2 is high
  if (is.null(data.filter)) { # no data filtering, use all OCO data
    sel.dat = oco.dat    
  } else if (data.filter[1] == 'QF') {
    sel.dat = oco.dat %>% filter(qf <= data.filter[2])
  } else if (data.filter[1] == 'WL') {
    sel.dat = oco.dat %>% filter(wl <= data.filter[2])
  } else {
    cat('get.recp.info(): Incorrect @param "data.filter"... please check\n')
  } # end if

  # round lat, lon for each sounding, fix bug, DW, 07/31/2018
  sel.dat = sel.dat %>% mutate(lat = signif(lat, 6), lon = signif(lon, 7))


  # ------------------- Step 2. SET UP the STILT receptors ----------------- #
  # if generate trajec with horizontal error component,
  # use the same lat.lon from original trajec, DW, 07/31/2018
  trajfile = list.files(file.path(output_wd, 'by-id'), 'X_traj.rds', 
                        recursive = T, full.names = T)
  tropomiTF = unique(tropomi_info$tropomiTF)

  if ( length(trajfile) > 0 & !run_trajec ) {

    cat('Found existing trajectories...\n') # if trajec data exists
    trajname  = basename(trajfile)
    recp.info = ident.to.info(ident = trajname, stilt.ver = 2)

    if (tropomiTF) {  # for TROPOMI runs, use its overpass time
      recp.info = recp.info %>% mutate(run_time = as.POSIXct(as.character(timestr), 
                                                            'UTC', format = '%Y%m%d%H%M')) %>% 
                  dplyr::rename(lati = recp.lat, long = recp.lon) %>% arrange(lati)
    
    } else {
      
      # get OCO overpass time by merging with observed XCO2 and compute run_time
      sel.dat$find.lat = signif(sel.dat$lat, max(nchar(recp.info$recp.lat)) - 1)
      sel.dat$find.lon = signif(sel.dat$lon, max(nchar(recp.info$recp.lon)) - 1)
      recp.info = recp.info %>% 
                  dplyr::select('lati' = 'recp.lat', 'long' = 'recp.lon') %>%
                  left_join(sel.dat, by = c('lati' = 'find.lat', 'long' = 'find.lon')) %>%
                  mutate(run_time = as.POSIXct(substr(id, 1, 14), 'UTC', format = '%Y%m%d%H%M%S')) %>% 
                  arrange(lati)
      
    } # end if 


  } else {

    ### if generate trajec without error component
    if (selTF) {
      
      # place denser receptors within lat range with high XCO2
      recp.indx = c(seq(lon.lat$minlat,  peak.lat[1],     1 / num.bg),
                    seq(peak.lat[1],     peak.lat[2],     1 / num.peak),
                    seq(peak.lat[1],     lon.lat$maxlat,  1 / num.bg))

      # select lat, lon based on OCO-2 soundings and data quality
      sel.lat  = sel.dat$lat
      sort.lat = sort(sel.lat)
      recp.lat = sort.lat[findInterval(recp.indx, sort.lat)]
      match.index = unique(match(recp.lat, sel.lat))
      recp.info = sel.dat[match.index, ]

    } else {   # if no requirement, simulate all soundings, no selection
      recp.info = sel.dat
    }  # end if selTF

    # compute simulation timing, yyyy-mm-dd HH:MM:SS (UTC), aka 'receptors' info
    # that are used for each simulation and match Ben's code
    recp.info = recp.info %>% 
                mutate(run_times_utc = as.POSIXct(substr(id, 1, 14), 'UTC', 
                                                  format = '%Y%m%d%H%M%S')) %>% 
                dplyr::select(run_time = run_times_utc, lati = lat, long = lon) %>% arrange(lati)

    # subset receptor data frame or find the nearest lat..
    if (!is.null(recp.num)) recp.info = recp.info[min(recp.num) : max(recp.num), ]
    if (!is.null(find.lat)) recp.info = recp.info[findInterval(find.lat, recp.info$lati), ]

    if (tropomiTF) {  # for TROPOMI runs, use its overpass time

      cat('get.recp.info(): tropomiTF == TRUE, use TROPOMI overpass time as receptor times...\n')
      tinfo = tropomi_info[1, ]
      ttime = seq(as.POSIXct(as.character(tinfo$overpass.start.time), 'UTC', format = '%Y%m%d%H%M%S'), 
                  as.POSIXct(as.character(tinfo$overpass.end.time), 'UTC', format = '%Y%m%d%H%M%S'), 
                  length = nrow(recp.info))
      recp.info$run_time = ttime
    } # end if

  } # end if trajec file existed

  ## add release height
  recp.info$zagl = list(agl)

  # return receptor info
  return(recp.info)
} # end of subroutine
