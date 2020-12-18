# subroutine for running X-STILT trajec
# Written by DW, update on 05/19/2017 (add GDAS-driven STILT)
# ADD fixed agl run, 07/07/2017, DW
# clear up, 04/30/2017, DW

# read variables from output of 'create_namelist_trajec.r'
# use variables instead of namelist, DW, 07/25/2018
# for generating trajec with horizontal error compoennt,
# use the same lat.lon from original trajec, DW, 07/31/2018
# add run_trajec, if run_trajec == T, always recalculate 'recp.info', 
#   rather than grab from exising trajec, DW, 10/30/2018

# add options for using TROPOMI overpass hour, but still keep the finer lat/lon
# soundings that match OCO-2/3, DW, 10/09/2020

get.recp.info <- function(timestr, oco.ver, oco.path, lon.lat, selTF, recp.indx,
                          recp.num, find.lat, agl, run_trajec, outdir = NULL, 
                          data.filter = c('QF', 0), tropomiTF = F, 
                          tropomi.path = NULL, tropomi.fn = NULL){

  # ------------------- Step 1. READ IN OCO-2 LITE FILES ------------------- #
  oco.dat <- grab.oco(oco.path, timestr, lon.lat, oco.ver) %>% 
             dplyr::select(-c('lons', 'lats', 'vertices')) %>% unique()
  if (nrow(oco.dat) == 0) cat('No sounding found over this overpass for this region\n')
  
  # filter by quality flag too, more receptors when XCO2 is high
  if (is.null(data.filter)) { # no data filtering, use all OCO data
    sel.dat <- oco.dat    
  } else if (data.filter[1] == 'QF') {
    sel.dat <- oco.dat %>% filter(qf <= data.filter[2])
  } else if (data.filter[1] == 'WL') {
    sel.dat <- oco.dat %>% filter(wl <= data.filter[2])
  } else {
    cat('get.recp.info(): Incorrect @param "data.filter"... please check\n')
  } # end if

  # round lat, lon for each sounding, fix bug, DW, 07/31/2018
  sel.dat <- sel.dat %>% mutate(lat = signif(lat, 6), lon = signif(lon, 7))


  # ------------------- Step 2. SET UP the STILT receptors ----------------- #
  # if generate trajec with horizontal error component,
  # use the same lat.lon from original trajec, DW, 07/31/2018
  trajfile <- list.files(file.path(outdir, 'by-id'), 'X_traj.rds', 
                         recursive = T, full.names = T)

  if ( length(trajfile) > 0 & !run_trajec ) {

    cat('Found existing trajectories...\n') # if trajec data exists
    trajname  <- basename(trajfile)
    recp.info <- ident.to.info(ident = trajname, stilt.ver = 2)

    if (tropomiTF) {  # for TROPOMI runs, use its overpass time
      recp.info <- recp.info %>% mutate(run_time = as.POSIXct(as.character(timestr), 
                                                              'UTC', '%Y%m%d%H%M')) %>% 
                   dplyr::rename(lati = recp.lat, long = recp.lon) %>% arrange(lati)
    
    } else {
      
      # get OCO overpass time by merging with observed XCO2 and compute run_time
      sel.dat$find.lat <- signif(sel.dat$lat, max(nchar(recp.info$recp.lat)) - 1)
      sel.dat$find.lon <- signif(sel.dat$lon, max(nchar(recp.info$recp.lon)) - 1)
      recp.info <- recp.info %>% 
                  dplyr::select('lati' = 'recp.lat', 'long' = 'recp.lon') %>%
                  left_join(sel.dat, by = c('lati' = 'find.lat', 'long' = 'find.lon')) %>%
                  mutate(run_time = as.POSIXct(substr(id, 1, 14), '%Y%m%d%H%M%S', tz = 'UTC')) %>% 
                  arrange(lati)
      
    } # end if 


  } else {

    ### if generate trajec without error component
    if (selTF) {
      # select lat, lon based on OCO-2 soundings and data quality
      sel.lat  <- sel.dat$lat
      sort.lat <- sort(sel.lat)
      recp.lat <- sort.lat[findInterval(recp.indx, sort.lat)]
      match.index <- unique(match(recp.lat, sel.lat))
      recp.info <- sel.dat[match.index, ]

    } else {   # if no requirement, simulate all soundings, no selection
      recp.info <- sel.dat
    }  # end if selTF

    # compute simulation timing, yyyy-mm-dd HH:MM:SS (UTC), aka 'receptors' info
    # that are used for each simulation and match Ben's code
    recp.info <- recp.info %>% mutate(run_times_utc = as.POSIXct(substr(id, 1, 14), 
                                                                 format = '%Y%m%d%H%M%S',
                                                                 tz = 'UTC')) %>% 
                               dplyr::select('run_time' = 'run_times_utc', 
                                             'lati' = 'lat', 'long' = 'lon') %>% 
                               arrange(lati)


    # subset receptor data frame or find the nearest lat..
    if (!is.null(recp.num)) recp.info <- recp.info[min(recp.num) : max(recp.num), ]
    if (!is.null(find.lat)) recp.info <- recp.info[findInterval(find.lat, recp.info$lati), ]

    if (tropomiTF) {  # for TROPOMI runs, use its overpass time

      cat('get.recp.info(): tropomiTF == TRUE, use TROPOMI instead of OCO overpass time...\n')
      tropomi.dat <- grab.tropomi.co(lon.lat = lon.lat, tropomi.fn = tropomi.fn) %>% na.omit()
      #print(as.character(mean(tropomi.dat$time_utc)))
      
      tropomi.dat$date = as.POSIXct(as.character(tropomi.dat$time_utc), 'UTC', 
                                    format = '%Y%m%d%H%M%S')
      tropomi.time <- seq(min(tropomi.dat$date), max(tropomi.dat$date), 
                          length = nrow(recp.info))
      recp.info$run_time <- tropomi.time
    }

  } # end if trajec file existed

  ## add release height
  recp.info$zagl <- list(agl)

  # return receptor info
  return(recp.info)
} # end of subroutine
