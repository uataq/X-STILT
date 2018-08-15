# subroutine for running X-STILT trajec
# Written by DW, update on 05/19/2017 (add GDAS-driven STILT)
# ADD fixed agl run, 07/07/2017, DW
# clear up, 04/30/2017, DW

# read variables from output of 'create_namelist_trajec.r'
# use variables instead of namelist, DW, 07/25/2018
# for generating trajec with horizontal error compoennt,
# use the same lat.lon from original trajec, DW, 07/31/2018

get.recp.info <- function(timestr, oco2.path, oco2.ver, lon.lat, selTF,
  recp.indx, recp.num, find.lat, agl, plotTF = F, run_trajec, run_hor_err,
  trajpath = NULL, stilt.ver = 2){

  # ------------------- Step 1. READ IN OCO-2 LITE FILES ------------------- #
  source('r/dependencies.r') # source all functions
  oco2 <- grab.oco2(oco2.path, timestr, lon.lat)
  # filter by quality flag too, more receptors when XCO2 is high
  if (oco2.ver == 'b7rb') sel.oco2 <- oco2 %>% filter(qf == 0)
  if (oco2.ver == 'b8r')  sel.oco2 <- oco2 %>% filter(wl <= 1)

  # round lat, lon for each sounding, fix bug, DW, 07/31/2018
  sel.oco2 <- sel.oco2 %>% mutate(lat = signif(lat, 6), lon = signif(lon, 7))

  # whether plotting XCO2 from OCO-2
  if (plotTF) {
    zoom <- 8
    font.size <- rel(1.0)
    col <- def.col()
    col.range <- seq(380, 420, 2)
    m1 <- ggplot.map(map = 'ggmap', center.lat = lon.lat$citylat,
          center.lon = lon.lat$citylon, zoom = zoom)[[1]] + theme_bw()

    # add observed XCO2
    c1 <- m1 + geom_point(data = sel.oco2, aes(lon, lat, colour = xco2))
    c1 <- c1 + labs(x = 'LONGITUDE [degE]', y = 'LATITUDE [degN]')
    c1 <- c1 + labs(title = paste('OCO-2 XCO2 [ppm] for', site, 'on',
      substr(timestr, 1, 8)))
    c1 <- c1 + scale_colour_gradientn(name = 'OCO-2 XCO2 [ppm]',
      colours = col, breaks = col.range, labels = col.range)

    # add themes
    c2 <- c1 + theme(legend.position = 'bottom',
      legend.text = element_text(size = font.size),
      legend.key = element_blank(), legend.key.height = unit(0.5, 'cm'),
      legend.key.width = unit(3, 'cm'),
      axis.title.y = element_text(size = font.size, angle = 90),
      axis.title.x = element_text(size = font.size, angle = 0),
      axis.text = element_text(size = font.size),
      axis.ticks = element_line(size = font.size),
      title = element_text(size = font.size))
    picname <- paste0('ggmap_xco2_', site,'_', substr(timestr, 1, 8), '.png')
    ggsave(c2, filename = picname, width = 6, height = 6)
  }  # end if plotTF

  # ------------------- Step 2. SET UP the STILT receptors ----------------- #
  # if generate trajec with horizontal error component,
  # use the same lat.lon from original trajec, DW, 07/31/2018
  if (run_trajec & run_hor_err) {

    trajfile  <- list.files(path = trajpath, pattern = 'X_traj.rds',
      recursive = T)
    trajfile  <- file.path(trajpath, trajfile)
    trajname  <- basename(trajfile)
    recp.info <- ident.to.info(ident = trajname, stilt.ver = stilt.ver)

    # get overpass time by merging with observed XCO2 and compute run_time
    sel.oco2$find.lat <- signif(sel.oco2$lat, max(nchar(recp.info$recp.lat)) - 1)
    sel.oco2$find.lon <- signif(sel.oco2$lon, max(nchar(recp.info$recp.lon)) - 1)
    recp.info <- recp.info %>%
      dplyr::select('lati' = 'recp.lat', 'long' = 'recp.lon') %>%
      left_join(sel.oco2, by = c('lati' = 'find.lat', 'long' = 'find.lon')) %>%
      mutate(run_times_utc = as.POSIXct(substr(id, 1, 14), '%Y%m%d%H%M%S',
        tz = 'UTC'))

  } else {

    ### if generate trajec without error component
    if(selTF){
      # select lat, lon based on OCO-2 soundings and data quality
      sel.lat  <- sel.oco2$lat
      sort.lat <- sort(sel.lat)
      recp.lat <- sort.lat[findInterval(recp.indx, sort.lat)]
      match.index <- unique(match(recp.lat, sel.lat))
      recp.info <- sel.oco2[match.index, ]

    }else{   # if no requirement, simulate all soundings, no selection
      recp.info <- sel.oco2
    }  # end if selTF

    # compute simulation timing, yyyy-mm-dd HH:MM:SS (UTC), aka 'receptors' info
    # that are used for each simulation and match Ben's code
    recp.info <- recp.info %>%
      mutate(run_times_utc = as.POSIXct(substr(id, 1, 14), '%Y%m%d%H%M%S',
        tz = 'UTC'))
  } # end if run_hor_err

  recp.info <- recp.info[order(recp.info$lat), ]
  recp.info <- recp.info %>%
    dplyr::select('run_time' = 'run_times_utc', 'lati' = 'lat', 'long' = 'lon')

  # subset receptor data frame or find the nearest lat..
  if (!is.null(recp.num)) recp.info <- recp.info[min(recp.num):max(recp.num), ]
  if (!is.null(find.lat)) recp.info <- recp.info[findInterval(find.lat, recp.info$lati), ]

  ## add release height
  recp.info$zagl <- agl

  # return receptor info
  return(recp.info)
} # end of subroutine
