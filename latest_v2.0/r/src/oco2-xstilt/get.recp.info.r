# subroutine for running X-STILT trajec
# Written by DW, update on 05/19/2017 (add GDAS-driven STILT)
# ADD fixed agl run, 07/07/2017, DW
# clear up, 04/30/2017, DW

# read variables from output of 'create_namelist_trajec.r'
# use variables instead of namelist, DW, 07/25/2018

get.recp.info <- function(timestr, oco2.path, lon.lat, selTF, recp.indx,
  recp.num, find.lat, agl, plotTF){

  # ------------------- Step 1. READ IN OCO-2 LITE FILES --------------------- #
  YYYYMMDD  <- substr(timestr, 1, 8)
  YYMMDD    <- substr(timestr, 3, 8)
  oco2.file <- list.files(pattern = YYMMDD, path = oco2.path)
  oco2.dat  <- nc_open(file.path(oco2.path, oco2.file))

  # grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  # may need to reverse later
  oco2.lat    <- ncvar_get(oco2.dat, 'latitude')
  oco2.lon    <- ncvar_get(oco2.dat, 'longitude')
  oco2.foot   <- ncvar_get(oco2.dat, 'Sounding/footprint')
  oco2.level  <- ncvar_get(oco2.dat, 'levels')
  xco2        <- ncvar_get(oco2.dat, 'xco2')
  xco2.uncert <- ncvar_get(oco2.dat, 'xco2_uncertainty')  # posterior error

  # grabbing time for STILT receptors
  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id   <- as.character(ncvar_get(oco2.dat, 'sounding_id'))
  YYYY <- substr(id, 1, 4)
  MM   <- substr(id, 5, 6)
  DD   <- substr(id, 7, 8)
  HH   <- substr(id, 9, 10)
  mm   <- substr(id, 11, 12)
  ss   <- substr(id, 13, 14)
  ms   <- substr(id, 15, 15)
  f    <- substr(id, 16, 16)
  YYYYMMDDHHmmss <- substr(id, 1, 14)

  # 0:Nadir, 1:Glint, 2:Target, 3: Transition
  OM <- ncvar_get(oco2.dat, 'Sounding/operation_mode')
  LF <- ncvar_get(oco2.dat, 'Sounding/land_fraction')   # > 80%: land, < 20%: sea
  QF <- ncvar_get(oco2.dat, 'xco2_quality_flag')
  WL <- ncvar_get(oco2.dat, 'warn_level')

  # operation modes:
  mode <- rep(NULL, length(OM))
  mode[LF > 80 & OM == 0] <- 'Land_Nadir'
  mode[LF < 20 & OM == 0] <- 'Sea_Nadir'
  mode[LF > 80 & OM == 1] <- 'Land_Glint'
  mode[LF < 20 & OM == 1] <- 'Sea_Glint'
  mode[LF > 80 & OM == 2] <- 'Land_Target'
  mode[LF < 20 & OM == 2] <- 'Sea_Target'
  mode[LF > 80 & OM == 3] <- 'Land_Transition'
  mode[LF < 20 & OM == 3] <- 'Sea_Transition'

  ## SELECT REGION,  minlat, maxlat, minlon, maxlon
  region.filter <- oco2.lat >= lon.lat[3] & oco2.lat < lon.lat[4] &
                   oco2.lon >= lon.lat[1] & oco2.lon < lon.lat[2]

  # filter by quality flag too, more receptors when XCO2 is high
  if (oco2.ver == 'b7rb') data.filter <- QF == 0
  if (oco2.ver == 'b8r')  data.filter <- WL <= 1
  sel.lat  <- oco2.lat[region.filter & data.filter]
  sel.lon  <- oco2.lon[region.filter & data.filter]
  sel.xco2 <- xco2[region.filter & data.filter]
  sel.qf   <- QF[region.filter & data.filter]
  sel.wl   <- WL[region.filter & data.filter]
  sel.id   <- id[region.filter & data.filter]

  sel.foot <- oco2.foot[region.filter & data.filter]  # OCO-2 footprint
  sel.time <- substr(id, 1, 10)[region.filter & data.filter]
  sel.YYYY <- YYYY[region.filter & data.filter]
  sel.MM   <- MM[region.filter & data.filter]
  sel.DD   <- DD[region.filter & data.filter]
  sel.HH   <- HH[region.filter & data.filter]
  sel.mm   <- mm[region.filter & data.filter]  # min
  sel.ss   <- ss[region.filter & data.filter]  # seconds

  sel.mode <- mode[region.filter & data.filter]
  cat('Operational Modes:', unique(sel.mode), '\n\n')

  # whether plotting XCO2 from OCO-2
  if (plotTF) {
    xco2.obs <- data.frame(lat = as.numeric(sel.lat), lon = as.numeric(sel.lon),
      xco2 = as.numeric(sel.xco2))

    zoom <- 8; font.size <- rel(1.0)
    col <- c('black', '#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8',
      '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')
    col.range <- seq(380, 420, 2)
    sitemap <- get_map(location = c(lon = lon.lat[5], lat = lon.lat[6]),
      zoom = zoom, maptype = 'roadmap')  # plot google map
    m1 <- ggmap(sitemap) + theme_bw()

    # add observed XCO2
    c1 <- m1 + geom_point(data = xco2.obs, aes(x = lon, y = lat, colour = xco2))
    c1 <- c1 + labs(x = 'LONGITUDE [degE]', y = 'LATITUDE [degN]')
    c1 <- c1 + labs(title=paste('OCO-2 XCO2 [ppm] for', site, 'on', YYYYMMDD))
    c1 <- c1 + scale_colour_gradientn(name = 'OCO-2 XCO2 [ppm]', colours = col,
      breaks = col.range, labels = col.range)

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
    picname <- paste0('ggmap_xco2_', site,'_', YYYYMMDD, '.png')
    ggsave(c2, filename = picname, width = 6, height = 6)
  }  # end if plotTF

  # ------------------- Step 2. SET UP the STILT receptors ------------------- #
  # e.g., Time, Lat, Lon, and Height based on observation
  yr  <- as.numeric(sel.YYYY) - 2000
  mon <- as.numeric(sel.MM)
  day <- as.numeric(sel.DD)
  hr  <- as.numeric(sel.HH)
  min <- as.numeric(sel.mm)
  sec <- as.numeric(sel.ss)

  # select lat, lon based on OCO-2 soundings and data quality
  if(selTF){

    ### 'lat','lon',&'agl' can be a VECTOR of the same length
    sort.lat <- sort(sel.lat)
    recp.lat <- sort.lat[findInterval(recp.indx, sort.lat)]
    match.index <- unique(match(recp.lat, sel.lat))

    # change order for longitude, and time
    recp.lat <- signif(sel.lat[match.index], 6)
    recp.lon <- signif(sel.lon[match.index], 7)
    recp.yr  <- yr [match.index]; recp.mon <- mon[match.index]
    recp.day <- day[match.index]; recp.hr  <- hr [match.index]
    recp.min <- min[match.index]; recp.sec <- sec[match.index]

  }else{   # if no requirement, simulate all soundings, no selection
    recp.lat <- sel.lat; recp.lon <- sel.lon
    recp.yr  <- yr;  recp.mon <- mon
    recp.day <- day; recp.hr  <- hr
    recp.min <- min; recp.sec <- sec
  }  # end if selTF

  # grab simulation timing, yyyy-mm-dd HH:MM:SS (UTC), aka 'receptors' info
  run_times_utc <- paste0(recp.yr + 2000, '-',
    formatC(recp.mon, width = 2, flag = 0), '-',
    formatC(recp.day, width = 2, flag = 0), ' ',
    formatC(recp.hr, width = 2, flag = 0), ':',
    formatC(recp.min, width = 2, flag = 0), ':',
    formatC(recp.sec, width = 2, flag = 0))

  # Expand the run times, latitudes, and longitudes to form the unique receptors
  # that are used for each simulation and match Ben's code
  recp.info <- data.frame(run_time = run_times_utc, lati = recp.lat,
    long = recp.lon, stringsAsFactors = F)
  recp.info <- recp.info[order(recp.info$lati), ]

  # subset receptor data frame, if needed
  if (!is.null(recp.num))
    recp.info <- recp.info[min(recp.num):max(recp.num), ]

  if (!is.null(find.lat))
    recp.info <- recp.info[findInterval(find.lat, recp.info$lati), ]

  recp.info$zagl <- agl
  return(recp.info)    # return receptor info
} # end of subroutine
