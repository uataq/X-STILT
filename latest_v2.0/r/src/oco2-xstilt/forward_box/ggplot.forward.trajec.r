### subroutine to draw STILT particle locations on 2D or 3D plot
# need joinPolys() and point.in.polygon() functions from PBSmapping
# td for the threshold for cutting polygons
# by Dien Wu, 06/21/2017

# convert kernel density to percentile, DW, 11/15/2017
# return latitude range, DW, 11/16/2017
# use readRDS instead of getr (have changed in Trajecmulti), DW, 07/29/2018

ggplot.forward.trajec <- function(ident, trajpath = outpath, site, timestr,
  oco2.path, oco2.ver, zoom = 8, lon.lat, font.size = rel(1.2), td = 0.05,
  clean.side = c('north','south', 'both')[3]){

  # call grab.oco2() to read in observations and compute overpass time
  # lon.lat used for grabbing OCO2 should be wider, e.g., by 2 deg +
  mod.lon.lat <- lon.lat
  mod.lon.lat[c(1, 3)] <- mod.lon.lat[c(1, 3)] - 2
  mod.lon.lat[c(2, 4)] <- mod.lon.lat[c(2, 4)] + 2
  obs <- grab.oco2(oco2.path, timestr, mod.lon.lat)
  obs.datestr <- as.POSIXlt(as.character(obs$time),
    format = '%Y-%m-%d %H:%M:%S', tz = 'UTC')

  # get overpass durations, and further allow for fwe more mins
  # allow for 2 more minutes, convert to sec
  min.xtime <- min(obs.datestr) - 2 * 60
  max.xtime <- max(obs.datestr) + 2 * 60

  # read in forward trajec
  cat('\n\nggplot.forward.trajec(): reading forward trajec, it takes time...\n')
  trajdat <- NULL
  for(f in 1:length(ident)){
    tmp.trajdat <- readRDS(paste0(trajpath, ident[f]))
    trajdat <- rbind(trajdat, tmp.trajdat)
  }

  # before subsetting trajec, grab box receptor
  recp.trajdat <- trajdat %>% filter(time == min(time))
  llon <- c(min(recp.trajdat$lon), max(recp.trajdat$lon))
  llat <- c(min(recp.trajdat$lat), max(recp.trajdat$lat))
  box.recp <- data.frame(lon = c(llon[1], llon[2], llon[2], llon[1]),
                         lat = c(llat[1], llat[1], llat[2], llat[2]))

  # only allow for an hour duration, and compute release times in POSIXct format
  info <- ident.to.info(ident = ident, aglTF = F)[[1]]
  info$rel.date <- as.POSIXlt(info$timestr, format = '%Y%m%d%H%M', tz = 'UTC')

  # get first release time
  min.datestr <- min(info$rel.date)

  # compute POSIXct time each particle and crop traj based on overpass duration
  sel.trajdat <- trajdat %>% mutate(datestr = min.datestr + time * 60) %>%
    filter(datestr >= min.xtime & datestr <= max.xtime)

  # load map
  mm <- ggplot.map(map = 'ggmap', center.lon = unique(info$recp.lon),
    center.lat = unique(info$recp.lat), zoom = zoom)

  # draw receptor box
  m1 <- mm[[1]] + geom_polygon(data = box.recp, aes(lon, lat), linetype = 3,
    fill = 'gray30', alpha = 0.5)

  ### calculate 2D kernel density and Normalized by the max density
  cat('ggplot.forward.trajec(): calculating kernel density...\n')
  dens <- kde2d(sel.trajdat$lon, sel.trajdat$lat, h = c(0.1, 0.1), n = 100)
  densf <- data.frame(expand.grid(lon = dens$x, lat = dens$y),
    prob = as.vector(dens$z)) %>% mutate(norm.prob = prob / max(prob))

  # plot kernel density on map
  lab.norm <- c(td, seq(0, 1, 0.1)[-1])
  p1 <- m1 + geom_contour(data = densf, aes(x = lon, y = lat, z = norm.prob,
      colour = ..level..), breaks = lab.norm, size = 1.3) +
    scale_colour_gradient(name = 'Normalized\nKernel\nDensity',
      low = 'lightblue', high = 'purple', breaks = lab.norm, labels = lab.norm,
      limits = c(0, max(lab.norm)))

  ### get plotting info, i.e., different kernel density levels
  kd.info <- ggplot_build(p1)$data[[5]]
  uni.level <- unique(kd.info$level)

  # assign density levels to each selected particle
  sel.trajdat <- sel.trajdat %>% mutate(dens.level = 0)
  for(u in 1:length(uni.level)){
    tmp.kd <- kd.info %>% filter(level == uni.level[u])
    tmp.index <- point.in.polygon(sel.trajdat$lon, sel.trajdat$lat,
      tmp.kd$x, tmp.kd$y)
    sel.trajdat[tmp.index > 0, 'dens.level'] <- uni.level[u]
  }

  p1 <- p1 + geom_point(data = sel.trajdat, aes(lon, lat, colour = dens.level),
    size = 0.4, alpha = 0.5)

  # add observed soundings
  max.y <- 406; min.y <- 394
  if(timestr >= '20160101'){max.y <- 410; min.y <- 400}

  p2 <- p1 + geom_point(data = obs, aes(lon, lat, fill = xco2), shape = 21,
      colour = 'gray90') +
    scale_fill_gradientn(colours = def.col(), name = 'XCO2',
      limits = c(min.y, max.y), breaks = seq(min.y, max.y, 2),
      labels = seq(min.y, max.y, 2))

  # compute outmost boundary
  wide.bound <- kd.info %>% filter(level == td)
  bound.traj <- wide.bound %>% dplyr::select(X = x, Y = y) %>%
    mutate(POS = 1:nrow(wide.bound), PID = 1)

  # if diff is too large, meaning there are more than one polygon
  if (length(which(abs(diff(bound.traj$Y)) > td)) > 0)
    bound.traj <- bound.traj[1:which(abs(diff(bound.traj$Y)) > td),]

  if (length(which(abs(diff(bound.traj$X)) > td)) > 0)
    bound.traj <- bound.traj[1:which(abs(diff(bound.traj$X)) > td),]

  p3 <- p2 + geom_polygon(data = bound.traj, aes(X, Y), colour = 'gray10',
    linetype = 1, fill = NA, size = 0.9, alpha = 0.5)

  # further select OCO2 soundings over OCO-2 track
  sel.obs <- obs %>%
    filter(lat <= max(bound.traj$Y) & lat >= min(bound.traj$Y) &
           lon <= max(bound.traj$X) & lon >= min(bound.traj$X))

  if (nrow(sel.obs) == 0) {

    cat('ggplot.forward.trajec(): no intersection with obs..return ggplot\n')
    return(p3)

  } else {

    # find qualified obs
    hpts <- chull(x = sel.obs$lon, y = sel.obs$lat); hpts <- c(hpts, hpts[1])
    narrow.bound <- sel.obs[hpts,]
    bound.obs <- data.frame(X = narrow.bound$lon, Y = narrow.bound$lat,
      PID = rep(2, nrow(narrow.bound)), POS = 1:nrow(narrow.bound))

    # find overlapping region and find all points in overlapping polygon
    joint.bound <- joinPolys(bound.traj, bound.obs) # in PBS Mapping package
    pol.bound   <- joint.bound

    # get obs that fall into the polluted lat range
    pol.index <- point.in.polygon(obs$lon, obs$lat, pol.bound$X, pol.bound$Y)
    pol.obs   <- obs[pol.index >0,]

    # plot overlap polygon and polluted obs
    p3 <- p3 + geom_polygon(data = pol.bound, aes(X, Y), colour = 'gray50',
        fill = 'gray40', alpha = 0.5) +
      geom_point(data = pol.obs, aes(lon, lat, fill = xco2),
        shape = 21, colour = 'gray50')

    p4 <- p3 +
      annotate('text', x = unique(info$recp.lon) + mm[[3]],
        y = unique(info$recp.lat) + mm[[2]] - 0.05, label = site, size = 6) +
      annotate('point', x = info$recp.lon + mm[[3]],
        y = info$recp.lat + mm[[2]], size = 2, shape = 17)

    p5 <- p4 + theme(legend.position = 'right',
      legend.key.width = unit(0.5, 'cm'),
      legend.key.height = unit(1.5, 'cm'),
      legend.text = element_text(size = font.size), legend.key = element_blank(),
      axis.title.y = element_text(size = font.size,angle = 90),
      axis.title.x = element_text(size = font.size,angle = 0),
      axis.text = element_text(size = font.size),
      axis.ticks = element_line(size = font.size),
      title = element_text(size = font.size))

    # compute mean background for both northern and southern
    if (nrow(pol.obs) > 0) {

      # screen obs data
      if (oco2.ver == 'b7rb') scr.obs <- obs %>% filter(qf == 0)
      if (oco2.ver == 'b8r') scr.obs <- obs %>% filter(wl <= 1)

      # get Enhanced latitude range
      min.lat <- min(pol.obs$lat)
      max.lat <- max(pol.obs$lat)

      # allow for some uncertainty, to avoid including high XCO2
      dlat <- max.lat - min.lat
      min.lat <- min.lat - dlat * 0.2
      max.lat <- max.lat + dlat * 0.2
      cat(paste('Enhanced lat range:', signif(min.lat, 4), '-',
        signif(max.lat, 4),'N\n'))

      north.obs <- scr.obs %>%
        filter(lat > max.lat & lat < min(lon.lat[4], max.lat + 1))
      south.obs <- scr.obs %>%
        filter(lat < min.lat & lat > max(lon.lat[3], min.lat - 1))

      north.bg <- mean(north.obs$xco2)
      south.bg <- mean(south.obs$xco2)

      # calculate background:
      if (clean.side == 'both') {
        mean.bg <- c(north.bg + south.bg)/2
        sd.bg <- sd(c(north.obs$xco2, south.obs$xco2))

      } else if (clean.side == 'north') {
        mean.bg <- north.bg; sd.bg <- sd(north.obs$xco2)

      } else if(clean.side == 'south') {
        mean.bg <- south.bg; sd.bg <- sd(south.obs$xco2)

      } else {
        cat('Incorrect input of clean.side\n')
      } # end if clean.side

      cat(paste('North:', signif(north.bg, 5),'ppm; South:',
          signif(south.bg, 5),'ppm; Final:', signif(mean.bg, 5),'ppm..\n\n'))

    }else{  # if no intersection
      cat('ggplot.forward.trajec(): No intersection with OCO-2 track...\n')
      north.bg <- NA; south.bg <- NA; mean.bg  <- NA
      min.lat  <- NA; max.lat  <- NA
    }  # end if nrow(pol.obs)

    title <- paste('Forward-time urban plume for overpass on', timestr)
    p5 <- p5 + labs(x = 'LONGITUDE', y = 'LATITUDE', title = title)

    picname <- paste0('urban_plume_forward_', site, '_', timestr, '.png')
    picname <- file.path(trajpath, picname)
    ggsave(p5, filename = picname, width = 12, height = 12)

    # also, return the max min latitude ranges for polluted range
    bg.info <- data.frame(timestr, north.bg, south.bg, final.bg = mean.bg,
      min.lat, max.lat)
    return(bg.info)
  } # end if

}  # end of subroutine
