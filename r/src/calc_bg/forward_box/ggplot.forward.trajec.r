### subroutine to draw STILT particle locations on 2D or 3D plot
# need joinPolys() and point.in.polygon() functions from PBSmapping
# td for the threshold for cutting polygons
# by Dien Wu, 06/21/2017

# convert kernel density to percentile, DW, 11/15/2017
# return latitude range, DW, 11/16/2017
# use readRDS instead of getr (have changed in Trajecmulti), DW, 07/29/2018
# add customized data filtering, DW, 08/20/2018

# add perc for adjustable background extension length, DW, 08/21/2018
# e.g., perc = 0.2 -> extend derived enhanced lat range by 20% on both side

# add dlat range for background range, DW, 08/24/2018
# add plotting for latitude series, DW, 09/20/2018

# add dlat range for background range, DW, 08/24.2018
# add numbers of soundings used for background, DW, 09/05/2018
# add background uncertainty (including spread sd + retrieval err), DW, 09/07/2018
# bug fixed if there's no intersection between overpass and forward plume, DW, 10/03/2018
# updates for b9r, ALWAYS use obs with QF = 0 for plotting and bg estimates, 
#   get rid of data.filter, DW, 10/30/2018
# merge two plots together if there's an intersection between overpass and plume,
#   DW, 10/30/2018 

ggplot.forward.trajec <- function(ident, trajpath, site, timestr, oco2.path, 
                                  oco2.ver, met, zoom = 7, lon.lat, 
                                  font.size = rel(1.2), td = 0.05, 
                                  bg.dlat = 0.5, perc = 0.2, 
                                  clean.side = c('north', 'south', 'both')[3]){

  # call grab.oco2() to read in observations and compute overpass time
  # lon.lat used for grabbing OCO2 should be wider, e.g., by 2 deg +
  mod.lon.lat <- lon.lat
  mod.lon.lat[, c('minlon', 'minlat')] <- mod.lon.lat[, c('minlon', 'minlat')] - 2
  mod.lon.lat[, c('maxlon', 'maxlat')] <- mod.lon.lat[, c('maxlon', 'maxlat')] + 2

  obs <- grab.oco2(oco2.path, timestr, mod.lon.lat, oco2.ver) %>% filter(qf == 0)
  if (nrow(obs) == 0) {
    cat('*** NO observed data with QF = 0 for this event...returning NA ***\n')
    return()
  }
  
  obs.datestr <- as.POSIXlt(as.character(obs$time),
                            format = '%Y-%m-%d %H:%M:%S', tz = 'UTC')

  # get satellite crossing durations, and further allow for fwe more mins
  # allow for 2 more minutes, convert to sec
  min.xtime <- min(obs.datestr) - 2 * 60
  max.xtime <- max(obs.datestr) + 2 * 60

  # read in forward trajec
  cat('\n\nggplot.forward.trajec(): reading forward trajec, it takes time...\n')
  trajdat <- NULL
  for (f in 1 : length(ident)) {

    if (!file.exists(paste0(trajpath, ident[f]))) {
      cat('***** NO forward trajec found...returning NA *****\n')
      return()
    }

    tmp.trajdat <- readRDS(paste0(trajpath, ident[f]))
    trajdat <- rbind(trajdat, tmp.trajdat)
  } # end for f

  # before subsetting trajec, grab box receptor
  recp.trajdat <- trajdat %>% filter(time == min(time))
  llon <- c(min(recp.trajdat$lon), max(recp.trajdat$lon))
  llat <- c(min(recp.trajdat$lat), max(recp.trajdat$lat))
  box.recp <- data.frame(lon = c(llon[1], llon[2], llon[2], llon[1]),
                         lat = c(llat[1], llat[1], llat[2], llat[2]))

  # only allow for an hour duration, and compute release times in POSIXct format
  info <- ident.to.info(ident = ident, stilt.ver = 1, aglTF = F)[[1]]
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
                      prob = as.vector(dens$z)) %>% 
           mutate(norm.prob = prob / max(prob))

  # plot kernel density on map
  lab.norm <- c(td, seq(0, 1, 0.1)[-1])
  p1 <- m1 + geom_contour(data = densf, aes(x = lon, y = lat, z = norm.prob,
                          colour = ..level..), breaks = lab.norm, size = 1.3) +
             scale_colour_gradient(name = 'Normalized\nKernel\nDensity',
                                   low = 'lightblue', high = 'purple', 
                                   breaks = lab.norm, labels = lab.norm,
                                   limits = c(0, max(lab.norm)))

  ### get plotting info, i.e., different kernel density levels
  kd.info <- ggplot_build(p1)$data[[5]]
  uni.level <- unique(kd.info$level)

  # assign density levels to each selected particle
  sel.trajdat <- sel.trajdat %>% mutate(dens.level = 0)
  for (u in 1:length(uni.level)) {
    tmp.kd <- kd.info %>% filter(level == uni.level[u])
    tmp.index <- point.in.polygon(sel.trajdat$lon, sel.trajdat$lat,
                                  tmp.kd$x, tmp.kd$y)
    sel.trajdat[tmp.index > 0, 'dens.level'] <- uni.level[u]
  } # end for u

  p1 <- p1 + geom_point(data = sel.trajdat, aes(lon, lat, colour = dens.level),
                        size = 0.2, alpha = 0.5)

  # add observed soundings
  max.y <- ceiling(max(obs$xco2)); min.y <- floor(min(obs$xco2))

  # only plot SCREENED data, QF == 0
  p2 <- p1 + geom_point(data = obs, aes(lon, lat, fill = xco2), shape = 21, 
                        colour = 'gray90') +
             scale_fill_gradientn(colours = def.col(), name = 'XCO2',
                                  limits = c(min.y, max.y), 
                                  breaks = seq(min.y, max.y, 2),
                                  labels = seq(min.y, max.y, 2))

  # compute outmost boundary
  wide.bound <- kd.info %>% filter(level == td)
  bound.traj <- wide.bound %>% dplyr::select(X = x, Y = y) %>%
                               mutate(POS = 1:nrow(wide.bound), PID = 1)

  # if diff is too large, meaning there are more than one polygon
  abs.dy <- abs(diff(bound.traj$Y))
  abs.dx <- abs(diff(bound.traj$X))

  # select the polygon we need, bug fixed, DW, 08/21/2018
  if (length(which(abs.dy > td)) > 0 | length(which(abs.dx > td)) > 0){
    cutoff.yindx <- which(abs.dy > td)
    cutoff.xindx <- which(abs.dx > td)

    y.vec <- c(1, cutoff.yindx, nrow(bound.traj)); dy.indx <- diff(y.vec)
    x.vec <- c(1, cutoff.xindx, nrow(bound.traj)); dx.indx <- diff(x.vec)

    # now select polygons
    plg.yindx <- which(dy.indx == max(dy.indx))  # polygon index
    plg.xindx <- which(dx.indx == max(dx.indx))  # polygon index
    plg.xy <- intersect(
      seq(y.vec[plg.yindx], y.vec[plg.yindx + 1], 1),
      seq(x.vec[plg.xindx], x.vec[plg.xindx + 1], 1)
    )
    bound.traj <- bound.traj[plg.xy, ]
  }  # end if

  # add bounrdary of city plume
  p3 <- p2 + geom_polygon(data = bound.traj, aes(X, Y), colour = 'gray10',
                          linetype = 1, fill = NA, size = 0.9, alpha = 0.5)

  title.p <- paste('Forward-time urban plume for overpass on', timestr)
  p3 <- p3 + labs(x = 'LONGITUDE', y = 'LATITUDE', title = title.p)

  # further select OCO2 soundings over OCO-2 track
  sel.obs <- obs %>% filter(lat <= max(bound.traj$Y) & lat >= min(bound.traj$Y) &
                            lon <= max(bound.traj$X) & lon >= min(bound.traj$X))

  if (nrow(sel.obs) == 0) {

    cat('ggplot.forward.trajec(): no intersection with screened obs..return NA\n')
    
    # store plot even if there is no intersection
    picname <- file.path(trajpath, paste0('forward_plume_', site, '_', timestr, 
                                          '_', met, '_', oco2.ver, '.png'))
    ggsave(p3, filename = picname, width = 10, height = 10)
    return()

  } else {

    # find qualified obs
    hpts <- chull(x = sel.obs$lon, y = sel.obs$lat); hpts <- c(hpts, hpts[1])
    narrow.bound <- sel.obs[hpts,]
    bound.obs <- data.frame(X = narrow.bound$lon, Y = narrow.bound$lat,
                            PID = rep(2, nrow(narrow.bound)), 
                            POS = 1:nrow(narrow.bound))

    # find overlapping region and find all points in overlapping polygon
    joint.bound <- joinPolys(bound.traj, bound.obs) # in PBS Mapping package
    pol.bound   <- joint.bound

    # get obs that fall into the polluted lat range
    pol.index <- point.in.polygon(obs$lon, obs$lat, pol.bound$X, pol.bound$Y)
    pol.obs   <- obs[pol.index > 0,]

    # compute mean background for both northern and southern
    if (nrow(pol.obs) > 0) {

      # plot overlap polygon and polluted obs (only screened data), DW, 08/20/2018
      p3 <- p3 + geom_polygon(data = pol.bound, aes(X, Y), colour = 'gray50',
                              fill = 'gray40', alpha = 0.5) +
                 geom_point(data = pol.obs, aes(lon, lat, fill = xco2), 
                           shape = 21, colour = 'gray50') + 
                 annotate('text', x = unique(info$recp.lon) + mm[[3]],
                                  y = unique(info$recp.lat) + mm[[2]] - 0.05, 
                                  label = site, size = 6) +
                 annotate('point', x = info$recp.lon + mm[[3]],
                                   y = info$recp.lat + mm[[2]], 
                                   size = 2, shape = 17)

      p4 <- p3 + theme(legend.position = 'right',
                       legend.key.width = unit(0.5, 'cm'),
                       legend.key.height = unit(1.5, 'cm'),
                       legend.text = element_text(size = font.size),
                       legend.key = element_blank(),
                       axis.title.y = element_text(size = font.size, angle = 90),
                       axis.title.x = element_text(size = font.size, angle = 0),
                       axis.text = element_text(size = font.size),
                       axis.ticks = element_line(size = font.size),
                       title = element_text(size = font.size))

      # get Enhanced latitude range
      pol.min.lat <- min(pol.obs$lat)
      pol.max.lat <- max(pol.obs$lat)

      # allow for some uncertainty, to avoid including high XCO2
      dlat <- pol.max.lat - pol.min.lat
      pol.min.lat <- pol.min.lat - dlat * perc
      pol.max.lat <- pol.max.lat + dlat * perc
      cat(paste('Enhanced lat range:', signif(pol.min.lat, 4), '-',
                signif(pol.max.lat, 4),'N\n'))

      # north part
      north.min.lat <- pol.max.lat
      north.max.lat <- min(lon.lat$maxlat, pol.max.lat + bg.dlat)
      north.obs <- obs %>% filter(lat > north.min.lat & lat < north.max.lat)
      north.bg <- north.obs %>% dplyr::summarize(mean = mean(xco2)) %>% as.numeric()

      # south part
      south.min.lat <- max(lon.lat$minlat, pol.min.lat - bg.dlat)
      south.max.lat <- pol.min.lat
      south.obs <- obs %>% filter(lat > south.min.lat & lat < south.max.lat)
      south.bg <- south.obs %>% dplyr::summarize(mean = mean(xco2)) %>% as.numeric()

      # calculate background:
      if (clean.side == 'both') {

        clean.obs <- rbind(north.obs, south.obs)
        clean.min.lat <- south.min.lat 
        clean.max.lat <- north.max.lat
        cat(paste('Background lat range:',  signif(north.min.lat, 4), '-',
                                            signif(north.max.lat, 4), 'N + ', 
                                            signif(south.min.lat, 4), '-',
                                            signif(south.max.lat, 4), 'N\n'))

      } else if (clean.side == 'north') {

        clean.obs <- north.obs 
        clean.min.lat <- north.min.lat 
        clean.max.lat <- north.max.lat 
        cat(paste('Background lat range:', signif(clean.min.lat, 4), '-',
                                           signif(clean.max.lat, 4), '\n'))

      } else if (clean.side == 'south') {
        
        clean.obs <- south.obs 
        clean.min.lat <- south.min.lat 
        clean.max.lat <- south.max.lat 
        cat(paste('Background lat range:', signif(clean.min.lat, 4), '-',
                                           signif(clean.max.lat, 4), '\n'))

      } else {
        cat('Incorrect input of clean.side\n');return()
      } # end if clean.side

      mean.bg   <- mean(clean.obs$xco2, na.rm = T)
      sd.spread <- sd(clean.obs$xco2, na.rm = T)

      # include retrieval error in background uncert, DW, 09/06/2018
      sd.retriv <- sqrt(mean(clean.obs$xco2.uncert ^ 2))
      bg.sd  <- sqrt(sd.spread ^ 2 + sd.retriv ^ 2)
      num.bg <- nrow(clean.obs)

      cat(paste('North:', signif(north.bg, 5),'ppm; South:',
                          signif(south.bg, 5),'ppm; Final:', 
                          signif(mean.bg, 5),'ppm..\n\n'))

      # ---------------------------------------------------------------------- #
      # plot latitude series as well 
      bg.dat <- data.frame(x = seq(clean.min.lat, clean.max.lat, 0.1), 
                           ymin = mean.bg - bg.sd, 
                           ymax = mean.bg + bg.sd, 
                           y = mean.bg)

      l1 <- ggplot() + theme_bw() +
        geom_ribbon(data = bg.dat, aes(x, ymin = ymin, ymax = ymax),
                    colour = 'limegreen', fill = 'limegreen', alpha = 0.3) + 
        geom_point(data = obs, aes(lat, xco2, fill = as.factor(2)),
                   colour = 'white', shape = 24, size = 3) + 
        geom_point(data = pol.obs, aes(lat, xco2, fill = as.factor(4)),
                   colour = 'white', shape = 24, size = 3) + 
        geom_line(data = bg.dat, 
                  aes(x, y, colour = as.factor(6), linetype = as.factor(6)), 
                  size = 1.1)

      # add smooth spline, make sure clean.obs has data, DW, 10/30/2018
      # if clean.obs or north.obs/south.obs is missing, do not plot smooth spline
      if (clean.side != 'both' & nrow(clean.obs) > 0) 
      l1 <- l1 + geom_smooth(data = clean.obs, 
                      aes(lat, xco2, colour = as.factor(5), linetype = as.factor(5)), 
                      size = 1, se = F)

      if (clean.side == 'both' & nrow(north.obs) * nrow(south.obs) > 0) 
      l1 <- l1 + geom_smooth(data = north.obs, 
                      aes(lat, xco2, colour = as.factor(5), linetype = as.factor(5)), 
                      size = 1.3, se = F) + 
                 geom_smooth(data = south.obs, 
                      aes(lat, xco2, colour = as.factor(5), linetype = as.factor(5)), 
                      size = 1.3, se = F)

      lab <- c('1' = 'All OBS', '2' = 'Screened OBS (QF=0)', 
               '3' = 'All enhanced OBS', 
               '4' = 'Screened enhanced OBS (QF=0)', 
               '5' = 'Smooth splines of screened OBS\nover backgorund range',
               '6' = 'Final overpass-specific background')
      
      fill.val <- c('1' = 'gray80', '2' = 'black', '3' = 'pink', '4' = 'brown', 
                    '5' = 'deepskyblue', '6' = 'darkgreen')
      col.val <- c('1' = 'gray80', '2' = 'black', '3' = 'pink', 
                   '4' = 'brown', '5' = 'deepskyblue', '6' = 'darkgreen')
      lt.val <- c('6' = 4, '5' = 1)
      title.l <- paste('Demostration of overpass-specific background [ppm] for',
                       site, 'on', timestr)

      l2 <- l1 + scale_fill_manual(name = NULL, values = fill.val, labels = lab) + 
                 scale_colour_manual(name = NULL, values = col.val, labels = lab) + 
                 scale_linetype_manual(name = NULL, values = lt.val, labels = lab) +
                 labs(x = 'LATITUDE [deg]', y = 'OBS with QF = 0 [ppm]', title = title.l) + 
                 scale_x_continuous(breaks = seq(-90, 90, 0.4), 
                                    labels = seq(-90, 90, 0.4), 
                                    limits = c(lon.lat$minlat, lon.lat$maxlat)) + 
                 scale_y_continuous(breaks = seq(min.y, max.y, 2), 
                                    labels = seq(min.y, max.y, 2), 
                                    limits = c(min.y, max.y))

      l3 <- l2 + theme(legend.position = 'bottom', 
                       legend.key.width = unit(2, 'cm'),
                       legend.key.height = unit(0.5, 'cm'), 
                       legend.text = element_text(size = font.size),
                       legend.key = element_blank(), 
                       panel.grid.minor=element_blank(),
                       axis.title.y = element_text(size = font.size, angle = 90), 
                       axis.title.x = element_text(size = font.size, angle = 0), 
                       axis.text = element_text(size = font.size), 
                       axis.ticks = element_line(size = font.size),
                       title = element_text(size = font.size)) + 
                guides(fill = guide_legend(nrow = 2, byrow = F),
                       colour = guide_legend(nrow = 2, byrow = TRUE))

      # merge map of forward plume and latitude series, DW, 10/30/2018
      p5 <- ggarrange(plotlist = list(p4, l3), heights = c(3, 2), nrow = 2, 
                      labels = c('a)', 'b)'))
      width <- 10; height <- 13
    } else {  # if no intersection

      cat('ggplot.forward.trajec(): No intersection with OCO-2 track...return NA\n')
      north.bg <- NA; south.bg <- NA; mean.bg  <- NA; bg.sd <- NA
      sd.spread <- NA; sd.retriv <- NA; num.bg <- NA
      pol.min.lat <- NA; pol.max.lat <- NA; north.min.lat <- NA 
      north.max.lat <- NA; south.min.lat <- NA; south.max.lat <- NA

      p5 <- p3  # pass the last fig to p5 and plot it later
      width <- 10; height <- 10
    }  # end if nrow(pol.obs)

    picname <- file.path(trajpath, paste0('forward_plume_', site, '_', timestr, 
                                          '_', met, '_', oco2.ver, '.png'))
    ggsave(p5, filename = picname, width = width, height = height)

    # also, return the max min latitude ranges for polluted range
    bg.info <- data.frame(timestr, north.bg, south.bg, final.bg = mean.bg,
                          final.bg.sd = bg.sd, sd.spread, sd.retriv, num.bg,
                          pol.min.lat, pol.max.lat, clean.side, north.min.lat, 
                          north.max.lat, south.min.lat, south.max.lat)
    return(bg.info)
  } # end if

}  # end of subroutine
