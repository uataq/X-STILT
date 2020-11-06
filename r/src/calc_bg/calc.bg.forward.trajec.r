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
# minor update to the script in using OCO-3 data, DW, 06/30/2020

if (F) {  # for debugging

  timestr = all.timestr[tt]
  sensor = 'TROPOMI'
  oco.ver = 'VEarlyR'
  oco.path = oco.path[tt]
  qfTF = T
  zoom = 7
  td = 0.05
  bg.dlat = 0.5
  perc = 0.2
  clean.side = 'both'

}


calc.bg.forward.trajec <- function(traj.path, site, timestr, sensor, oco.ver, 
                                   oco.path, qfTF = T, met, zoom = 7, td = 0.05,
                                   bg.dlat = 0.5, perc = 0.2, 
                                   clean.side = c('north', 'south', 'both')[3], 
                                   api.key, font.size = rel(1.0)){
  
  library(ggpubr)
  register_google(key = api.key)
  
  # figure name --
  picname <- file.path(traj.path, paste0('forward_plume_', site, '_', timestr, '_', 
                                         met, '_', sensor, '.png'))
  
  cat(paste('\n\n# ---------------------------------#
  calc.bg.forward.trajec(): loading XCO2 obs and forward trajec for', timestr, 
    ', it takes time...\n'))

  # call grab.oco() to read in observations and compute overpass time
  # lon.lat used for grabbing OCO2 should be wider, e.g., by 2 deg +
  lon.lat <- get.lon.lat(site, dlat = 2, dlon = 2)
  obs <- grab.oco(oco.path, timestr, lon.lat, oco.ver) %>% 
         dplyr::select(-c(lons, lats, vertices)) %>% unique() 
  rownames(obs) <- seq(1, nrow(obs))

  # if this overpass is a SAM, DW, 08/04/2020
  if ('Land_SAM' %in% obs$mode) { 
    
    # figure out the different tracks using time and latitude difference
    diff.lat <- abs(diff((obs %>% arrange(time))$lat))
    track.indx <- c(0, which(diff.lat > 0.5)) + 1
    track.indx <- track.indx[track.indx != nrow(obs)]
    obs <- obs %>% mutate(group = findInterval(rownames(obs), track.indx))
    width = 10
  } else {
    obs <- obs %>% mutate(group = 1)   # only one group for normal non-SAM
    width = 8
  } # end if

  if (qfTF) obs <- obs %>% filter(qf == 0)

  # sanity check 
  #o1 <- ggplot(data = obs) + geom_point(aes(lon, lat, color = as.factor(group)))


  if (nrow(obs) == 0) {
    cat('*** NO observed data with QF = 0 for this event...returning NA ***\n')
    return() 
  }
  obs.datestr <- as.POSIXlt(as.character(obs$time), 'UTC', format = '%Y-%m-%d %H:%M:%S')

  # get satellite overpassing durations for selecting forward-time trajec
  min.xtime <- min(obs.datestr) - 1 * 60  # allow for 1 more minutes, convert to sec
  max.xtime <- max(obs.datestr) + 1 * 60

  # ---------------------------------------------------------------------------- 
  # read in forward trajec
  traj.files <- list.files(file.path(traj.path, paste0('out_forward_', timestr, '_', site)), 
                          '.rds', full.names = T, recursive = T)

  if (length(traj.files) == 0) {
    cat('*** NO forward trajec found...returning NA ***\n'); return() }

  trajdat <- NULL
  for (f in traj.files) {
    if (!file.exists(f)) { 
      cat('*** NO forward trajec found...returning NA ***\n'); return() }
    trajdat <- rbind(trajdat, readRDS(f))
  } # end for f

  # before subsetting trajec, grab box receptor
  recp.trajdat <- trajdat %>% filter(time == min(time))
  llon <- c(min(recp.trajdat$lon), max(recp.trajdat$lon))
  llat <- c(min(recp.trajdat$lat), max(recp.trajdat$lat))
  box.recp <- data.frame(lon = c(llon[1], llon[2], llon[2], llon[1]),
                         lat = c(llat[1], llat[1], llat[2], llat[2]))


  # only allow for an hour duration, and compute release times in POSIXct format
  info <- ident.to.info(ident = basename(traj.files), stilt.ver = 1, aglTF = F)[[1]]
  info$rel.date <- as.POSIXlt(info$timestr, format = '%Y%m%d%H%M', tz = 'UTC')

  # get first release time and select trajec based on overpass duration
  min.datestr <- min(info$rel.date)
  sel.trajdat <- trajdat %>% mutate(datestr = min.datestr + time * 60) %>%
                             filter(datestr >= min.xtime & datestr <= max.xtime)

  if (nrow(sel.trajdat) == 0) 
    stop('*** NO trajectories exist during the overpass time\nlikely caused by incomplete traj files, please check...\n')


  # ---------------------------------------------------------------------------- 
  # load map
  #clon = unique(info$recp.lon); clat = unique(info$recp.lat)
  clon = mean(sel.trajdat$lon); clat = mean(sel.trajdat$lat)
  mm <- ggplot.map(map = 'ggmap', center.lon = clon, center.lat = clat, zoom = zoom)

  ### calculate 2D kernel density and Normalized by the max density
  cat('calc.bg.forward.trajec(): calculating kernel density...\n')
  dens <- kde2d(sel.trajdat$lon, sel.trajdat$lat, h = c(0.1, 0.1), n = 100)
  densf <- data.frame(expand.grid(lon = dens$x, lat = dens$y),
                      prob = as.vector(dens$z)) %>% 
           mutate(norm.prob = prob / max(prob))

  # plot kernel density on map and draw receptor box
  lab.norm <- c(td, seq(0, 1, 0.1)[-1])
  p1 <- mm[[1]] + 
        geom_contour(data = densf, aes(lon, lat, z = norm.prob, colour = ..level..), 
                     breaks = lab.norm, size = 1.3) +
        scale_colour_gradient(name = 'Normalized\nKernel\nDensity',
                              low = 'lightblue', high = 'purple', 
                              breaks = lab.norm, labels = lab.norm,
                              limits = c(0, max(lab.norm)))

  ### get plotting info, i.e., different kernel density levels
  kd.info <- ggplot_build(p1)$data[[4]]
  uni.level <- unique(kd.info$level)

  # assign density levels to each selected particle
  sel.trajdat <- sel.trajdat %>% mutate(dens.level = 0)
  for (u in 1:length(uni.level)) {
      tmp.kd <- kd.info %>% filter(level == uni.level[u])
      tmp.index <- sp::point.in.polygon(sel.trajdat$lon, sel.trajdat$lat,
                                        tmp.kd$x, tmp.kd$y)
      sel.trajdat[tmp.index > 0, 'dens.level'] <- uni.level[u]
  } # end for u


  # add air parcels and observations 
  max.y <- ceiling(max(obs$xco2)); min.y <- floor(min(obs$xco2))
  uni.modes <- unique(obs$mode)
  if (length(uni.modes) > 1) uni.modes <- toString(uni.modes)

  title.p <- paste('Forward-time urban plume for', uni.modes, 'on', timestr)
  p2 <- p1 + geom_point(data = sel.trajdat, aes(lon, lat, colour = dens.level),
                        size = 0.2, alpha = 0.4) + 
        geom_polygon(data = box.recp, aes(lon, lat), color = 'darkgreen', 
                     size = 1.2, linetype = 3, fill = NA) +
        geom_point(data = obs, aes(lon, lat, fill = xco2), shape = 21, colour = 'gray90') +
        scale_fill_gradientn(colours = def.col(), name = 'XCO2',
                             limits = c(min.y, max.y), 
                             breaks = seq(min.y, max.y, 2),
                             labels = seq(min.y, max.y, 2)) + 
        labs(x = 'LONGITUDE', y = 'LATITUDE', title = title.p)
  

  # ---------------------- compute outmost boundary -------------------------
  wide.bound <- kd.info %>% filter(level == td)
  bound.traj <- wide.bound %>% dplyr::select(X = x, Y = y) %>%
                               mutate(POS = 1:nrow(wide.bound), PID = 1)

  # if diff is too large, meaning there are more than one polygon
  abs.dy <- abs(diff(bound.traj$Y))
  abs.dx <- abs(diff(bound.traj$X))

  # there could be multiple polygons, so select the one we need, DW, 08/21/2018
  if (length(which(abs.dy > td)) > 0 | length(which(abs.dx > td)) > 0){
      
      cutoff.yindx <- which(abs.dy > td)
      cutoff.xindx <- which(abs.dx > td)

      y.vec <- c(1, cutoff.yindx, nrow(bound.traj)); dy.indx <- diff(y.vec)
      x.vec <- c(1, cutoff.xindx, nrow(bound.traj)); dx.indx <- diff(x.vec)

      # now select polygons
      plg.yindx <- which(dy.indx == max(dy.indx))  # polygon index
      plg.xindx <- which(dx.indx == max(dx.indx))  # polygon index

      # it is possible that there are more than 1 polygon that matches, 
      # use the first one that found for now, DW, 07/01/2020
      if (length(plg.xindx) * length(plg.yindx) > 1) {
        for (xi in plg.xindx) {
          for (yi in plg.yindx) {
            plg.xy <- intersect(seq(y.vec[yi], y.vec[yi + 1], 1),
                                seq(x.vec[xi], x.vec[xi + 1], 1))
            if (length(plg.xy) > 0) break   # if meet the criteria then break the loop
          }  # end for yi
        }  # end for xi
      } else plg.xy <- intersect(seq(y.vec[plg.yindx], y.vec[plg.yindx + 1], 1),
                                 seq(x.vec[plg.xindx], x.vec[plg.xindx + 1], 1))

      bound.traj <- bound.traj[plg.xy, ]
  }  # end if

  p3 <- p2 + geom_polygon(data = bound.traj, aes(X, Y), colour = 'gray10',
                          linetype = 1, fill = NA, size = 0.9, alpha = 0.5)

  # further select OCO2 soundings over OCO-2 track
  sel.obs <- obs %>% filter(lat <= max(bound.traj$Y) & lat >= min(bound.traj$Y) &
                            lon <= max(bound.traj$X) & lon >= min(bound.traj$X))

  if (nrow(sel.obs) == 0) {  # store plot even if there is no intersection
    cat('calc.bg.forward.trajec(): *** no intersection with screened obs\ntry adjusting @param td...return NA\n')
    ggsave(p3, filename = picname, width = width, height = 10); return()
  } # end if nrow(sel.obs) == 0


  # ---------------------------------------------------------------------------- 
  # calculate background since intersection between swath and trajec exists ----
  # compute convex hull
  hpts <- chull(x = sel.obs$lon, y = sel.obs$lat); hpts <- c(hpts, hpts[1])
  narrow.bound <- sel.obs[hpts,]
  bound.obs <- data.frame(X = narrow.bound$lon, Y = narrow.bound$lat,
                          PID = rep(2, nrow(narrow.bound)), 
                          POS = 1:nrow(narrow.bound))

  # find overlapping region and find all points in overlapping polygon
  joint.bound <- PBSmapping::joinPolys(bound.traj, bound.obs) # in PBS Mapping package
  pol.bound   <- joint.bound

  # compute mean background for both northern and southern
  if (!is.null(pol.bound)) {

    # get obs that fall into the polluted lat range
    pol.index <- point.in.polygon(obs$lon, obs$lat, pol.bound$X, pol.bound$Y)
    pol.obs   <- obs[pol.index > 0,]

    # plot overlap polygon and polluted obs (only screened data), DW, 08/20/2018
    p4 <- p3 + geom_polygon(data = pol.bound, aes(X, Y), colour = 'gray20',
                            fill = 'gray60', alpha = 0.5) +
               geom_point(data = pol.obs, aes(lon, lat, fill = xco2), 
                          shape = 21, colour = 'gray50') + 
               annotate('text', x = unique(info$recp.lon) + mm[[3]],
                                y = unique(info$recp.lat) + mm[[2]] - 0.05, 
                                label = site, size = 5) +
               annotate('point', x = info$recp.lon + mm[[3]], size = 2, 
                                 y = info$recp.lat + mm[[2]], shape = 17) + 
               theme(legend.position = 'right', legend.key = element_blank(),
                      legend.key.width = unit(0.5, 'cm'),
                      legend.key.height = unit(1.2, 'cm'),
                      legend.text = element_text(size = font.size),
                      axis.title.y = element_text(size = font.size, angle = 90),
                      axis.title.x = element_text(size = font.size, angle = 0),
                      axis.text = element_text(size = font.size),
                      axis.ticks = element_line(size = font.size),
                      title = element_text(size = font.size))


    # get urban enhanced lat range for each track group for SAM, DW, 08/04/2020
    pol.loc <- pol.obs %>% group_by(group) %>% 
               dplyr::summarize(pol.min.lat = min(lat), pol.max.lat = max(lat), 
                                dlat = pol.max.lat - pol.min.lat, 

                                # adjust enhanced lat range if buffer exists (i.e., dlat * perc)
                                pol.min.lat = pol.min.lat - dlat * perc,
                                pol.max.lat = pol.max.lat + dlat * perc)

    ### calculate multiple background based on 'group' column in obs
    # multiple groups for SAMs; only one group for non-SAM, DW, 08/04/2020
    pol.info <- obs %>% group_by(group) %>% left_join(pol.loc, by = 'group') %>%
                mutate(# add lat boundaries for northern and southern parts outside urban plume
                       north.min.lat = pol.max.lat, 
                       north.max.lat = min(lon.lat$maxlat, pol.max.lat + bg.dlat), 
                       south.min.lat = max(lon.lat$minlat, pol.min.lat - bg.dlat), 
                       south.max.lat = pol.min.lat, 

                       # before calculating the backgorund, label the type of obs as
                       # e.g., north.clean and south.clean for north/south background obs
                       #       urban for urban enhanced obs 
                       label = ifelse(lat > north.min.lat & lat < north.max.lat, 'north.clean', 
                               ifelse(lat > south.min.lat & lat < south.max.lat, 'south.clean', 'urban'))
                      ) %>% ungroup()
    
    # sanity check 
    #o2 <- ggplot(data = pol.info) + geom_point(aes(lon, lat, color = label))

    # calculate background for either nothern or southern background for each track group
    if (clean.side == 'north') {
      clean.obs <- pol.info %>% group_by(group) %>% filter(label == 'north.clean')
      bg.info <- clean.obs %>% dplyr::summarize(clean.min.lat = min(lat), 
                                                clean.max.lat = max(lat), 
                                                bg.mean = mean(xco2), 
                                                bg.sd.spread = sd(xco2), 
                                                bg.sd.retrv = sqrt(mean(xco2.uncert ^ 2))) %>% 
                               left_join(clean.obs %>% tally(), by = 'group')

    } else if (clean.side == 'south') {
      clean.obs <- pol.info %>% group_by(group) %>% filter(label == 'south.clean') 
      bg.info <- clean.obs %>% dplyr::summarize(clean.min.lat = min(lat), 
                                                clean.max.lat = max(lat), 
                                                bg.mean = mean(xco2), 
                                                bg.sd.spread = sd(xco2), 
                                                bg.sd.retrv = sqrt(mean(xco2.uncert ^ 2))) %>% 
                               left_join(clean.obs %>% tally(), by = 'group')

    } else if (clean.side == 'both') {
      clean.obs <- pol.info %>% group_by(group) %>% filter(label != 'urban') 
      bg.info <- clean.obs %>% dplyr::summarize(clean.min.lat = min(lat), 
                                                clean.max.lat = max(lat), 
                                                bg.mean = mean(xco2), 
                                                bg.sd.spread = sd(xco2), 
                                                bg.sd.retrv = sqrt(mean(xco2.uncert ^ 2))) %>% 
                               left_join(clean.obs %>% tally(), by = 'group')

    } else stop('Wrong input of @param clean.side')

    if (nrow(clean.obs) == 0) { 
      cat('calc.bg.forward.trajec(): NO clean OBS found, try another `clean.side`...\n')
      ggsave(p3, filename = picname, width = width, height = 10)
      return()
    }

    # include retrieval error in background uncert, DW, 09/06/2018
    bg.info <- bg.info %>% mutate(bg.sd = sqrt(bg.sd.spread ^ 2 + bg.sd.retrv ^ 2), 
                                  clean.side = clean.side, timestr = timestr) %>% 
               relocate(timestr, .before = group) %>% left_join(pol.loc, by = 'group')


    # ---------------------------------------------------------------------- #
    # plot latitude series as well 
    lab <- c('1' = 'All OBS', '2' = 'Screened OBS (QF=0)', 
             '3' = 'All enhanced OBS', '4' = 'Screened enhanced OBS (QF=0)', 
             '5' = 'Smooth splines of screened OBS\nover backgorund range',
             '6' = 'Overpass-specific MEAN background +\nuncertainty in lightgreen ribbon')
    
    fill.val <- c('1' = 'gray80', '2' = 'black', '3' = 'pink', '4' = 'brown', 
                  '5' = 'deepskyblue', '6' = 'darkgreen')
    col.val <- c('1' = 'gray80', '2' = 'black', '3' = 'pink', '4' = 'brown', 
                 '5' = 'deepskyblue', '6' = 'darkgreen')
    lt.val <- c('6' = 4, '5' = 1, '7' = 1)

    l1 <- ggplot() + theme_bw() + facet_wrap(~group, nrow = 1) + 
          geom_rect(data = bg.info, aes(xmin = clean.min.lat, xmax = clean.max.lat, 
                                        ymin = bg.mean - bg.sd, ymax = bg.mean + bg.sd), 
                    alpha = 0.3, fill = 'limegreen', color = 'limegreen') + 
          geom_point(data = obs %>% filter(qf == 0), aes(lat, xco2, fill = as.factor(2)),
                     colour = 'white', shape = 24, size = 2) +   # plot screened obs
          geom_point(data = pol.obs %>% filter(qf == 0), aes(lat, xco2, fill = as.factor(4)),
                     colour = 'white', shape = 24, size = 2) + 
          geom_hline(data = bg.info, aes(yintercept = bg.mean, colour = as.factor(6)), 
                     linetype = 4, size = 1.1)

    # remove smooth splines for background obs
    title.l <- paste('Overpass-specific background [ppm] for', site, 'on', timestr)
    l2 <- l1 + scale_fill_manual(name = NULL, values = fill.val, labels = lab) + 
               scale_colour_manual(name = NULL, values = col.val, labels = lab) + 
               #scale_linetype_manual(name = NULL, values = lt.val, labels = lab) +
               labs(x = 'LATITUDE [deg]', y = 'OBS XCO2 [ppm]', title = title.l) + 
               scale_x_continuous(breaks = seq(-90, 90, 0.4), 
                                  labels = seq(-90, 90, 0.4)) +  
                                  #limits = c(lon.lat$minlat, lon.lat$maxlat)) + 
               scale_y_continuous(breaks = seq(min.y, max.y, 2), 
                                  labels = seq(min.y, max.y, 2), 
                                  limits = c(min.y, max.y))

    l3 <- l2 + theme(legend.position = 'bottom', legend.key = element_blank(), 
                     legend.key.width = unit(1.8, 'cm'),
                     legend.key.height = unit(0.5, 'cm'), 
                     legend.text = element_text(size = font.size),
                     panel.grid.minor = element_blank(),
                     axis.title.y = element_text(size = font.size, angle = 90), 
                     axis.title.x = element_text(size = font.size, angle = 0), 
                     axis.text = element_text(size = font.size), 
                     axis.ticks = element_line(size = font.size),
                     title = element_text(size = font.size)) + 
               guides(fill = guide_legend(nrow = 2, byrow = F),
                      colour = guide_legend(nrow = 2, byrow = TRUE), linetype = F)

    # merge map of forward plume and latitude series, DW, 10/30/2018
    pl <- ggarrange(plotlist = list(p4, l3), heights = c(3, 1.5), nrow = 2, 
                    labels = c('a)', 'b)'))
    bg.info <- bg.info %>% dplyr::select(-c(clean.min.lat, clean.max.lat))

    ggsave(pl, filename = picname, width = width, height = 10)
    return(bg.info)

  } else {  # if no enhanced obs

    cat('calc.bg.forward.trajec(): No intersection with OCO-2 track...return NA\n')
    ggsave(p3, filename = picname, width = width, height = 10)
    return()
  }  # end if is.null(pol.bound)

}  # end of subroutine









if (F) {

    # original scripts --------
    pol.min.lat <- min(pol.obs$lat); pol.max.lat <- max(pol.obs$lat)
    dlat <- pol.max.lat - pol.min.lat       # enhanced lat range in deg
    pol.min.lat <- pol.min.lat - dlat * perc
    pol.max.lat <- pol.max.lat + dlat * perc
    cat(paste('Enhanced lat range:', signif(pol.min.lat, 4), '-', signif(pol.max.lat, 4),'N\n'))

    # north part
    north.min.lat <- pol.max.lat
    north.max.lat <- min(lon.lat$maxlat, pol.max.lat + bg.dlat)
    north.obs <- obs %>% filter(lat > north.min.lat & lat < north.max.lat)
    north.bg  <- north.obs %>% dplyr::summarize(mean = mean(xco2)) %>% as.numeric()

    # south part
    south.min.lat <- max(lon.lat$minlat, pol.min.lat - bg.dlat)
    south.max.lat <- pol.min.lat
    south.obs <- obs %>% filter(lat > south.min.lat & lat < south.max.lat)
    south.bg  <- south.obs %>% dplyr::summarize(mean = mean(xco2)) %>% as.numeric()

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
      cat('calc.bg.forward.trajec(): Incorrect input of clean.side\n');return()
    } # end if clean.side


   # arrange results
    mean.bg <- mean(clean.obs$xco2, na.rm = T)
    sd.spread <- sd(clean.obs$xco2, na.rm = T)

    # include retrieval error in background uncert, DW, 09/06/2018
    sd.retrv <- sqrt(mean(clean.obs$xco2.uncert ^ 2))
    bg.sd  <- sqrt(sd.spread ^ 2 + sd.retrv ^ 2)
    num.bg <- nrow(clean.obs)

    # ---------------------------------------------------------------------- #
    # plot latitude series as well 
    bg.dat <- data.frame(x = seq(clean.min.lat, clean.max.lat, 0.1), 
                         ymin = mean.bg - bg.sd, ymax = mean.bg + bg.sd, 
                         y = mean.bg)


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

north.bg <- NA; south.bg <- NA; mean.bg  <- NA; bg.sd <- NA
    sd.spread <- NA; sd.retrv <- NA; num.bg <- NA
    pol.min.lat <- NA; pol.max.lat <- NA; north.min.lat <- NA 
    north.max.lat <- NA; south.min.lat <- NA; south.max.lat <- NA

  # also, return the max min latitude ranges for polluted range
  bg.info <- data.frame(timestr, north.bg, south.bg, final.bg = mean.bg,
                        final.bg.sd = bg.sd, sd.spread, sd.retrv, num.bg,
                        pol.min.lat, pol.max.lat, clean.side, north.min.lat, 
                        north.max.lat, south.min.lat, south.max.lat)

}