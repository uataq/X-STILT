# script to plot footprint with observed XCO2 on spatial maps,
# written by Dien Wu, 06/18/2018

# add the sum of foot in plotted region, DW, 07/18/2018
# last update, DW, 08/22/2018
# if no OCO-2 path, do not plot observed XCO2, DW, 10/29/2018

ggmap.xfoot.obs <- function(mm, lon.lat, site, oco2.ver, oco2.path = NULL, 
                            facet.nrow, nhrs, dpar, foot.sf, zisf, met, 
                            stilt.ver, timestr, font.size = rel(0.9), recp.lon, 
                            recp.lat, foot, min.foot.sig = 1E-6, 
                            max.foot.sig = 1E-2, titleTF = T, sumTF = T, 
                            qfTF = T, picname, storeTF = T, width = 12, 
                            height = 8, anthromesTF = F, anthro.path = NULL, 
                            leg.pos = c('bottom', 'right')[1]){

  col <- def.col()
  m1 <- mm[[1]] + theme_bw() + coord_equal(1.1)

  # grab observations using map lat/lon
  map.ext <- data.frame(minlon = min(mm[[1]]$data$lon),
                        maxlon = max(mm[[1]]$data$lon),
                        minlat = min(mm[[1]]$data$lat),
                        maxlat = max(mm[[1]]$data$lat))

  if (!is.null(oco2.path)) {
    cat('Reading OCO-2 data according to the spatial domain of ggmap...\n')
    obs <- grab.oco2(ocopath = oco2.path, timestr, lon.lat = map.ext, oco2.ver)
    qf.obs <- obs %>% filter(qf == 0)
  }

  # select footprints using map.ext
  sel.foot <- foot %>% filter(lon >= map.ext$minlon & lon <= map.ext$maxlon &
                              lat >= map.ext$minlat & lat <= map.ext$maxlat &
                              foot >= min.foot.sig)

  title <- paste0('Spatial time-integrated weighted column footprint (', nhrs, 
                  ' hours; ', dpar, ' dpar; ', foot.sf, '; ziscale = ', zisf,
                  '; met = ', met, ')\nusing STILT version', stilt.ver, 
                  ' for overpass on ', timestr, ' over ', site,
                  '\nOnly large footprints > ', min.foot.sig, ' are displayed')
  if (titleTF == F) title <- NULL

  p1 <- m1 + labs(title = title, x = 'LONGITUDE [E]', y = 'LATITUDE [N]')

  # if there is 1+ receptors
  if ('fac' %in% colnames(sel.foot)){
    if (length(unique(sel.foot$fac)) > 1){

      # summing the values within map domain
      sum.foot <- sel.foot %>% group_by(fac) %>% dplyr::summarize(sum = sum(foot))

      # receptor locations and add receptors on map
      sel.recp <- data.frame(lon = recp.lon, lat = recp.lat,
                             fac = unique(sel.foot$fac)) %>% 
                  full_join(sum.foot, by = 'fac') %>%
                  mutate(x = map.ext$maxlon - 0.5, y = map.ext$maxlat - 0.5)
      print(sel.recp)

      p1 <- p1 +
        geom_point(data = sel.recp, aes(lon, lat), colour = 'purple', size = 3,
                   shape = 1) +
        facet_wrap(~ fac, nrow = facet.nrow) +
        geom_text(data = sel.recp, aes(lon + 0.5, lat), colour = 'purple', 
                  size = 4, label = 'receptor', fontface = 1) 

      if (sumTF) p1 <- p1 + geom_text(data = sel.recp,
                                      aes(x, y, label = signif(sum, 3)), 
                                      fontface = 2, size = 4)
    } # end if fac

  } else if (sumTF) {

    # summing the values within map domain
    p1 <- p1 +
      annotate('text', x = map.ext$maxlon - (map.ext$maxlon - map.ext$minlon) /10,
                       y = map.ext$maxlat - (map.ext$maxlat - map.ext$minlat) /10,
                       label = signif(sum(sel.foot$foot), 3), fontface = 2, size = 5)
  } # end if

  # plot observed XCO2, add footprint raster layer
  if (!is.null(oco2.path)) {
    if (qfTF) {
      p1 <- p1 + geom_point(data = qf.obs, aes(lon, lat, colour = xco2), size = 0.4)
      min.y <- floor(min(qf.obs$xco2, na.rm = T))
      max.y <- ceiling(max(qf.obs$xco2, na.rm = T))

    } else {
      p1 <- p1 + geom_point(data = obs, aes(lon, lat, colour = xco2), size = 0.4)
      min.y <- floor(min(obs$xco2, na.rm = T))
      max.y <- ceiling(max(obs$xco2, na.rm = T))
    }  # end if qfTF 
  }

  lab <- 10 ^ seq(-10, 2, 1)
  p2 <- p1 + 
    geom_raster(data = sel.foot, aes(lon + mm[[3]], lat + mm[[2]], fill = foot),
                alpha = 0.8) +
    scale_fill_gradientn(limits = c(min.foot.sig, max.foot.sig), 
                         name = 'FOOTPRINT\nppm/(umol/m2/s)', trans = 'log10', 
                         colours = col, breaks = lab, labels = lab) +
    scale_alpha_manual(values = c('all' = 0.5, 'screened' = 1.0))

  if ('fac' %in% colnames(sel.foot)){
    if (length(unique(sel.foot$fac)) > 1)
      p2 <- p2 + facet_wrap(~ fac, nrow = facet.nrow) +
                 theme(strip.text = element_text(size = font.size))
  }

  p3 <- p2 + theme(legend.position = leg.pos,
                   legend.text = element_text(size = font.size),
                   legend.key = element_blank(), 
                   legend.key.width = unit(width/10, 'cm'),
                   legend.key.height = unit(height/10, 'cm'),
                   axis.title.y = element_text(size = font.size, angle = 90),
                   axis.title.x = element_text(size = font.size, angle = 0),
                   axis.text = element_text(size = font.size),
                   axis.ticks = element_line(size = font.size),
                   title = element_text(size = font.size))

  if (!is.null(oco2.path)) {
    breaks <- seq(min.y, max.y, 2)
    p3 <- p3 + scale_colour_gradientn(name = 'OBS XCO2\n(ppm):', colours = col,
                                      limits = c(min.y, max.y), breaks = breaks, 
                                      labels = breaks) + 
              guides(colour = guide_colourbar(order = 2))
  }

  if (leg.pos == 'bottom') 
    p4 <- p3 + guides(fill = guide_legend(order = 1, nrow = 2, byrow = T))
  if (leg.pos == 'right') 
    p4 <- p3 + guides(fill = guide_legend(order = 1, ncol = 1))

  if (anthromesTF) {
    atm <- ggmap.anthromes(lon.lat, anthro.path, mm = NULL, site = site)
    urban.atm <- atm %>% filter(class < 20)  # for urban and dense settlements

    # convert data frame to shapefile and then find footprints in shapefile
    dfToRaster <- function(x) {
      coordinates(x) <- ~lon + lat
      gridded(x) <- TRUE
      projection(x) <- CRS("+proj=longlat +datum=WGS84")
      return(raster(x))
    }

    f.atm <- dfToRaster(urban.atm) %>% rasterToPolygons() %>% fortify()
    p4 <- p4 +
      geom_polygon(data = f.atm, aes(long, lat, group = group), fill = NA) +
      geom_path(data = f.atm, aes(long, lat, group = group), colour = 'darkorange', 
                linetype = 2, size = 1.0)
  } # end if anthromes

  print(picname)
  if (storeTF) 
    ggsave(p4, file = picname, width = width, height = height, dpi = 300)

  return(p4)
}
