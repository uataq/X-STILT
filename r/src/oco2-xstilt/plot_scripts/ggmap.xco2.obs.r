# script to plot XCO2 contribution with observed XCO2 on spatial maps,
# written by Dien Wu, 06/18/2018

# last update for lite b9 data, DW, 10/19/2018

ggmap.xco2.obs <- function(mm, lon.lat, site, oco2.ver, oco2.path, facet.nrow, 
                           nhrs, dpar, foot.sf, zisf, met, stilt.ver, timestr, 
                           font.size = rel(0.9), recp.lon, recp.lat, xco2, 
                           min.xco2.sig = 1E-6, max.xco2.sig = 2, titleTF = T, 
                           sumTF = T, qfTF = T, storeTF = T, picname, width = 12, 
                           height = 8, leg.pos = c('bottom', 'right')[1]){

  col <- def.col()
  m1 <- mm[[1]] + theme_bw() + coord_equal(1.1)

  # grab observations using map lat/lon
  map.ext <- data.frame(minlon = min(mm[[1]]$data$lon),
                        maxlon = max(mm[[1]]$data$lon),
                        minlat = min(mm[[1]]$data$lat),
                        maxlat = max(mm[[1]]$data$lat))
  cat('Reading OCO-2 data according to the spatial domain of ggmap...\n')
  obs <- grab.oco2(ocopath = oco2.path, timestr, lon.lat = map.ext, oco2.ver)
  qf.obs <- obs %>% filter(qf == 0)

  # select xco2 using map.ext
  sel.xco2 <- xco2 %>% filter(lon >= map.ext$minlon & lon <= map.ext$maxlon &
                              lat >= map.ext$minlat & lat <= map.ext$maxlat &
                              xco2 >= min.xco2.sig)

  title <- paste0('Spatial contribution of XCO2 [ppm] (', nhrs, ' hours; dpar = ',
                  dpar, '; smooth factor = ', foot.sf, '; ziscale = ', zisf, 
                  '; met = ', met, ')\nusing STILT version', stilt.ver, 
                  ' for overpass on ', timestr,  ' over ', site, 
                  '\nOnly large XCO2 enhancements > ', min.xco2.sig, 
                  ' ppm are displayed')
  if (titleTF == F) title <- NULL

  # plot and label receptors and add title/xy axis
  p1 <- m1 + labs(title = title, x = 'LONGITUDE [E]', y = 'LATITUDE [N]')

  if ('fac' %in% colnames(sel.xco2)){
    if (length(unique(sel.xco2$fac)) > 1){

      # summing the values within map domain
      sum.xco2 <- sel.xco2 %>% group_by(fac) %>% dplyr::summarize(sum = sum(xco2))

      # receptor locations and add receptors on map
      sel.recp <- data.frame(lon = recp.lon, lat = recp.lat,
                             fac = unique(xco2$fac)) %>% 
                  full_join(sum.xco2, by = 'fac') %>%
                  mutate(x = map.ext$maxlon - (map.ext$maxlon - map.ext$minlon) / 10,
                         y = map.ext$maxlat - (map.ext$maxlon - map.ext$minlon) / 10)
      print(sel.recp)

      p1 <- p1 +
        geom_point(data = sel.recp, aes(lon, lat), size = 3, colour = 'purple',
                   shape = 1) +
        facet_wrap(~ fac, nrow = facet.nrow) +
        geom_text(data = sel.recp, aes(lon + 0.5, lat), size = 4,
                  colour = 'purple', label = 'receptor', fontface = 1)

      if (sumTF) p1 <- p1 + geom_text(data = sel.recp,
          aes(x, y, label = signif(sum, 3)), fontface = 2, size = 4)
    }  # end if
  }  # end if

  # plot observed XCO2, add xco2 raster layer
  lab <- 10 ^ seq(-10, 2, 1)

  if (qfTF) {
    p2 <- p1 + geom_point(data = qf.obs, aes(lon, lat, colour = xco2), size = 0.4)
    min.y <- floor(min(qf.obs$xco2, na.rm = T))
    max.y <- ceiling(max(qf.obs$xco2, na.rm = T))

  } else {
    p2 <- p1 + geom_point(data = obs, aes(lon, lat, colour = xco2), size = 0.4)
    min.y <- floor(min(obs$xco2, na.rm = T))
    max.y <- ceiling(max(obs$xco2, na.rm = T))
  }  # end if qfTF

  p2 <- p2 + 
    geom_raster(data = sel.xco2, aes(lon + mm[[3]], lat + mm[[2]], fill = xco2),
                alpha = 0.8) +
    scale_fill_gradientn(limits = c(min.xco2.sig, max.xco2.sig), 
                         name = 'SIM\nXCO2 [ppm]', colours = col, 
                         trans = 'log10', breaks = lab, labels = lab)

  if ('fac' %in% colnames(sel.xco2))
    if (length(unique(sel.xco2$fac)) > 1)
      p2 <- p2 + facet_wrap(~fac, nrow = facet.nrow) +
                 theme(strip.text = element_text(size = font.size))

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

  breaks <- seq(min.y, max.y, 2); limits <- c(min.y, max.y)
  p4 <- p3 + scale_colour_gradientn(name = 'OBS\nXCO2 [ppm]:', colours = col,
                                    limits = limits, breaks = breaks, 
                                    labels = breaks) +
             guides(colour = guide_colourbar(order = 2))
                              

  if (leg.pos == 'bottom') 
    p4 <- p4 + guides(fill = guide_legend(order = 1, nrow = 2, byrow = T))
  if (leg.pos == 'right') 
    p4 <- p4 + guides(fill = guide_legend(order = 1, ncol = 1))

  print(picname)
  if (storeTF) ggsave(p4, file = picname, width = width, height = height)

  return(p4)
}
