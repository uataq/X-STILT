# script to plot footprint with observed XCO2 on spatial maps,
# written by Dien Wu, 06/18/2018

# add the sum of foot in plotted region, DW, 07/18/2018

ggmap.xfoot.obs <- function(mm, lon.lat, site, oco2.path, facet.nrow = NULL,
  facet.ncol = NULL, nhrs.back, dpar, stilt.ver, foot.str = NULL, zisf = NULL,
  timestr, font.size = rel(0.9), recp.lon, recp.lat, foot, picname, storeTF,
  width = 12, height = 8){

  col <- def.col()
  m1 <- mm[[1]] + theme_bw() + coord_equal(1.0)

  # grab observations using map lat/lon
  map.ext <- c(min(mm[[1]]$data$lon), max(mm[[1]]$data$lon),
    min(mm[[1]]$data$lat), max(mm[[1]]$data$lat))

  cat('Reading OCO-2 data according to the spatial domain of ggmap...\n')
  obs <- grab.oco2(ocopath = oco2.path, timestr = timestr, lon.lat = map.ext)

  # select footprints using map.ext
  library(dplyr)
  sel.foot <- foot %>% filter(lon >= map.ext[1] & lon <= map.ext[2] &
    lat >= map.ext[3] & lat <= map.ext[4])
  sum.foot <- sel.foot %>% group_by(fac) %>% dplyr::summarize(sum = sum(foot))

  title <- paste0('Spatial time-integrated weighted column footprint (',
    nhrs.back, ' hours back; ', dpar, ' dpar; ', foot.str, '; ziscale = ', zisf,
    ')\nusing STILT version', stilt.ver, ' for overpass on ', timestr, ' for ',
    site, '\nSmall footprints < 1E-6 are displayed in gray')

  p1 <- m1 + labs(title = title, x = 'LONGITUDE [E]', y = 'LATITUDE [N]')

  # if there is 1+ receptors
  if (length(unique(sel.foot$fac)) > 1){

    # receptor locations and add receptors on map
    sel.recp <- data.frame(lon = recp.lon + mm[[3]], lat = recp.lat + mm[[2]],
      fac = unique(sel.foot$fac)) %>% full_join(sum.foot, by = 'fac')
    sel.recp$x <- map.ext[2] - 0.8
    sel.recp$y <- map.ext[4] - 0.5
    print(sel.recp)

    p1 <- p1 +
      geom_point(data = sel.recp, aes(lon, lat), colour = 'purple', size = 2) +
      facet_wrap(~ fac, scales = 'free', nrow = facet.nrow, ncol = facet.ncol) +
      geom_text(data = sel.recp, aes(lon + 0.4, lat),
        colour = 'purple', size = 4, label = 'receptor', fontface = 2) +
      facet_wrap(~ fac, scales = 'free', nrow = facet.nrow, ncol = facet.ncol) +
      geom_text(data = sel.recp, aes(x, y, label = signif(sum, 3)), fontface = 2,
        size = 4)
  }  # end if

  # plot observed XCO2, add footprint raster layer
  p2 <- p1 +
    geom_point(data = obs, aes(lon, lat, colour = xco2), size = 0.3) +
    geom_raster(data = sel.foot, aes(lon + mm[[3]], lat + mm[[2]], fill = foot),
      alpha = 0.8) +
    scale_fill_gradientn(limits = c(1E-6, 1E-2), name = 'Xfoot',
      trans = 'log10', colours = col, breaks = c(1E-6, 1E-4, 1E-2, 0.1, 1.0),
      labels = c(1E-6, 1E-4, 1E-2, 0.1, 1.0))

  if (length(unique(sel.foot$fac)) > 1){
    p2 <- p2 + facet_wrap(~ fac, scales = 'free', nrow = facet.nrow,
      ncol = facet.ncol)
  }

  p3 <- p2 + theme(legend.position = 'bottom',
    legend.text = element_text(size = font.size),
    legend.key = element_blank(), legend.key.width = unit(width/8, 'cm'),
    legend.key.height = unit(height/20, 'cm'),
    axis.title.y = element_text(size = font.size, angle = 90),
    axis.title.x = element_text(size = font.size, angle = 0),
    axis.text = element_text(size = font.size),
    axis.ticks = element_line(size = font.size),
    title = element_text(size = font.size),
    strip.text.x = element_text(size = font.size))

  max.y  <- ceiling(max(obs$xco2))
  min.y  <- floor(min(obs$xco2))
  breaks <- seq(min.y, max.y, 2)
  limits <- c(min.y, max.y)

  p4 <- p3 + scale_colour_gradientn(name = 'OBS XCO2:', colours = col,
    limits = limits, breaks = breaks, labels = breaks)
  p4 <- p4 + guides(colour = guide_colourbar(order = 2),
    fill = guide_legend(order = 1, nrow = 1))

  print(picname)
  if (storeTF) ggsave(p4, file = picname, width = width, height = height)

  return(p4)
}
