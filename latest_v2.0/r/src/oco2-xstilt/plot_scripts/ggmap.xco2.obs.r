# script to plot XCO2 contribution with observed XCO2 on spatial maps,
# written by Dien Wu, 06/18/2018

ggmap.xco2.obs <- function(mm, lon.lat, site, facet.nrow, nhrs.back, dpar, sf,
  zisf, stilt.ver, timestr, font.size = rel(0.9), recp.lon, recp.lat, obs, xco2,
  picname, storeTF, width = 12, height = 8){

  col <- def.col()
  m1 <- mm[[1]] + theme_bw() + coord_equal(1.0)

  # grab observations using map lat/lon
  map.ext <- c(min(mm[[1]]$data$lon), max(mm[[1]]$data$lon),
    min(mm[[1]]$data$lat), max(mm[[1]]$data$lat))
  cat('Reading OCO-2 data according to the spatial domain of ggmap...\n')
  obs <- grab.oco2(ocopath = oco2.path, timestr, lon.lat = map.ext) # grab obs

  # select xco2 using map.ext
  library(dplyr)
  sel.xco2 <- xco2 %>% filter(
    lon >= map.ext[1] & lon <= map.ext[2] &
    lat >= map.ext[3] & lat <= map.ext[4])
  sum.xco2 <- sel.xco2 %>% group_by(fac) %>% dplyr::summarize(sum = sum(xco2))

  title <- paste('Spatial contribution of XCO2 [ppm] (',
    nhrs.back, ' hours back; dpar =', dpar, '; smooth factor =', sf,
    '; ziscale =', zisf, ')\nusing STILT version', stilt.ver, 'for overpass on',
    timestr, 'over', site,
    '\nSmall XCO2 enhancements < 1E-6 ppm are displayed in gray')

  # plot and label receptors and add title/xy axis
  p1 <- m1 + labs(title = title, x = 'LONGITUDE [E]', y = 'LATITUDE [N]')

  if (length(unique(xco2$fac)) > 1){

    # receptor locations and add receptors on map
    sel.recp <- data.frame(lon = recp.lon, lat = recp.lat,
      fac = unique(xco2$fac)) %>% full_join(sum.xco2, by = 'fac')
    sel.recp$x <- map.ext[2] - 0.8
    sel.recp$y <- map.ext[4] - 0.5
    print(sel.recp)

    p1 <- p1 +
      geom_point(data = sel.recp, aes(lon, lat), size = 1, colour = 'purple') +
      facet_wrap(~ fac, nrow = facet.nrow) +
      geom_text(data = sel.recp, aes(lon + 1, lat), size = 2, colour = 'purple',
        label = 'receptor', fontface = 2) +
      facet_wrap(~fac, nrow = facet.nrow) +
      geom_text(data = sel.recp, aes(x, y, label = signif(sum, 3)),
        fontface = 2, size = 4)
  }  # end if

  # plot observed XCO2, add xco2 raster layer
  p2 <- p1 +
    geom_point(data = obs, aes(lon, lat, colour = xco2),size = 0.3) +
    geom_raster(data = sel.xco2, aes(lon + mm[[3]], lat + mm[[2]], fill = xco2),
      alpha = 0.8) +
    scale_fill_gradientn(limits = c(1E-6, max(sel.xco2$xco2)),
      name = 'SIM XCO2 [ppm]', colours = col, trans = 'log10',
      breaks = c(1E-6, 1E-4, 1E-2, 0.1, 1.0),
      labels = c(1E-6, 1E-4, 1E-2, 0.1, 1.0))

  if (length(unique(sel.xco2$fac)) > 1){
    p2 <- p2 + facet_wrap(~fac, nrow = facet.nrow) +
      theme(strip.text = element_text(size = font.size))
  }

  p3 <- p2 + theme(legend.position = 'bottom',
    legend.text = element_text(size = font.size),
    legend.key = element_blank(), legend.key.width = unit(width/8, 'cm'),
    legend.key.height = unit(height/20, 'cm'),
    axis.title.y = element_text(size = font.size, angle = 90),
    axis.title.x = element_text(size = font.size, angle = 0),
    axis.text = element_text(size = font.size),
    axis.ticks = element_line(size = font.size),
    title = element_text(size = font.size))

  max.y  <- ceiling(max(obs$xco2)); min.y  <- floor(min(obs$xco2))
  breaks <- seq(min.y, max.y, 2); limits <- c(min.y, max.y)

  p4 <- p3 + scale_colour_gradientn(name = 'OBS XCO2:', colours = col,
    limits = limits, breaks = breaks, labels = breaks)
  p4 <- p4 + guides(colour = guide_colourbar(order = 2),
    fill = guide_legend(order = 1, nrow = 1))

  print(picname)
  if (storeTF) ggsave(p4, file = picname, width = width, height = height)

  return(p4)
}
