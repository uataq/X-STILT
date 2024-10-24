# script to plot footprint with observed XCO2 on spatial maps,
# written by Dien Wu, 06/18/2018

# add the sum of foot in plotted region, DW, 07/18/2018
# last update, DW, 08/22/2018
# if no OCO-2 path, do not plot observed XCO2, DW, 10/29/2018

# mm can be generated by ggplot.map()
ggmap.xfoot.obs <- function(mm, site, oco.ver = NULL, oco.path = NULL, timestr, 
                            recp.lon = NULL, recp.lat = NULL, foot, 
                            min.foot.sig = 1E-8, max.foot.sig = 1E-2, qfTF = T, 
                            dxco2 = 2, title = NULL, picname, storeTF = T, 
                            width = 9, height = 7, leg.pos = 'bottom', 
                            scale.coord = 1.2, facet.nrow = 1, 
                            foot.unit = 'ppm/(umol/m2/s)', 
                            font.size = rel(0.8)){
  
  library(ggpubr)
  col <- def.col()
  m1 <- mm[[1]] + theme_bw() + coord_equal(scale.coord)

  # grab observations using map lat/lon
  map.ext <- data.frame(minlon = min(m1$data$lon), maxlon = max(m1$data$lon),
                        minlat = min(m1$data$lat), maxlat = max(m1$data$lat))

  if (!is.null(oco.path)) {
    cat('Reading satellite data according to the spatial domain of ggmap...\n')
    obs <- grab.oco(oco.path, timestr, lon.lat = map.ext)
    if (qfTF) obs <- obs %>% filter(qf == 0)
  }

  # select footprints using map.ext
  sel.foot <- foot %>% filter(lon >= map.ext$minlon & lon <= map.ext$maxlon &
                              lat >= map.ext$minlat & lat <= map.ext$maxlat &
                              foot >= min.foot.sig)
  p1 <- m1 + labs(x = 'LONGITUDE [E]', y = 'LATITUDE [N]')

  # if there is 1+ receptors
  if ('fac' %in% colnames(sel.foot)){
    if (length(unique(sel.foot$fac)) > 1){

      # summing the values within map domain
      sum.foot <- sel.foot %>% group_by(fac) %>% dplyr::summarize(sum = sum(foot))

      # receptor locations and add receptors on map
      sel.recp <- data.frame(lon = recp.lon, lat = recp.lat, fac = unique(sel.foot$fac)) %>% 
                  full_join(sum.foot, by = 'fac') %>%
                  mutate(x = map.ext$maxlon - 0.5, y = map.ext$maxlat - 0.5)
      print(sel.recp)

      p1 <- p1 + geom_point(data = sel.recp, aes(lon, lat), colour = 'purple', 
                            size = 3, shape = 1) +
                 facet_wrap(~ fac, nrow = facet.nrow) +
                 geom_text(data = sel.recp, aes(lon + 0.5, lat), 
                           colour = 'purple', size = 4, label = 'receptor', 
                           fontface = 1) 
    } # end if fac
  } 

  lab <- 10 ^ seq(-20, 3, 1)
  p2 <- p1 + geom_tile(data = sel.foot, 
                       aes(lon + mm[[3]], lat + mm[[2]], fill = foot),
                       alpha = 0.8) +
             scale_fill_gradientn(limits = c(min.foot.sig, max.foot.sig), 
                                  name = paste0('FOOTPRINT\n', foot.unit), 
                                  trans = 'log10', colours = col, 
                                  breaks = lab,labels = lab)

  if ('fac' %in% colnames(sel.foot))
    if (length(unique(sel.foot$fac)) > 1)
      p2 <- p2 + facet_wrap(~ fac, nrow = facet.nrow) +
                 theme(strip.text = element_text(size = font.size))

  # draw a rectangle around the site_
  #lon.lat <- get.lon.lat(site, dlon = 0.25, dlat = 0.25)
  #d <- data.frame(x = c(lon.lat$minlon, rep(lon.lat$maxlon, 2), lon.lat$minlon),
  #                y = c(rep(lon.lat$minlat, 2), rep(lon.lat$maxlat, 2)))
  #d <- rbind(d, d[1,])
  p3 <- p2 + #geom_path(data = d, aes(x, y), size = 0.4, colour = 'gray30') +
             theme(legend.position = leg.pos, legend.key = element_blank(), 
                   legend.text = element_text(size = font.size),
                   legend.key.width = unit(width/8, 'cm'),
                   legend.key.height = unit(height/14, 'cm'),
                   axis.title.y = element_text(size = font.size, angle = 90),
                   axis.title.x = element_text(size = font.size, angle = 0),
                   axis.text = element_text(size = font.size),
                   axis.ticks = element_line(size = font.size),
                   title = element_text(size = font.size))

  if (!is.null(oco.path)) {
    p3 <- p3 + geom_point(data = obs, aes(lon, lat, colour = xco2), size = 0.4)
    min.y <- floor(min(obs$xco2, na.rm = T))
    max.y <- ceiling(max(obs$xco2, na.rm = T))
    #min.y = 410; max.y = 424
    breaks <- seq(min.y, max.y, dxco2)
    p3 <- p3 + scale_colour_gradientn(name = 'OBS XCO2\n(ppm):', colours = col,
                                      limits = c(min.y, max.y), breaks = breaks, 
                                      labels = breaks) + 
               guides(colour = guide_colourbar(order = 2))
  }

  #if (leg.pos %in% c('bottom', 'top')) 
  #  p4 <- p3 + guides(fill = guide_legend(order = 1, nrow = 1))
  #if (leg.pos %in% c('left', 'right')) 
  #  p4 <- p3 + guides(fill = guide_legend(order = 1, ncol = 1))
  if (!is.null(title)) p3 <- annotate_figure(p3, top = title)
  if (storeTF) ggsave(p3, file = picname, width = width, height = height, dpi = 300)
 
  return(p3)
}









if (F) {

  l1 <- ggplot() + theme_bw() + labs(x = 'XCO2', y = 'LATITUDE') + 
        geom_point(data = obs, aes(xco2, lat), size = 1.6, shape = 17) + 
        ylim(c(map.ext$minlat + 0.1, map.ext$maxlat - 0.1))
        
  if (length(recp.lat) * length(recp.lon) == 1) {
    recp.obs <- obs %>% filter(abs(lat - recp.lat) < 1E-4, 
                               abs(lon - recp.lon) < 1E-4)
    l1 <- l1 + geom_point(data = recp.obs, aes(xco2, lat), size = 2, 
                          shape = 17, color = 'orange')
  } 

  #p1.cy <- ggplot_gtable(ggplot_build(p4))
  #l1.cy <- ggplot_gtable(ggplot_build(l1))
  #l1.cy$heights <- p1.cy$heights
  #pp <- gridExtra::grid.arrange(p1.cy, l1.cy, ncol = 2, widths = c(3, 1))
  
  pp <- ggarrange(p4, l1, ncol = 2, widths = c(2.5, 1))#, common.legend = T)
  if (leg.pos %in% c('bottom', 'top'))
    pp <- ggarrange(p4, l1, ncol = 2, widths = c(2.5, 1), common.legend = T)

  #if (!is.null(title)) pp <- annotate_figure(pp, top = title)
}
