# script to plot observed XCO2 on map, DW
# add section to plot observed XCO2 with quality flag on latitude series, 07/03/2018
# use grab.oco2 to read observations, DW, 08/07/2018
# update for v9 lite data, QF filtering controled by qfTF, DW, 10/19/2018

ggmap.obs.xco2 <- function(site, timestr, oco2.ver, oco2.path, lon.lat, workdir,
                           plotdir = file.path(workdir, 'plot'), zoom = 8, 
                           qfTF = F, box.dlat = 0.5, box.dlon = 0.5){

  library(ggmap); library(ggplot2)
  obs.all <- grab.oco2(ocopath = oco2.path, timestr, lon.lat, oco2.ver)

  # plot google map
  alpha <- 1; font.size <- rel(0.9); col <- def.col()
  m1 <- ggplot.map(map = 'ggmap', zoom = zoom, center.lat = lon.lat$citylat,
                   center.lon = lon.lat$citylon)[[1]]

  if (qfTF) {
    qf.obs <- obs.all %>% filter(qf == 0)
    c1 <- m1 + geom_point(data = qf.obs, aes(lon, lat, colour = xco2), size = 0.6)
    min.y <- min(qf.obs$xco2, na.rm = T)
    max.y <- max(qf.obs$xco2, na.rm = T)

  } else {
    c1 <- m1 + geom_point(data = obs.all, aes(lon, lat, colour = xco2), size = 0.6)
    min.y <- min(obs.all$xco2, na.rm = T)
    max.y <- max(obs.all$xco2, na.rm = T)
  }  # end if qfTF

  c1 <- c1 + theme_bw() + 
    scale_colour_gradientn(name = 'OCO-2 XCO2 [ppm]', colours = col,
                           limits = c(max(390, min.y), max.y), 
                           breaks = seq(380, 420, 2), labels = seq(380, 420, 2)) +
    labs(x = 'LONGITUDE', y = 'LATITUDE') +
    labs(title = paste('OCO-2 XCO2 [ppm] for', site, 'on', timestr))

  # draw a rectangle around the city
  d <- data.frame(x = c(lon.lat$citylon - box.dlon, 
                        rep(lon.lat$citylon + box.dlon, 2), 
                        lon.lat$citylon - box.dlon),
                  y = c(rep(lon.lat$citylat - box.dlat, 2), 
                        rep(lon.lat$citylat + box.dlat, 2)))

  d <- rbind(d, d[1,])
  c2 <- c1 + geom_path(data = d, aes(x, y), linetype = 2, size = 0.9,
                       colour = 'gray50')

  c2 <- c2 + theme(legend.position = 'bottom',
                   legend.text = element_text(size = font.size),
                   legend.key = element_blank(), 
                   legend.key.height = unit(0.5, 'cm'),
                   legend.key.width = unit(2, 'cm'),
                   axis.title.y = element_text(size = font.size, angle = 90),
                   axis.title.x = element_text(size = font.size, angle = 0),
                   axis.text = element_text(size = font.size),
                   axis.ticks = element_line(size = font.size),
                   title = element_text(size = font.size))

  # plot on latitude-series
  l1 <- ggplot() + theme_bw() + labs(y = 'OCO-2 XCO2 [ppm]', x = 'LATITUDE')
  if (qfTF) {
    l1 <- l1 + geom_point(data = qf.obs, aes(lat, xco2), size = 3, 
                          colour = 'black', shape = 17)
  } else {
    l1 <- l1 + geom_point(data = obs.all, aes(lat, xco2), colour = 'gray70', 
                          shape = 17, size = 3) + 
               geom_point(data = qf.obs, aes(lat, xco2), size = 3, 
                          colour = 'black', shape = 17)
  }  # end if qfTF

  # merge map and latseries
  library(ggpubr)
  merge.plot <- ggarrange(plotlist = list(c2, l1), nrow = 2, heights = c(2, 1))

  picname <- paste0('ggmap_xco2_', site, '_', timestr, '.png')
  picfile <- file.path(plotdir, picname)
  print(picfile)
  ggsave(merge.plot, filename = picfile, width = 11, height = 15)
}
