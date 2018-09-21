# script to plot observed XCO2 on map, DW
# add section to plot observed XCO2 with quality flag on latitude series, 07/03/2018
# use grab.oco2 to read observations, DW, 08/07/2018

ggmap.obs.xco2 <- function(site, timestr, oco2.path, lon.lat, workdir,
  plotdir = file.path(workdir, 'plot'), zoom = 8){

  library(ggmap); library(ggplot2)
  obs.all <- grab.oco2(ocopath = oco2.path, timestr, lon.lat)

  # plot google map
  alpha <- 1; font.size <- rel(1.2); col <- def.col()

  m1 <- ggplot.map(map = 'ggmap', zoom = zoom, center.lat = lon.lat$citylat,
    center.lon = lon.lat$citylon)[[1]] + theme_bw()
  c1 <- m1 + geom_point(data = obs.all, aes(lon, lat, colour = xco2)) +
    scale_colour_gradientn(name = 'OCO-2 XCO2 [ppm]', colours = col,
      limits = c(max(390, min(obs.all$xco2, na.rm = T)),
        max(obs.all$xco2, na.rm = T)), breaks = seq(380, 420, 2),
      labels = seq(380, 420, 2)) +
    labs(x = 'Longitude', y = 'Latitude') +
    labs(title = paste('OCO-2 XCO2 [ppm] for', site, 'on', timestr))

  # draw a rectangle around the city
  d <- data.frame(
    x = c(lon.lat$citylon - 0.5, rep(lon.lat$citylon + 0.5, 2), lon.lat$citylon - 0.5),
    y = c(rep(lon.lat$citylat - 0.5, 2), rep(lon.lat$citylat + 0.5, 2)))

  d <- rbind(d, d[1,])
  c2 <- c1 + geom_path(data = d, aes(x, y), linetype = 2, size = 1.2,
    colour = 'gray50')

  c2 <- c2 + theme(legend.position = 'bottom',
    legend.text = element_text(size = font.size),
    legend.key = element_blank(), legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(3, 'cm'),
    axis.title.y = element_text(size = font.size, angle = 90),
    axis.title.x = element_text(size = font.size, angle = 0),
    axis.text = element_text(size = font.size),
    axis.ticks = element_line(size = font.size),
    title = element_text(size = font.size))

  # plot on latitude-series
  l1 <- ggplot() + theme_bw() +
    geom_point(data = obs.all, aes(lat, xco2), colour = 'gray70', shape = 17,
      size = 3) +
    geom_point(data = obs.all[obs.all$qf == 0, ], aes(lat, xco2), size = 3,
      colour = 'black', shape = 17)

  # merge map and latseries
  library(ggpubr)
  merge.plot <- ggarrange(plotlist = list(c2, l1), nrow = 2, heights = c(2, 1))

  picname <- paste0('ggmap_xco2_', site, '_', timestr, '.png')
  picfile <- file.path(plotdir, picname)
  print(picfile)
  ggsave(merge.plot, filename = picfile, width = 11, height = 15)
}
