# subroutine to read OCO-2 SIF, DW, 06/06/2018
# 'thred.count' for the sounding # thredshold

ggmap.obs.sif <- function(site, timestr, sif.path, lon.lat, workdir,
                          plotdir = file.path(workdir, 'plot'), zoom = 8, 
                          box.dlon = 0.5, box.dlat = 0.5){

  sel.sif <- grab.sif(sif.path, timestr, lon.lat)
  if (length(sel.sif) == 0) {cat('NO SIF file matched, returning NA...\n'); return()}
  if (nrow(sel.sif) == 0) {cat('NO SIF retrieved within lon.lat, returning NA...\n'); return()}

  # plot center
  alpha <- 1; font.size <- rel(0.9); col <- def.col()

  # plot google map
  m1 <- ggplot.map(map = 'ggmap', zoom = zoom, center.lat = lon.lat$citylat,
                   center.lon = lon.lat$citylon)[[1]] + theme_bw()

  melt.sif <- sel.sif %>%
    dplyr::select('timestr', 'lat', 'lon', 'sif757', 'sif771', 'avg.sif')
  melt.sif <- melt(melt.sif, id.var = c('timestr', 'lat', 'lon'))

  title <- paste('OCO-2 SIF [W/m2/sr/µm] for', site, 'on', timestr)
  c1 <- m1 + labs(x = 'Longitude', y = 'Latitude', title = title) +
             geom_point(data = melt.sif, aes(lon, lat, colour = value), size = 0.9) +
             facet_wrap(~variable, ncol = 2) +
             scale_colour_gradientn(name = 'OCO-2 SIF', colours = col,
                                    limits = c(-1, max(2.5, max(melt.sif$value))),
                                    breaks = seq(-4, 4, 0.5), 
                                    labels = seq(-4, 4, 0.5))

  # draw a rectangle around the city
  d <- data.frame(x = c(lon.lat$citylon - box.dlon, 
                        rep(lon.lat$citylon + box.dlon, 2), 
                        lon.lat$citylon - box.dlon),
                  y = c(rep(lon.lat$citylat - box.dlat, 2), 
                        rep(lon.lat$citylat + box.dlat, 2)))

  d <- rbind(d, d[1,])
  c2 <- c1 + geom_path(data = d, aes(x, y), linetype = 2, size = 0.9,
                       colour = 'gray50') +
    theme(legend.position = 'bottom',
          legend.text = element_text(size = font.size),
          legend.key = element_blank(), legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(2, 'cm'),
          axis.title.y = element_text(size = font.size, angle = 90),
          axis.title.x = element_text(size = font.size, angle = 0),
          axis.text = element_text(size = font.size),
          axis.ticks = element_line(size = font.size),
          title = element_text(size = font.size),
          strip.text = element_text(size = font.size))

  picname <- paste0('ggmap_sif_', site, '_', timestr, '.png')
  picfile <- file.path(plotdir, picname)
  print(picfile)
  ggsave(c2, filename = picfile, width = 12, height = 13)

}
