# script to plot observed XCO2 on map, DW
# add section to plot observed XCO2 with quality flag on latitude series, 07/03/2018

ggmap.obs.xco2 <- function(site, timestr, oco2.path, lon.lat, workdir,
  plotdir = file.path(workdir, 'plot'), zoom = 8){

  library(ggmap); library(ggplot2)

  oco2.file <- list.files(pattern = substr(timestr, 3, 8), path = oco2.path)
  oco2.dat <- nc_open(file.path(oco2.path, oco2.file))

  ## grabbing OCO-2 levels, lat, lon
  # level 1 to 20, for space-to-surface, level 20 is the bottom level
  oco2.level <- ncvar_get(oco2.dat, "levels")
  oco2.lat <- ncvar_get(oco2.dat, "latitude")
  oco2.lon <- ncvar_get(oco2.dat, "longitude")
  xco2.obs <- ncvar_get(oco2.dat, "xco2")
  xco2.obs[xco2.obs == -999999] <- NA
  xco2.obs.uncert <- ncvar_get(oco2.dat, "xco2_uncertainty")

  wl <- ncvar_get(oco2.dat, "warn_level")
  qf <- ncvar_get(oco2.dat, "xco2_quality_flag")
  foot<-ncvar_get(oco2.dat, "Sounding/footprint")

  # YYYY MM DD HH mm ss m (millisecond) f (footprint)
  id <- as.character(ncvar_get(oco2.dat, "sounding_id"))
  sec <- ncvar_get(oco2.dat, "time")
  time <- as.POSIXct(sec, origin="1970-01-01 00:00:00", tz="UTC")

  # select regions, lon.lat: c(minlat, maxlat, minlon, maxlon)
  region.index <- oco2.lon >= lon.lat[1] & oco2.lon <= lon.lat[2] &
                  oco2.lat >= lon.lat[3] & oco2.lat <= lon.lat[4]

  # get OCO-2 overpassing time
  sel.id   <- as.numeric(id[region.index])
  sel.time <- as.character(time[region.index])
  sel.hr   <- as.numeric(substr(sel.time, 12, 13))
  sel.min  <- as.numeric(substr(sel.time, 15, 16))
  sel.foot <- as.numeric(foot[region.index])
  sel.wl   <- as.numeric(wl[region.index])
  sel.qf   <- as.numeric(qf[region.index])
  sel.lat  <- as.numeric(oco2.lat[region.index])
  sel.lon  <- as.numeric(oco2.lon[region.index])
  sel.xco2.obs <- as.numeric(xco2.obs[region.index])
  sel.xco2.obs.uncert <- as.numeric(xco2.obs.uncert[region.index])

  obs.all <- data.frame(id = sel.id, time = sel.time, hr = sel.hr, min = sel.min,
    lat = sel.lat, lon = sel.lon, xco2 = sel.xco2.obs, foot = sel.foot,
    xco2.uncert = sel.xco2.obs.uncert, wl = sel.wl, qf = sel.qf)

  # plot center
  obs.lat <- lon.lat[6]
  obs.lon <- lon.lat[5]
  alpha <- 1; font.size <- rel(1.2)
  col <- c('black', '#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8',
          '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')

  # plot google map
  sitemap <- get_map(location = c(lon = obs.lon, lat = obs.lat), zoom = zoom,
    maptype = 'roadmap')
  m1 <- ggmap(sitemap) + theme_bw()
  c1 <- m1 + geom_point(data = obs.all, aes(lon, lat, colour = xco2)) +
    scale_colour_gradientn(name = 'OCO-2 XCO2 [ppm]', colours = col,
      limits = c(max(390, min(obs.all$xco2, na.rm = T)),
        max(obs.all$xco2, na.rm = T)), breaks = seq(380, 420, 2),
      labels = seq(380, 420, 2)) +
    labs(x = 'Longitude', y = 'Latitude') +
    labs(title = paste('OCO-2 XCO2 [ppm] for', site, 'on', timestr))

  # draw a rectangle around the city
  d <- data.frame(
    x = c(lon.lat[5] - 0.5, lon.lat[5] + 0.5, lon.lat[5] + 0.5, lon.lat[5] - 0.5),
    y = c(lon.lat[6] - 0.5, lon.lat[6] - 0.5, lon.lat[6] + 0.5, lon.lat[6] + 0.5))
  d <- rbind(d, d[1,])
  c2 <- c1 + geom_path(data = d, aes(x, y), linetype = 2, size = 1.2, colour = 'gray50')

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
