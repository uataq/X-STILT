# subroutine to read OCO-2 SIF, DW, 06/06/2018
# 'thred.count' for the sounding # thredshold

ggmap.obs.sif <- function(site, timestr, sif.path, lon.lat, workdir,
  plotdir = file.path(workdir, 'plot'), zoom = 9){

  # get Sif file name
  sif.file <- list.files(pattern = paste0('LtSIF_', substr(timestr, 3, 8)),
    path = sif.path)

  if (length(sif.file) == 0) {
    warnings('NO SIF file found for this OCO-2 overpass...'); next

  } else {  # if file found...

    library(ncdf4)
    sif.dat <- nc_open(file.path(sif.path, sif.file))
    time <- ncvar_get(sif.dat, "time")
    lat  <- ncvar_get(sif.dat, "latitude")
    lon  <- ncvar_get(sif.dat, "longitude")

    sif757 <- ncvar_get(sif.dat, "SIF_757nm") # unit in W/m2/sr/µm
    sif771 <- ncvar_get(sif.dat, "SIF_771nm")
    igbp   <- ncvar_get(sif.dat, "IGBP_index")

    sif <- data.frame(timestr = timestr, lat = as.numeric(lat),
      lon = as.numeric(lon), sif757 = as.numeric(sif757),
      sif771 = as.numeric(sif771), igbp = as.numeric(igbp))

    # select SIF in given region and scale SIF_771nm with sacling factor of
    # 1.35 (used in Luus et al., 2017) to calculate an averaged SIF
    sel.sif <- sif %>%
      filter(lon >= lon.lat[1] & lon <= lon.lat[2] &
             lat >= lon.lat[3] & lat <= lon.lat[4]) %>%
      mutate(avg.sif = (sif757 + sif771 * 1.35)/2)

    # assign months and seasons
    sel.sif <- sel.sif %>%  mutate(mon = substr(timestr, 5, 6),
        season =
          ifelse(mon %in% c('12', '01', '02'), "WINTER",
          ifelse(mon %in% c('03', '04', '05'), "SPRING",
          ifelse(mon %in% c('06', '07', '08'), "SUMMER", "FALL"))))

    # plot center
    obs.lat <- lon.lat[6]
    obs.lon <- lon.lat[5]
    alpha <- 1; font.size <- rel(1.1)
    col <- c('black', '#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8',
            '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')

    # plot google map
    sitemap <- get_map(location = c(lon = obs.lon, lat = obs.lat), zoom = zoom,
      maptype = 'roadmap')
    m1 <- ggmap(sitemap) + theme_bw()

    melt.sif <- sel.sif %>%
      dplyr::select('timestr', 'lat', 'lon', 'sif757', 'sif771', 'avg.sif')
    melt.sif <- melt(melt.sif, id.var = c('timestr', 'lat', 'lon'))

    title <- paste('OCO-2 SIF [W/m2/sr/µm] for', site, 'on', timestr)
    c1 <- m1 + labs(x = 'Longitude', y = 'Latitude', title = title) +
      geom_point(data = melt.sif, aes(lon, lat, colour = value), size = 1.0) +
      facet_wrap(~variable, ncol = 2) +
      scale_colour_gradientn(name = 'OCO-2 SIF', colours = col,
        limits = c(-1, max(2.5, max(melt.sif$sif))),
        breaks = seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5))

    # draw a rectangle around the city
    d <- data.frame(
      x = c(lon.lat[5] - 0.5, lon.lat[5] + 0.5, lon.lat[5] + 0.5, lon.lat[5] - 0.5),
      y = c(lon.lat[6] - 0.5, lon.lat[6] - 0.5, lon.lat[6] + 0.5, lon.lat[6] + 0.5))
    d <- rbind(d, d[1,])
    c2 <- c1 +
      geom_path(data = d, aes(x, y), linetype = 2, size = 1.2, colour = 'gray50')

    c2 <- c2 + theme(legend.position = 'bottom',
      legend.text = element_text(size = font.size),
      legend.key = element_blank(), legend.key.height = unit(0.5, 'cm'),
      legend.key.width = unit(3, 'cm'),
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

 }  # end if length()

 #return(sel.sif)
}
