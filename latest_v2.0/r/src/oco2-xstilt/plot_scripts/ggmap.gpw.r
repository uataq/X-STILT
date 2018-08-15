# subroutine to read gridded population data, DW, 07/07/2018

ggmap.gpw <- function(lon.lat, pop.path, mm = NULL, site = NULL,
  picpath = NULL, font.size = rel(1.0), width = 10, height = 11){

  # read data
  library(raster)
  pop.file <- list.files(path = pop.path, pattern = '.tif')
  pop.dat <- raster(file.path(pop.path, pop.file))

  # crop population data
  sel.pop.dat <- crop(pop.dat,
    extent(lon.lat$minlon, lon.lat$maxlon, lon.lat$minlat, lon.lat$maxlat))
  pop.df <- raster::as.data.frame(sel.pop.dat, xy = T)
  colnames(pop.df) <- list('lon', 'lat', 'pop')

  # plot
  if (!is.null(mm)) {

    g1 <- mm[[1]] + theme_bw() + coord_equal(1.1) +
      geom_raster(data = pop.df, aes(lon + mm[[3]], lat + mm[[2]], fill = pop),
      alpha = 0.7)

    col <- def.col()
    g2 <- g1 + labs(x = "LONGITUDE", y = "LATITUDE") +
      scale_fill_gradientn(name = 'Population', trans = 'log10',
      colors = col, breaks = c(10, 1E2, 1E3, 1E4, 1E5))

    g3 <- g2 + theme(legend.position = 'bottom',
      legend.text = element_text(size = font.size),
      legend.key = element_blank(), legend.key.width = unit(2, 'cm'),
      legend.key.height = unit(0.6, 'cm'),
      axis.title.y = element_text(size = font.size, angle = 90),
      axis.title.x = element_text(size = font.size, angle = 0),
      axis.text = element_text(size = font.size),
      axis.ticks = element_line(size = font.size),
      title = element_text(size = font.size),
      strip.text.x = element_text(size = font.size))

    picname <- paste0('population_', site, '.png')
    picname <- file.path(picpath, picname)
    print(picname)
    ggsave(g3, file = picname, width = width, height = height)

  } else {
    g3 <- NULL
  } # end if mm

  return(list(pop.df, g3))
}
