# subroutine to get emission errors, based on the spread and mean of 3 emission datasets
# DW

# odiac.file.2008 for the ODIAC emission (specifically for year 2008) to carry out emission error
# ffdas.file for FFDAS emissions (year 2008)
# edgar.file for EDGAR emissions (year 2008)
# emiss.file for the ODIAC emission (during overpass month, e.g., 2014 DEC emission)

cal.emiss.err <- function (site, timestr, odiac.file.2008, edgar.file = NA, 
                           ffdas.file = NA, emiss.file, overwrite, plotTF = F) {

  library(raster)
  odiac.2008  <- raster(odiac.file.2008)
  crop.extent <- extent(odiac.2008)
  if (crop.extent[3] < 1E-5 & crop.extent[3] > 0) crop.extent[3] <- 0
  if (crop.extent[1] < 1E-5 & crop.extent[1] > 0) crop.extent[1] <- 0

  # aggregate 1kmx1km odiac to 0.1degree (fac = 12) to match resolution of others
  cat('cal.emiss.err(): Aggregating 1km ODIAC to 0.1deg...\n')
  sel.odiac.2008 <- aggregate(odiac.2008, fact = 12)

  # load EDGAR and convert from kg/m2/s to umol/m2/s
  if (is.na(edgar.file)) stop('cal.emiss.err(): Missing EDGAR file, stop\n')
  edgar <- raster(edgar.file) * 1E3 / 44 * 1E6

  # fix (0, 360) to (-180, 180) for EDGAR's longitude
  if (extent(edgar)[2] > 181) {
    east <- crop(edgar, extent(0, 180, -90, 90))
    west <- crop(edgar, extent(180, 360, -90, 90))

    # then change extent of west to negative long
    extent(west) <- c(-180, 0, -90, 90)
    edgar <- merge(west, east)
  }
  sel.edgar <- crop(edgar, crop.extent)  # select EDGAR emissions

  # load FFDAS and convert unit from kgC/m2/y to umol/m2/s
  if (is.na(ffdas.file)) cat('cal.emiss.err(): Missing FFDAS file\n')
  ffdas <- raster(ffdas.file) * 1E3 / 12 / 366 / 24 / 3600 * 1E6
  ffdas[ffdas < 0] <- 0
  sel.ffdas <- crop(ffdas, crop.extent)

  # calculate the spread and mean of three emission grids
  stack.emiss <- stack(sel.odiac.2008, sel.edgar, sel.ffdas)
  sd.emiss <- calc(stack.emiss, sd) # in ppm
  mean.emiss <- calc(stack.emiss, mean) # in ppm
  frac.emiss.err <- sd.emiss / mean.emiss  # fractional uncertainty

  # !!!! for absolute emission error, 
  # remember to use the CURRENT emission for the OVERPASS MONTH, NOT for year 2008
  emiss.current <- raster(emiss.file)
  emiss.current.agg <- aggregate(emiss.current, fact = 12)
  abs.emiss.err <- frac.emiss.err * emiss.current.agg
  #plot(log10(abs.emiss.err))

  filename <- file.path(dirname(emiss.file), paste0(site, '_emiss_err.nc'))
  if (overwrite)
    writeRaster(x = abs.emiss.err, filename = filename, overwrite = overwrite)

  # plotting 
  if (plotTF) {
    df <- raster::as.data.frame(frac.emiss.err, xy = T)
    colnames(df) <- list('lon', 'lat', 'abs.uncert')
    df2 <- raster::as.data.frame(emiss.current.agg, xy = T)
    colnames(df2) <- list('lon', 'lat', 'emiss')
    df <- df %>% left_join(df2, by = c('lon', 'lat'))
    # select fractional uncertainty over large ODIAC

    m1 <- ggplot.map(map = 'black', minlon = crop.extent[1], 
                     maxlon = crop.extent[2], minlat = crop.extent[3], 
                     maxlat = crop.extent[4])
    e1 <- m1 + 
      geom_raster(data = df[df$emiss > 1, ], aes(lon, lat, fill = abs.uncert * 100)) + 
      scale_fill_gradient(low = 'yellow', high = 'red', 
                          name = 'Fractional emission error [%]')
  } # end if plotTF

  return(filename)
}
# end of subroutine
