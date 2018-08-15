#### subroutine to readin ODIAC emissions and then couple with STILT footprint
# as footprints have already been weighted by AK and PW,
# just multiple emission with 2D footprint map
# written by Dien Wu, 09/13/2016

# Updates:
# ADD ODIACv2016, flag 1 for using v2015a, flag 2 for v2016, DW 02/14/2017
# note that v2015a does not have emission for year 2015
# ADD PRD ODIAC emissions, DW 03/08/2017
# ADD TIMES hourly scaling factors for ODIACv2016, DW 03/08/2017
# Get rid of variable "odiac.vname",
# always preprocess and read ODIAC emission before call this function...

# version 2 modify based on Ben's code, DW
# can work with multiple receptors at a time, now, DW, 06/05/2018
# fix footprint lat/lon to lower lefts, as Ben uses centered lat/lon
# use raster rather than nc_open, DW, 06/19/2018
# add plotTF for plotting XCO2 contribution maps, DW, 06/20/2018
# store output contribution map into the same by-d directory,
#   remove store.path, DW, 07/26/2018
# remove foot.path, use full path as foot.file, DW, 07/26/2018

foot.odiacv3 <- function(foot.file, emiss.file, workdir, txtfile, lon.lat,
  plotTF = F){

  library(raster)

  # plot emissions
  if (plotTF) {
    emiss <- raster(emiss.file)
    emiss.df <- raster::as.data.frame(emiss, xy = T)
    colnames(emiss.df) <- list('lon', 'lat', 'emiss')
    emiss.df <- emiss.df %>% filter(emiss > 1)

    mm <- ggplot.map(map = 'ggmap', center.lat = lon.lat$citylat,
      center.lon = lon.lat$citylon, zoom = 8)

    # grab observations using map lat/lon
    map.ext <- c(min(mm[[1]]$data$lon), max(mm[[1]]$data$lon),
                 min(mm[[1]]$data$lat), max(mm[[1]]$data$lat))

    sel.emiss <- emiss.df %>% filter(
      lon >= map.ext[1] & lon <= map.ext[2] &
      lat >= map.ext[3] & lat <= map.ext[4])
    print(sel.emiss[sel.emiss$emiss >= 100, ])

    e1 <- mm[[1]] + coord_equal() +
      geom_raster(data = sel.emiss, aes(lon + mm[[3]], lat + mm[[2]],
        fill = emiss)) +
      scale_fill_gradientn(trans = 'log10', colours = def.col(),
        limits = c(1, 1E5))
    ggsave(plot = e1, filename = gsub('.nc', '.png', emiss.file),
      width = 8, hright = 8)
  }

  # if cannot find the correct format of nc file for emissions given selected area
  # and return ODIAC file name with path in front
  if (length(emiss.file) == 0) {
    cat('NO nc file found, check tif2nc.odiacv2() to create one...\n')
    return()

  } else {
    ## read in emissions
    emiss.dat <- raster(emiss.file)
    emiss.res <- res(emiss.dat)[1]
  }  # end if emiss.file

  # from foot.file, get receptor info
  receptor <- unlist(strsplit(gsub('_X_foot.nc', '', basename(foot.file)), "_"))
  receptor <- as.data.frame(matrix(receptor, byrow = T, ncol = 3),
    stringsAsFactors = F) %>% mutate_all(funs(as.numeric), colnames(receptor))
    # mutate_all() convert character to numberic
  colnames(receptor) <- list('timestr', 'lon', 'lat')

  order.index <- order(receptor$lat)
  receptor <- receptor[order.index, ]
  foot.file <- foot.file[order.index]

  receptor$xco2.ff <- NA

  # then loop over each receptor
  for (r in 1:nrow(receptor)) {

    # read in footprint
    foot.dat <- raster(file.path(foot.file[r]))
    crs(foot.dat) <- '+proj=longlat'
    #plot(log10(foot.dat))

    # NOW, foot and emiss should have the same dimension,
    # multiple them to get contribution map of CO2 enhancements
    #if (extent(emiss.dat) == extent(foot.dat)) {
    xco2.ff.sp <- raster::overlay(x = emiss.dat, y = foot.dat,
      fun = function(x, y){return(x * y)})    # spatial xco2.ff
    if (plotTF) plot(log10(xco2.ff.sp))

    #} else {
    #  cat(paste('Foot and emiss grid have different domain, please check...\n'))
    #  next
    #  if (r == nrow(receptor)) return()
    #} # end if check domain

    # sum the map to get the XCO2 enhancements,
    # note that AK and PW have been incorporated in footprint
    tmp.xco2.ff <- sum(getValues(xco2.ff.sp))
    receptor$xco2.ff[r] <- tmp.xco2.ff
    print(tmp.xco2.ff)

    ### store emission * column footprint = XCO2 contribution grid into .nc file
    # store into the same workding dir
    outfile <- file.path(dirname(foot.file[r]),
      gsub('foot', 'foot_emiss', basename(foot.file[r])))
    print(outfile)

    crs(xco2.ff.sp) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    writeRaster(xco2.ff.sp, outfile, overwrite = TRUE, format = "CDF",
      varname = "xco2", varunit = "PPM", xname = "lon", yname = "lat",
      longname = 'XCO2 enhancemnets due to ODIAC emission')
  }  # end for r

  # finally, write in a txt file
  write.table(x = receptor, file = txtfile, sep = ',', row.names = F, quote = F)

  return(receptor)
} # end of subroutine
