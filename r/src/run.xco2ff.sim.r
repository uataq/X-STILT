#### subroutine to readin ODIAC emissions and then couple with STILT footprint
#' as footprints have already been weighted by AK and PW,
#' just multiple emission with 2D footprint map, originated from 'foot.odiacv3'
#' @author Dien Wu, 09/13/2016

#' @Updates:
#' ADD ODIACv2016, flag 1 for using v2015a, flag 2 for v2016, DW 02/14/2017
#' note that v2015a does not have emission for year 2015
#' ADD PRD ODIAC emissions, DW 03/08/2017
#' ADD TIMES hourly scaling factors for ODIACv2016, DW 03/08/2017
#' Get rid of variable 'odiac.vname',
#'       always preprocess and read ODIAC emission before call this function...

#' version 2 modify based on Ben's code, DW
#'       can work with multiple receptors at a time, now, DW, 06/05/2018
#' fix footprint lat/lon to lower lefts, as Ben uses centered lat/lon
# '      use raster rather than nc_open, DW, 06/19/2018
#'
#' add plotTF for plotting XCO2 contribution maps, DW, 06/20/2018
#' store output contribution map into the same by-d directory,
# '      remove store.path, DW, 07/26/2018
#' 
#' remove foot.path, use full path as foot.file, DW, 07/26/2018
#' fix a minor bug in interpreting footprint filename, DW, 10/11/2018 
#' accommodate the diff in foot.nc filename, using two versions, DW, 01/25/2019
# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# simplify txtfile name, DW, 06/28/2020 

run.xco2ff.sim <- function(site = 'Riyadh', timestr = '2014100920', vname = '2018', 
                           tiff.path, outdir, foot.res, workdir, store.path, nhrs, 
                           oco.sensor, oco.ver, met, lon.lat, run_emiss_err, 
                           edgar.file = NA, ffdas.file = NA, plotTF = F, writeTF = T){

  # grab footprint files and get footprint domain
  foot.path <- file.path(outdir, 'by-id')
  #foot.patt <- paste0(signif(foot.res, 3), 'x', signif(foot.res, 3), '_foot.nc')
  foot.patt <- 'X_foot.nc'
  foot.files <- list.files(foot.path, foot.patt, recursive = T, full.names = T)
  foot.indx <- grep(signif(foot.res, 3), basename(foot.files))

  # in a previous version of X-STILT, foot.res is shown on the foot nc file, 
  # but, after the refactoring, no more foot.res, 
  # make this if statement, if foot was generated using a previous version
  # DW, 01/25/2019
  if (length(foot.indx) > 0) {
    foot.file <- foot.files[foot.indx]
  } else foot.file <- foot.files
  
  if (length(foot.file) == 0) {
    stop('run.xco2ff.sim(): NO footprint found...please check by-id...\n')
    return()
  }

  tmp.foot  <- stack(foot.file[1])
  foot.extent <- extent(tmp.foot)

  # call tif2nc.odiacv2() to subset and get emiss file name
  # get cropped ODIAC emission for the overpass month
  # moved from main script to this subroutine, DW, 10/21/2018
  emiss.file <- tif2nc.odiacv3(site, timestr, vname, workdir, foot.extent,
                               tiff.path, gzTF = F)

  # txt file name for outputting model results
  txtfile <- file.path(store.path, paste0(timestr, '_', site, '_XCO2ff_', abs(nhrs), 
                                          'hrs_', oco.sensor, oco.ver, '_', met, 
                                          '_odiac', vname, '.txt'))
  
  # add emission error file with absolute emission uncertainty and txtfile
  if (run_emiss_err) {

    # get emission files for ODIAC, FFDAS, EDGAR
    # use year 2008 emissions to calculate absolute emission errors
    odiac.file.2008 <- tif2nc.odiacv3(site, timestr = '20081229', vname, workdir,
                                      foot.extent, tiff.path, gzTF = F)

    # get absolute emission eror that can further be convolved with footprints
    # **** in 0.1 deg resolution, need to generate 0.1 deg res of footprint
    emiss.file <- cal.emiss.err(site, timestr, odiac.file.2008, edgar.file,
                                ffdas.file, emiss.file, overwrite = F)

    # txt file name for outputting model results
    txtfile <- file.path(store.path, paste0(timestr, '_', site, '_XCO2ff_emiss_err_', 
                                            abs(nhrs), 'hrs_', oco.sensor, oco.ver, 
                                            '_', met, '_odiac', vname, '.txt'))
  }  # end if run_emiss_err
                                        
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

    sel.emiss <- emiss.df %>% filter(lon >= map.ext[1] & lon <= map.ext[2] &
                                     lat >= map.ext[3] & lat <= map.ext[4])
    print(sel.emiss[sel.emiss$emiss >= 100, ])

    e1 <- mm[[1]] + coord_equal() +
          geom_raster(data = sel.emiss, aes(lon + mm[[3]], lat + mm[[2]],
                      fill = emiss)) +
          scale_fill_gradientn(trans = 'log10', colours = def.col(),
                               limits = c(1, 1E5))
    ggsave(plot = e1, filename = gsub('.nc', '.png', emiss.file),
           width = 8, hright = 8)
  } # end if plotTF


  # if cannot find the correct format of nc file for emissions given selected area
  # and return ODIAC file name with path in front
  if (length(emiss.file) == 0) {
    cat('run.xco2.sim(): NO nc file found, check ODIAC tiff file...\n')
    return()

  } else {
    ## read in emissions
    emiss.dat <- raster(emiss.file)
    emiss.res <- res(emiss.dat)[1]
  }  # end if emiss.file

  # from foot.file, get receptor info
  # fix a minor bug in interpreting footprint filename, DW, 10/11/2018 
  nbin <- str_count(basename(foot.file[1]), '_') + 1
  receptor <- unlist(strsplit(basename(foot.file), '_'))
  receptor <- as.data.frame(matrix(receptor, byrow = T, ncol = nbin),
                            stringsAsFactors = F) %>%
              dplyr::select('V1', 'V2', 'V3') %>% 
              # mutate_all() convert character to numberic
              mutate_all(funs(as.numeric), colnames(receptor)) %>% 
              rename(timestr = V1, lon = V2, lat = V3)

  order.index <- order(receptor$lat)
  receptor  <- receptor[order.index, ]
  foot.file <- foot.file[order.index]
  if (run_emiss_err) {receptor$xco2.ff.err <- NA} else {receptor$xco2.ff <- NA}

  # then loop over each receptor
  for (r in 1:nrow(receptor)) {

    # read in footprint
    foot.dat <- stack(foot.file[r])
    if (nlayers(foot.dat) > 1) foot.dat <- sum(foot.dat)
    #crs(foot.dat) <- '+proj=longlat'

    # NOW, foot and emiss should have the same dimension,
    # multiple them to get contribution map of CO2 enhancements
    xco2.ff.sp <- raster::overlay(x = emiss.dat, y = foot.dat,
                                  fun = function(x, y){return(x * y)}) # spatial xco2.ff
    if (plotTF) plot(log10(xco2.ff.sp))

    # sum the map to get the XCO2 enhancements,
    # note that AK and PW have been incorporated in footprint
    tmp.xco2.ff <- sum(getValues(xco2.ff.sp), na.rm = T)
    
    if (run_emiss_err) {
      receptor$xco2.ff.err[r] <- tmp.xco2.ff
    } else {
      receptor$xco2.ff[r] <- tmp.xco2.ff
    }
    print(tmp.xco2.ff)

    ### store emission * column footprint = XCO2 contribution grid into .nc file
    # store into the same workding dir
    outfile <- file.path(dirname(foot.file[r]),
                         gsub('foot', 'foot_emiss', basename(foot.file[r])))
    print(outfile)

    crs(xco2.ff.sp) <-
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    
    longname <- 'XCO2 enhancemnets due to ODIAC emission'; varname <- 'XCO2'
    if (run_emiss_err) 
      longname <- 'XCO2 error due to ODIAC emission error'; varname <- 'XCO2 error'

    # write raster in nc file
    writeRaster(xco2.ff.sp, outfile, overwrite = TRUE, format = 'CDF',
                varname = varname, varunit = 'PPM', xname = 'lon', yname = 'lat',
                longname = longname)
  }  # end for r

  # finally, write in a txt file
  if (writeTF)
    write.table(x = receptor, file = txtfile, sep = ',', row.names = F, quote = F)

  return(receptor)
} # end of subroutine
