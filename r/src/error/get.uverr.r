# subroutine to get required input variables for transport error runs in STILT
# Dien Wu, 08/31/2018

# It reads raob data ('grab.raob()') and interpolate met data ('cal.met.wind()')
# if no model-data comparison exists
#
# Will prescribe typical correlation lengthscale and timescale based on met used
# on e can also customize them as an input
#
#' @param run_wind_err T for re-interpolating raob and met winds using get.met.varh()
#' @param err.path path that store wind comparisons of met vs. NOAA RAOB data

get.uverr <- function(run_hor_err, site, timestr, xstilt_wd, run_wind_err = F,
                      simstep_namelist = NULL, raob.path = NULL, 
                      raob.format = c('fsl', 'ncdf'), nhrs = -72, met, met_path, 
                      met_file_format, lon.lat, agl = c(0, 3000), err.path, 
                      nfTF = F, siguverr = NULL, TLuverr = NULL, 
                      zcoruverr = NULL, horcoruverr = NULL) {
  
  forwardTF <- ifelse(nhrs > 0, TRUE, FALSE)
  
  if (run_hor_err) {
    cat('\n# ------ horizontal wind error component ----- #\n')

    # intput correlation lengthscale (in meter) and timescales (in mins)
    # correlation timescale, horizontal and vertical lengthscales
    if (met == 'gdas0p5') {
      TLuverr <- 1 * 60   # in mins
      zcoruverr <- 600    # in m
      horcoruverr <- 40   # in km

    } else if (met == 'gfs0p25') {
      TLuverr <- 1 * 60   # in mins
      zcoruverr <- 600    # in m
      horcoruverr <- 25   # in km

    } else if (met == 'gdas1') {
      TLuverr <- 2.4 * 60
      zcoruverr <- 700
      horcoruverr <- 97

    } else if (met == 'edas40') {
      TLuverr <- 1 * 60
      zcoruverr <- 600
      horcoruverr <- 40

    } else if (met == 'hrrr') {
      TLuverr <- 1 * 60
      zcoruverr <- 100
      horcoruverr <- 5

    } else {
      cat ('get.uverr(): NO default values found for your met fields...
            Please calculate and insert the time, horizontal and vertical 
            correlation scales (i.e., TLuverr, horcoruverr, zcoruverr), 
            see Wu et al. [2018] and Lin & Gerbig [2005]...\n')
    } # end if met

    # grab modeled winds, *** if no file found, this takes a long time to run
    dir.create(err.path, showWarnings = F, recursive = T)
    err.file <- file.path(err.path,  paste0(site, '_', met, '_rad_', timestr, '.txt'))

    if ( !file.exists(err.file) & run_wind_err == F ) {

      # if one does not want to run wind error analysis and no RAOB file found
      cat('get.uverr(): no wind error comparisons found; 
           make a conservative assumption for siguverr;
           use user-defined RMSE wind error OR a default value of 3 m/s...\n')
      err.stat <- NULL

      # make a conservative assumption about the wind error, or use prescribed #
      # < 2 m/s for GDAS 0.5deg for the Middle East, based on Wu et al., 2018, GMD
      if (is.null(siguverr)) siguverr <- 3
      err.stat$siguverr <- siguverr

    } else {
      
      # run the entire wind interpolation and get error statistics
      cat('get.uverr(): Need to calculate wind error according to RAOB data...it takes a while\n')
      met.raob <- cal.wind.err(err.file, met, met_path, met_file_format, xstilt_wd, 
                               err.path, site, timestr, simstep_namelist, raob.path, 
                               raob.nhrs = -72, raob.format, err.hgt = max(agl))

      # call get.SIGUVERR() to interpolate most near-field wind errors
      err.stat <- get.siguverr(met.raob, nfTF, forwardTF, lon.lat, nhrs, agl, timestr)
      
      if (is.null(err.stat)) { 
        cat('get.uverr(): no raob data found, use default value of 3 m/s \n')
        err.stat$siguverr <- siguverr
      } else siguverr <- err.stat$siguverr

    }  # end if file.exists

    cat(paste('SIGUVERR:', signif(siguverr, 3), 'm/s..\n'))

    err.stat$TLuverr   <- TLuverr
    err.stat$zcoruverr <- zcoruverr
    err.stat$horcoruverr <- horcoruverr
    err.stat$run_hor_err <- run_hor_err

  } else {
    cat('NO horizontal wind error component for generating trajec...\n')
    err.stat <- data.frame(run_hor_err = run_hor_err, siguverr = NA, 
                           u.bias = NA, v.bias = NA, ws.bias = NA, wd.bias = NA,
                           TLuverr = NA, zcoruverr = NA, horcoruverr = NA)

  } # end if run_hor_err

  return(err.stat)
}
