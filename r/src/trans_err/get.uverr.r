# subroutine to get required input variables for transport error runs in STILT
# Dien Wu, 08/31/2018

# It reads raob data ('grab.raob()') and interpolate met data ('cal.met.wind()')
# if no model-data comparison exists
#
# Will prescribe typical correlation lengthscale and timescale based on met used
# on e can also customize them as an input
#
# 'overwrite': T for re-interpolating raob and met winds (using get.ground.hgt.r)
# 'err.path': path that store wind comparisons of met vs. NOAA RAOB data

get.uverr <- function(run_hor_err, site, timestr, workdir, overwrite = F,
                      raob.path = NULL, raob.format = c('fsl', 'ncdf'), 
                      nhrs = -120, met = c('gdas1', 'gdas0p5', 'hrrr')[2], 
                      met.path, met.format, lon.lat, agl, nfTF = F, err.path, 
                      siguverr = NULL, TLuverr = NULL, zcoruverr = NULL, 
                      horcoruverr = NULL, met.files = NULL) {
  
  forwardTF <- ifelse(nhrs > 0, TRUE, FALSE)
  
  if (run_hor_err) {
    cat('+++ horizontal wind error component +++\n')

    # intput correlation lengthscale (in meter) and timescales (in mins)
    # correlation timescale, horizontal and vertical lengthscales
    if (met == 'gdas0p5') {
      TLuverr <- 1 * 60
      zcoruverr <- 600
      horcoruverr <- 40

    } else if (met == 'gdas1') {
      TLuverr <- 2.4 * 60
      zcoruverr <- 700
      horcoruverr <- 97

    } else if (met == 'edas40') {
      TLuverr <- 1 * 60
      zcoruverr <- 600
      horcoruverr <- 40

    } else {
      cat ('get.uverr(): NO default values found...
            Please input the time, horizontal and vertical 
            wind error covariance (TLuverr, horcoruverr, zcoruverr) 
            for your met fields...\n')
    } # end if met
    
    # grab modeled winds, *** if no file found, this takes a long time to run
    raob.file <- file.path(err.path, 
                           paste0(site, '_', met, '_rad_', timestr, '.txt'))

    if (!file.exists(raob.file) & overwrite == F) {

      # if one does not want to run wind error analysis and no RAOB file found
      cat('get.uverr(): no wind error comparisons found; 
           make a consevative assumption for siguverr;
           use user-defined RMSE wind error,
           OR a default value of 2 m/s from Wu et al. [2018]...\n')
      err.stat <- NULL

      # make a conservative assumption about the wind error, or use prescribed #
      # < 2 m/s for GDAS 0.5deg for the Middle East, based on Wu et al., 2018, GMD
      if (is.null(siguverr)) siguverr <- 2
      err.stat$siguverr <- siguverr

    } else {
      # run the entire wind interpolation and get error statistics
      met.raob <- cal.met.wind(filename = raob.file, met, met.path, met.format, 
                               met.files, workdir, site, timestr, overwrite, 
                               raob.path, nhrs = -120, raob.format)

      # call get.SIGUVERR() to interpolate most near-field wind errors
      err.stat <- get.siguverr(met.raob, nfTF, forwardTF, lon.lat, nhrs, agl, 
                               recp.time = timestr)
      if (!is.null(err.stat)) siguverr <- err.stat$siguverr
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
