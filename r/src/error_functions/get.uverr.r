# subroutine to get required input variables for transport error runs in STILT
# Dien Wu, 08/31/2018

# It reads raob data ('grab.raob()') and interpolate met data ('cal.met.wind()')
# if no model-data comparison exists
#
# Will prescribe typical correlation lengthscale and timescale based on met used
# on e can also customize them as an input
#
#' @param run_wind_err T for re-interpolating raob and met winds using get.met.varh()
#' @param err_path path that store wind comparisons of met vs. NOAA RAOB data

get.uverr = function(run_hor_err, site, timestr, xstilt_wd, run_wind_err = F,
                     simstep_namelist = NULL, raob_fn = NULL, nhrs = -72, 
                     met, met_path, met_file_format, lon_lat, maxagl = 3000, 
                     err_path, nfTF = F, siguverr = NULL, TLuverr = NULL, 
                     zcoruverr = NULL, horcoruverr = NULL) {
  
  forwardTF = ifelse(nhrs > 0, TRUE, FALSE)
  
  if (run_hor_err) {
    cat('\n# ------ horizontal wind error component ----- #\n')

    # intput correlation lengthscale (in meter) and timescales (in mins)
    # correlation timescale, horizontal and vertical lengthscales
    if (met == 'gdas0p5') {
      TLuverr = 1 * 60   # in mins
      zcoruverr = 600    # in m
      horcoruverr = 40   # in km

    } else if (met == 'gfs0p25') {
      TLuverr = 1 * 60   # in mins
      zcoruverr = 600    # in m
      horcoruverr = 25   # in km

    } else if (met == 'gdas1') {
      TLuverr = 2.4 * 60
      zcoruverr = 700
      horcoruverr = 97

    } else if (met == 'edas40') {
      TLuverr = 1 * 60
      zcoruverr = 600
      horcoruverr = 40

    } else if (met == 'hrrr') {
      TLuverr = 1 * 60
      zcoruverr = 100
      horcoruverr = 5

    } else {
      cat ('get.uverr(): NO default values found for your met fields...
            Please calculate and insert the time, horizontal and vertical 
            correlation scales (i.e., TLuverr, horcoruverr, zcoruverr), 
            see Wu et al. [2018] and Lin & Gerbig [2005]...\n')
    } # end if met

    # grab modeled winds, *** if no file found, this takes a long time to run
    dir.create(err_path, showWarnings = F, recursive = T)
    err_file = file.path(err_path, paste0(site, '_', met, '_raob_', 
                                          timestr, '.txt'))

    if ( !run_wind_err ) {
      if (file.exists(err_file)) {

          cat('get.uverr(): FOUND existing wind error comparisons...\n')
          met_raob = read.csv(err_file)

          # call get.SIGUVERR() to interpolate most near-field wind errors
          if (!is.null(met_raob)) {
            err_stat = get.siguverr(met_raob, nfTF, forwardTF, lon_lat, 
                                    nhrs, timestr)
            siguverr = err_stat$siguverr
          } else err_stat$siguverr = siguverr   # assign default error

      } else {

        # if one does not want to run wind error analysis and no RAOB file found
        cat('get.uverr(): no wind error comparisons found; 
             make a conservative assumption for siguverr;
             use user-defined RMSE wind error OR a default value of 3 m/s...\n')
        err_stat = NULL

        # make a conservative assumption about the wind error, or use prescribed #
        # < 2 m/s for GDAS 0.5deg for the Middle East, based on Wu et al., 2018
        if (is.null(siguverr)) siguverr = 3
        err_stat$siguverr = siguverr

      } # end if not running wind error

    } else {
      
      # run the entire wind interpolation and get error statistics
      cat('get.uverr(): run_wind_err = TRUE --> Need to calculate wind error according to RAOB data...it takes a while\n')
      met_raob = cal.wind.err(err_file, met, met_path, met_file_format, 
                              xstilt_wd, err_path, site, timestr, 
                              simstep_namelist, raob_fn, nhrs, maxagl)
    
      # call get.SIGUVERR() to interpolate most near-field wind errors
      if (!is.null(met_raob)) 
      err_stat = get.siguverr(met_raob, nfTF, forwardTF, lon_lat, nhrs, timestr)
      
      if (is.null(met_raob) | is.null(err_stat)) { 
        cat('get.uverr(): no raob data found, use default value of 3 m/s \n')
        err_stat$siguverr = siguverr
      } else siguverr = err_stat$siguverr

    }  # end if 

    cat(paste('SIGUVERR:', signif(siguverr, 3), 'm/s..\n'))
    err_stat$TLuverr = TLuverr
    err_stat$zcoruverr = zcoruverr
    err_stat$horcoruverr = horcoruverr
    err_stat$run_hor_err = run_hor_err

  } else {
    cat('NO horizontal wind error component for generating trajec...\n')
    err_stat = data.frame(run_hor_err = run_hor_err, siguverr = NA, 
                          u.bias = NA, v.bias = NA, ws.bias = NA, wd.bias = NA,
                          TLuverr = NA, zcoruverr = NA, horcoruverr = NA)

  } # end if run_hor_err

  return(err_stat)
}
