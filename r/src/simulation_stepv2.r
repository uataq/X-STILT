#' simulation_step runs STILT for the given timestep
#' @author Ben Fasoli
#'
#' Executes trajectory calculations with (and optionally without) transport
#' error and calculates kernel density derived footprint grids.
#'
#' For documentation, see https://uataq.github.io/stilt/
#'
#' @export

# modification for OCO-2/X-STILT, Dien Wu, 06/01/2018
# add input variables including 'oco2.path', 'ak.wgt', 'pwf.wgt'
# extra functions including 'get.oco2.info()', 'get.ground.hgt()',
#                           'get.wgt.funcv3()', 'wgt.trajec.footv3()'
# add 'run_foot', if false, then no need to run footprint. DW, 06/04/2018
# add STILTv1 and Trajecfoot(), DW, 07/17/2018
# add ziscale as list, remember to unlist, DW, 07/25/2018
# change the path of hymodelc executable, DW, 07/31/2018
# add PW and Ak weighting for trajec with error as well, DW, 10/21/2018
# add horizontal transport error here, DW, 10/21/2018 
# separate X-STILT modification fron STILT, DW, 01/20/2019

if (F) {
  X = 1
  r_run_time = receptors$run_time[X]
  r_lati = receptors$lati[X]
  r_long = receptors$long[X]
  r_zagl = receptors$zagl[X]
}
simulation_stepv2 <- function(before_traj = NULL, # can be a function for X-STILT
                              before_foot = NULL, # or can be a function for 
                                                  # weighting foot column for trajec
                              ak.wgt = NA, 
                              oco2.path = NULL, 
                              pwf.wgt = NA,

                              conage = 48,
                              cpack = 1,
                              delt = 0,
                              dxf = 1,
                              dyf = 1,
                              dzf = 0.01,
                              emisshrs = 0.01,
                              frhmax = 3,
                              frhs = 1,
                              frme = 0.1,
                              frmr = 0,
                              frts = 0.1,
                              frvs = 0.1, 
                              hnf_plume = T,
                              hscale = 10800,
                              horcoruverr = NA, 
                              horcorzierr = NA,
                              ichem = 0,
                              iconvect = 0, 
                              initd = 0,
                              isot = 0,
                              kbls = 1,
                              kblt = 1,
                              kdef = 1,
                              khmax = 9999,
                              kmix0 = 250,
                              kmixd = 3,
                              kmsl = 0, 
                              kpuff = 0,
                              krnd = 6,
                              kspl = 1,
                              kzmix = 1,
                              lib.loc = NULL,
                              maxdim = 1,
                              maxpar = max(10000, numpar),
                              met_file_format, 
                              met_loc,
                              mgmin = 2000,
                              n_hours = -24,
                              n_met_min = 1,
                              ncycl = 0,
                              ndump = 0,
                              ninit = 1,
                              nturb = 0,
                              numpar = 200,
                              outdt = 0,
                              outfrac = 0.9,
                              output_wd = file.path(stilt_wd, 'out'),
                              p10f = 1,
                              projection = '+proj=longlat',
                              qcycle = 0, 
                              r_run_time,
                              r_lati,
                              r_long,
                              r_zagl,
                              random = 1, 
                              rm_dat = T,
                              run_foot = T,
                              run_trajec = T,
                              siguverr = NA,
                              sigzierr = NA, 
                              smooth_factor = 1,
                              splitf = 1,
                              stilt_wd = getwd(),
                              time_integrate = F,
                              timeout = 3600,
                              tkerd = 0.18,
                              tkern = 0.18,
                              tlfrac = 0.1,
                              tluverr = NA,
                              tlzierr = NA, 
                              tratio = 0.9,
                              tvmix = 1,
                              varsiwant = c('time', 'indx', 'long', 'lati', 'zagl', 
                                            'sigw', 'tlgr', 'zsfc', 'icdx', 'temp',
                                            'samt', 'foot', 'shtf', 'tcld', 'dmas', 
                                            'dens', 'rhfr', 'sphu', 'solw', 'lcld', 
                                            'zloc', 'dswf', 'wout', 'mlht', 'rain', 
                                            'crai', 'pres'), 
                              veght = 0.5,
                              vscale = 200,
                              w_option = 0,
                              xmn,
                              xmx,
                              xres,
                              ymn,
                              ymx, 
                              yres = xres,
                              zicontroltf = 0,
                              ziscale = 0,
                              z_top = 25000,
                              zcoruverr = NA) {
  try({
    setwd(stilt_wd)
    
    # Vector style arguments passed as a list
    r_zagl <- unlist(r_zagl)
    varsiwant <- unlist(varsiwant)
    ziscale <- unlist(ziscale)

    # Validate arguments
    if (!run_trajec && !run_foot)
      stop('simulation_step(): Nothing to do, set run_trajec or run_foot to T')
    
    # Ensure dependencies are loaded for current node/process
    source(file.path(stilt_wd, 'r/dependencies.r'), local = T)
     
    # Creates subdirectories in out for each model run time. Each of these
    # subdirectories is populated with symbolic links to the shared datasets
    # below and a run-specific SETUP.CFG and CONTROL
    rundir_format <- paste0('%Y%m%d%H_', r_long, '_', r_lati, '_', 
                            ifelse(length(r_zagl) > 1, 'X', r_zagl))
    rundir  <- file.path(output_wd, 'by-id',
                         strftime(r_run_time, rundir_format, 'UTC'))
    dir.create(rundir, showWarnings = F, recursive = T)
    dir.create(file.path(output_wd, 'particles'), showWarnings = F, recursive = T)
    dir.create(file.path(output_wd, 'footprints'), showWarnings = F, recursive = T)
    message(paste('Running simulation ID:  ', basename(rundir)))
    
    # Calculate particle trajectories ------------------------------------------
    # run_trajec determines whether to try using existing trajectory files or to
    # recycle existing files
    output <- list()
    output$file <- file.path(rundir, paste0(basename(rundir), '_traj.rds'))

    if (run_trajec) {
      # Ensure necessary files and directory structure are established in the
      # current rundir
      if (!dir.exists(rundir)) dir.create(rundir)
      #link_files <- dir(file.path(stilt_wd, 'exe'))
      # change the path of hymodelc executable, DW, 07/31/2018
      # since we added another version of hymodelc (AER_NOAA_branch) under exe
      link_files <- dir(file.path(stilt_wd, 'exe', 'master'))
      suppressWarnings(
        #file.symlink(file.path(stilt_wd, 'exe', link_files), rundir)
        file.symlink(file.path(stilt_wd, 'exe', 'master', link_files), rundir)
      )
      
      # Find necessary met files
      met_files <- find_met_files(r_run_time, met_file_format, n_hours, met_loc)
      if (length(met_files) < n_met_min) {
        msg <- paste('Insufficient number of meteorological files found. Check',
                     'specifications in run_stilt.r')
        warning(msg)
        cat(msg, '\n', file = file.path(rundir, 'ERROR'))
        return()
      }
      
      # Execute particle trajectory simulation, and read results into data frame
      output$receptor <- list(run_time = r_run_time,
                              lati = r_lati,
                              long = r_long,
                              zagl = r_zagl)
      
      # check whether before_trajec() exist, if yes, call before_trajec
      # X-STILT needs to interpolate ground height before calc_trajectory
      # DW, 01/18/2019
      if (!is.null(before_trajec)) 
        output <- before_trajec(output, r_zagl = 5, run_trajec, varsiwant, conage, 
                                cpack, dxf, dyf, dzf, emisshrs, frhmax, frhs, 
                                frme, frmr, frts, frvs, hnf_plume = F, hscale, 
                                ichem, iconvect, initd, isot, kbls, kblt, kdef, 
                                khmax, kmix0, kmixd, kmsl, kpuff, krnd, kspl, 
                                kzmix, maxdim, maxpar, met_file_format, met_loc,
                                mgmin, ncycl, ndump, ninit, n_hours, outdt, 
                                outfrac, p10f, qcycle, random, splitf, tkerd, 
                                tkern, rm_dat, rundir, timeout, tlfrac, tratio,
                                tvmix, veght, vscale, w_option, zicontroltf, 
                                z_top)

      particle <- calc_trajectory(varsiwant, conage, cpack, delt, dxf, dyf, dzf,
                                  emisshrs, frhmax, frhs, frme, frmr, frts, frvs,
                                  hnf_plume, hscale, ichem, iconvect, initd, isot,
                                  ivmax, kbls, kblt, kdef, khmax, kmix0, kmixd,
                                  kmsl, kpuff, krnd, kspl, kzmix, maxdim, maxpar,
                                  met_files, mgmin, ncycl, ndump, ninit, numpar,
                                  nturb, n_hours, outdt, outfrac, output, p10f,
                                  qcycle, random, splitf, tkerd, tkern, rm_dat,
                                  timeout, tlfrac, tratio, tvmix, veght, vscale,
                                  0, w_option, zicontroltf, ziscale, z_top, 
                                  rundir)
      if (is.null(particle)) return()

      # Bundle trajectory configuration metadata with trajectory informtation
      output$particle <- particle
      output$params <- read_config(file = file.path(rundir, 'CONC.CFG'))

      # Optionally execute second trajectory simulations to quantify transport
      # error using parameterized correlation length and time scales
      xyerr <- write_winderr(siguverr, tluverr, zcoruverr, horcoruverr,
                             file.path(rundir, 'WINDERR'))
      zerr <- write_zierr(sigzierr, tlzierr, horcorzierr,
                          file = file.path(rundir, 'ZIERR'))
      winderrtf <- (!is.null(xyerr)) + 2 * !is.null(zerr)
      if (winderrtf > 0) {
        particle_error <- calc_trajectory(varsiwant, conage, cpack, delt, dxf,
                                          dyf, dzf, emisshrs, frhmax, frhs, frme,
                                          frmr, frts, frvs, hnf_plume, hscale,
                                          ichem, iconvect, initd, isot, ivmax,
                                          kbls, kblt, kdef, khmax, kmix0, kmixd,
                                          kmsl, kpuff, krnd, kspl, kzmix, maxdim,
                                          maxpar, met_files, mgmin, ncycl, ndump,
                                          ninit, numpar, nturb, n_hours, outdt,
                                          outfrac, output, p10f, qcycle, random,
                                          splitf, tkerd, tkern, rm_dat, timeout,
                                          tlfrac, tratio, tvmix, veght, vscale,
                                          winderrtf, w_option, zicontroltf,
                                          ziscale, z_top, rundir)
        if (is.null(particle_error)) return()
        output$particle_error <- particle_error
        output$particle_error_params <- list(siguverr = siguverr,
                                             tluverr = tluverr,
                                             zcoruverr = zcoruverr,
                                             horcoruverr = horcoruverr,
                                             sigzierr = sigzierr,
                                             tlzierr = tlzierr,
                                             horcorzierr = horcorzierr)
      }
    
      # Save output object to compressed rds file and symlink to out/particles
      saveRDS(output, output$file)
      invisible(file.symlink(output$file, file.path(output_wd, 'particles',
                                                    basename(output$file))))
      # Exit if not performing footprint calculations
      if (!run_foot) return(invisible(output$file))

    } else {
      # If user opted to recycle existing trajectory files, read in the recycled
      # file to a data_frame with an adjusted timestamp and index for the
      # simulation step. If none exists, report an error and proceed
      if (!file.exists(output$file)) {
        warning('simulation_step(): No _traj.rds file found in ', rundir,
                '\n    skipping this timestep and trying the next...')
        return()
      }

      # also need output to get xhgt info, DW, 06/01/2018
      output <- readRDS(output$file)
      particle <- output$particle
    }

    # if before_footprint exists, call it to weight footprint column for X-STILT 
    if (!is.null(before_footprint))
      particle <- before_footprint(output, rundir, oco2.path, ak.wgt, pwf.wgt)

    # Produce footprint --------------------------------------------------------
    # Aggregate the particle trajectory into surface influence footprints. This
    # outputs a .rds file, which can be read with readRDS() containing the
    # resultant footprint and various attributes
    
    # add foot.res on foot_file, DW, 09/15/2018
    foot_file <- file.path(rundir, paste0(basename(rundir), '_', 
                                          signif(xres, 3), 'x', signif(yres, 3), 
                                          '_foot.nc'))
    foot <- calc_footprint(particle, output = foot_file,
                           r_run_time = r_run_time,
                           smooth_factor = smooth_factor,
                           time_integrate = time_integrate,
                           xmn = xmn, xmx = xmx, xres = xres,
                           ymn = ymn, ymx = ymx, yres = yres)
    if (is.null(foot)) {
      msg <- 'No non-zero footprint values found within the footprint domain.'
      warning(msg)
      cat(msg, '\n', file = file.path(rundir, 'ERROR'))
      return()
    }
    
    # Symlink footprint to out/footprints
    invisible(file.symlink(foot_file, file.path(output_wd, 'footprints',
                                                basename(foot_file))))

    invisible(gc())
    return(foot)
  })
}
