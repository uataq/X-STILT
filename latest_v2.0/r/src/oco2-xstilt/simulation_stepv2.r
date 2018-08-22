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

# for debug --
if(F){
  X = 1; r_run_time = receptors$run_time; r_lati = receptors$lati
  r_long = receptors$long; r_zagl = receptors$zagl
}

simulation_stepv2 <- function(X, rm_dat = T, stilt_wd = getwd(), lib.loc = NULL,
                            ak.wgt = NA, conage = 48, cpack = 1, delt = 0,
                            dmassTF = F, dxf = 1, dyf = 1, dzf = 0.01,
                            emisshrs = 0.01, frhmax = 3, frhs = 1, frme = 0.1,
                            frmr = 0, frts = 0.1, frvs = 0.1, hnf_plume = T,
                            hscale = 10800, horcoruverr = NA, horcorzierr = NA,
                            ichem = 0, iconvect = 0, initd = 0, isot = 0,
                            kbls = 1, kblt = 1, kdef = 1, khmax = 9999,
                            kmix0 = 250, kmixd = 3, kmsl = 0, kpuff = 0,
                            krnd = 6, kspl = 1, kzmix = 1, maxdim = 1,
                            maxpar = 10000, met_file_format, met_loc,
                            mgmin = 2000, n_hours = -24, n_met_min = 1,
                            ncycl = 0, ndump = 0, ninit = 1, nturb = 0,
                            numpar = 200, oco2.path = NA, outdt = 0,
                            outfrac = 0.9, output_wd = file.path(stilt_wd,'out'),
                            p10f = 1, projection = '+proj=longlat', pwf.wgt = NA,
                            qcycle = 0, r_run_time, r_lati, r_long, r_zagl,
                            random = 1, run_foot = T, run_trajec = T,
                            siguverr = NA, sigzierr = NA, smooth_factor = 1,
                            splitf = 1, stilt.ver = 2, time_integrate = F,
                            timeout = 3600, tkerd = 0.18, tkern = 0.18,
                            tlfrac = 0.1, tluverr = NA, tlzierr = NA,
                            tratio = 0.9, tvmix = 1, varsiwant = NULL,
                            veght = 0.5, vscale = 200,
                            w_option = 0, xmn = -180, xmx = 180, xres = 0.1,
                            ymn = -90, ymx = 90, yres = xres, zicontroltf = 0,
                            ziscale = NULL, z_top = 25000, zcoruverr = NA) {
  try({
    # If using lapply or parLapply, receptors are passed as vectors and need to
    # be subsetted for the specific simulation index
    # add ziscale, DW
    print(r_run_time)
    if (length(r_run_time) > 1) {
      r_run_time <- r_run_time[X]
      r_lati  <- format(r_lati[X], digits = 4, nsmall = 4)
      r_long  <- format(r_long[X], digits = 4, nsmall = 4)
      r_zagl  <- r_zagl[X]
      ziscale <- ziscale[X]
    }

    # Column trajectories use r_zagl passed as a list of values for lapply and
    # parLapply but a vector in slurm_apply
    r_zagl  <- unlist(r_zagl)
    ziscale <- unlist(ziscale)  # a vector now

    # Ensure dependencies are loaded for current node/process
    source(file.path(stilt_wd, 'r/dependencies.r'), local = T)

    if (is.null(varsiwant)) {
      varsiwant <- c('time', 'indx', 'long', 'lati', 'zagl', 'sigw', 'tlgr',
                     'zsfc', 'icdx', 'temp', 'samt', 'foot', 'shtf', 'tcld',
                     'dmas', 'dens', 'rhfr', 'sphu', 'solw', 'lcld', 'zloc',
                     'dswf', 'wout', 'mlht', 'rain', 'crai')
    } else if (any(grepl('/', varsiwant))) {
      varsiwant <- unlist(strsplit(varsiwant, '/', fixed = T))
    }

    # Creates subdirectories in out for each model run time. Each of these
    # subdirectories is populated with symbolic links to the shared datasets
    # below and a run-specific SETUP.CFG and CONTROL
    rundir_format <- paste0('%Y%m%d%H_', r_long, '_', r_lati, '_',
                            ifelse(length(r_zagl) > 1, 'X', r_zagl))
    rundir  <- file.path(output_wd, 'by-id',
                         strftime(r_run_time, rundir_format, 'UTC'))
    uataq::br()
    message(paste('Running simulation ID:  ', basename(rundir)))

    # Calculate particle trajectories ------------------------------------------
    # run_trajec determines whether to try using existing trajectory files or to
    # recycle existing files
    output <- list()
    output$file <- file.path(rundir, paste0(basename(rundir), '_traj.rds'))

    if (run_trajec) {
      # Ensure necessary files and directory structure are established in the
      # current rundir
      dir.create(rundir)

      # change the path of hymodelc executable, DW, 07/31/2018
      # since we added another version of hymodelc (AER_NOAA_branch) under exe
      link_files <- dir(file.path(stilt_wd, 'exe', 'master'))
      file.symlink(file.path(stilt_wd, 'exe', 'master', link_files), rundir)

      # Execute particle trajectory simulation, and read results into data frame
      output$receptor <- list(run_time = r_run_time,
                              lati = r_lati,
                              long = r_long,
                              zagl = r_zagl)

      ## ---------------- modifications for OCO-2/X-STILT ------------------ ##
      # need modeled ground height to interpolate pres-hgt relation to
      # interpolate satellite weighting profiles from OCO-2 20 levels to model
      # levels, added by Dien Wu, 05/26/2018

      if (length(r_zagl) > 1) {
        cat('Column simulations, estimating modeled ground heights ...\n')

        # store trajec from 5mAGL in the same copy dir 'rundir' by calling
        # get.ground.height() that calls calc_trajectory() to estimate ground
        # height [m] along w. u-, v- and w- component instantaneous wind
        # given receptor lat/lon/time/agl=5 (near ground)
        # remove ziscale and zicontroltf from get.ground.hgt(), DW
        receptor <- output$receptor
        recp.var <- get.ground.hgt(varsiwant, conage, cpack, dxf, dyf, dzf,
                                   emisshrs, frhmax, frhs, frme, frmr, frts,
                                   frvs, hscale, ichem, iconvect, initd, isot,
                                   ivmax, kbls, kblt, kdef, khmax, kmix0, kmixd,
                                   kmsl, kpuff, krnd, kspl, kzmix, maxdim,
                                   maxpar, met_file_format, met_loc, mgmin,
                                   ncycl, ndump, ninit, n_hours, outdt, outfrac,
                                   p10f, qcycle, random, splitf, tkerd, tkern,
                                   rm_dat, receptor, rundir, timeout, tlfrac,
                                   tratio, tvmix, veght, vscale, w_option, z_top)
        # paste interpolated info to output$receptor
        output$receptor <- c(output$receptor, recp.var)
      }  # end if length(r_zagl) > 1
      ## ------------ END modifications for OCO-2/X-STILT ------------------ ##

      # Find necessary met files, if no prescribed met files found
      met_files <- find_met_files(r_run_time, met_file_format, n_hours, met_loc)

      if (length(met_files) < n_met_min) {
        warning('Insufficient amount of meteorological data found...')
        cat('Insufficient amount of meteorological data found. Check ',
            'specifications in run_stilt.r\n',
            file = file.path(rundir, 'ERROR'))
        return()
      } # end if length

      particle <- calc_trajectory(varsiwant, conage, cpack, delt, dxf, dyf, dzf,
                                  emisshrs, frhmax, frhs, frme, frmr, frts, frvs,
                                  hscale, ichem, iconvect, initd, isot, ivmax,
                                  kbls, kblt, kdef, khmax, kmix0, kmixd, kmsl,
                                  kpuff, krnd, kspl, kzmix, maxdim, maxpar,
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
                                          frmr, frts, frvs, hscale, ichem,
                                          iconvect, initd, isot, ivmax, kbls,
                                          kblt, kdef, khmax, kmix0, kmixd, kmsl,
                                          kpuff, krnd, kspl, kzmix, maxdim,
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
                                             horcorzierr = horcorzierr,
                                             ziscale = ziscale)
      }

      # Save output object to compressed rds file and symlink to out/particles
      # directory for convenience
      saveRDS(output, output$file)
      file.symlink(output$file, file.path(output_wd, 'particles',
                                          basename(output$file))) %>%
        invisible()

    } else {
      # If user opted to recycle existing trajectory files, read in the recycled
      # file to a data_frame with an adjusted timestamp and index for the
      # simulation step. If none exists, report an error and proceed
      if (!file.exists(output$file)) {
        warning('simulation_stepv2(): No _traj.rds file found in ', rundir,
                '\n    skipping this timestep and trying the next...')
        return()
      }

      # also need output to get xhgt info, DW, 06/01/2018
      output <- readRDS(output$file)
      particle <- output$particle
    }

    ## ---------------- modifications for OCO-2/X-STILT -------------------- ##
    # Weight footprint: call wgt.trajec.footv3() to weight trajec-level
    # footprint, added by Dien Wu, 06/01/2018
    if (length(r_zagl) > 1) {
      # check whether weighted trajec exists already,
      # directly grab from by-id directory, DW, 07/13/2018
      wgt.file <- file.path(rundir, paste0(basename(rundir), '_wgttraj.rds'))
      cat(paste0('wgttraj file:', wgt.file, '\n'))

      if (file.exists(wgt.file) > 0) {
        wgt.output <- readRDS(wgt.file)

      } else {
        # get OCO-2 profile first according to lat/lon of receptor, return a list
        oco2.info <- get.oco2.info(oco2.path, receptor = output$receptor)

        if (is.null(oco2.info)) {
          warning('NO OCO-2 info found for a given receptor lat/lon'); return()
        } # end if is.null

        # call weight.trajecfootv3() to start weighting trajec-level footprint,
        # before calculating 2D footprint; ak.wgt and pwf.wgt as weighting flags
        output$file <- wgt.file      # overwrite filename
        wgt.output <- wgt.trajec.footv3(output = output, oco2.info = oco2.info,
                                        ak.wgt = ak.wgt, pwf.wgt = pwf.wgt)
      }  # end if file.exists()

      ## returns the ak pw profiles at for all levels & the weighted footprint
      combine.prof <- wgt.output$wgt.prof	 # all vertical profs, ak, pw, apriori

      # overwrite 'particle' with weighted trajec
      particle <- wgt.output$particle
    }  # end if length(r_zagl) > 1
    ## ------------ END modifications for OCO-2/X-STILT -------------------- ##

    # add run_foot for whether to generate footprint
    if (run_foot) {

      ### ------- add Trajecfoot() if stilt.ver = 1, DW, 07/17/2018 -------- ##
      if (stilt.ver == 1) {

        foot_file <- file.path(rundir, paste0(basename(rundir), '_foot.nc'))
        # reform particle to match Trajecfoot
        ensemble <- particle %>%
          dplyr::select(time = time, lat = lati, lon = long, agl = zagl,
                        zi = zsfc, foot = foot, index = indx, dmass = dmas) %>%
          as.matrix()

        foottimes <- c(0, abs(n_hours))
        foot <- Trajecfoot(ident = NULL, part = ensemble, foottimes = foottimes,
                           dmassTF = dmassTF, lon.ll = xmn, lat.ll = ymn,
                           lon.res = xres, lat.res = yres,
                           numpix.x = (xmx - xmn) / xres,
                           numpix.y = (ymx - ymn) / yres) # [lat, lon]

        # change dims to [lon, lat] for write_footprint()
        foot <- aperm(foot, c(2, 1, 3))

        # Set footprint grid
        glong <- head(seq(xmn, xmx, by = xres), -1)
        glati <- head(seq(ymn, ymx, by = yres), -1)

        # Set footprint metadata and write to file
        write_footprint(foot, output = foot_file, glong = glong,
                        glati = glati, projection = '+proj=longlat',
                        xres = xres, yres = yres,
                        time_out = as.numeric(as.POSIXct(r_run_time, tz = "UTC")))
      } else {

        # Calculate near-field dilution height based on gaussian plume width
        # approximation and recalculate footprint sensitivity for cases when the
        # plume height is less than the PBL height scaled by veght
        if (hnf_plume)
          particle <- calc_plume_dilution(particle, numpar, r_zagl, veght)

        # Produce footprint --------------------------------------------------------
        # Aggregate the particle trajectory into surface influence footprints. This
        # outputs a .rds file, which can be read with readRDS() containing the
        # resultant footprint and various attributes
        foot_file <- file.path(rundir, paste0(basename(rundir), '_foot.nc'))
        foot <- calc_footprint(particle, output = foot_file,
                               r_run_time = r_run_time,
                               smooth_factor = smooth_factor,
                               time_integrate = time_integrate,
                               xmn = xmn, xmx = xmx, xres = xres,
                               ymn = ymn, ymx = ymx, yres = yres)
      }  # end of stilt.ver
      ### --- End of adding Trajecfoot() if stilt.ver = 1, DW, 07/17/2018 --- ##

      if (is.null(foot)) {
        warning('No non-zero footprint values found within the footprint domain.')
        cat('No non-zero footprint values found within the footprint domain.\n',
            file = file.path(rundir, 'ERROR'))
        return()
      } else {
        file.symlink(foot_file, file.path(output_wd, 'footprints',
                                          basename(foot_file))) #%>% invisible()
      } # end if is.null(foot)

      invisible(gc())
      return(foot)
    } else {
      cat('No need to generate footprint for this simulation.\n')
      return()
    }# end if run_foot

  })  # try()
}
