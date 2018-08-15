# subroutine to get ground height using trajwind()
# DW, 10/20/2017

# 'zagl': release height for interpolating ground hgt at receptor lat/lon

# Updates --
# add 0.5 deg GDAS, DW, 01/25/2018
# interpolate ground heights from multiple receptors,
# add vector forms of lat/lon/agl, DW, 05/02/2018
# version 2 for matching Ben's STILT-R version 2, call calc_trajectory instead,
# DW, 05/29/2018
# add met_files (prescribed or not), DW, 08/08/2018

get.ground.hgt <- function(varsiwant, conage, cpack, dxf, dyf, dzf, emisshrs,
                           frhmax, frhs, frme, frmr, frts, frvs, hscale, ichem,
                           iconvect, initd, isot, ivmax, kbls, kblt, kdef,
                           khmax, kmix0, kmixd, kmsl, kpuff, krnd, kspl, kzmix,
                           maxdim, maxpar, met_files, met_file_format, met_loc,
                           mgmin, ncycl, ndump, ninit, n_hours, outdt, outfrac,
                           p10f, qcycle, random, splitf, tkerd, tkern, rm_dat,
                           receptor, rundir, timeout, tlfrac, tratio, tvmix,
                           veght, vscale, w_option, z_top){

  # before run trajec, create new 'output'
  output  <- list()

  ### Do not change following input variables for this simulation!!
  # only release one particle, but turn on turbulance 'nturb'
  delt    <- 1
  nturb   <- T
  numpar  <- 1
  r_zagl  <- 5   # near-sfc release 5magl
  n_hours <- sign(n_hours)  # only allow for one hours back or forward

  # write output
  timestr <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  output$file <- file.path(rundir, paste0(timestr, '_', receptor$long, '_',
                                        receptor$lati, '_', r_zagl,'.traj.rds'))

  # replace 'receptor$zagl' with 'r_zagl'
  output$receptor <- list(run_time = receptor$run_time, lati = receptor$lati,
                          long = receptor$long, zagl = r_zagl)

  # get met files for + or - 1 hour
  if (is.null(met_files)) {
    met_files <- find_met_files(receptor$run_time, met_file_format,
                                n_hours, met_loc)
  } else {  # if met_files are given
    # select the 1km WRF for interpolating ground height, DW, 08/08/2018
    # 1st always with highest resolution met fields, use the 2nd in this case
    # need further modifications !!!
    met_files <- met_files[1]
  }

  # Execute particle trajectory simulation, and read results into data frame
  if (!file.exists(output$file)) {
    cat('get.ground.hgt(): trajec NOT found...\n')
    # for simply getting ground heights, no need to add ziscale or other wind errors
    particle <- calc_trajectory(varsiwant, conage, cpack, delt, dxf, dyf, dzf,
                                emisshrs, frhmax, frhs, frme, frmr, frts, frvs,
                                hscale, ichem, iconvect, initd, isot, ivmax,
                                kbls, kblt, kdef, khmax, kmix0, kmixd, kmsl,
                                kpuff, krnd, kspl, kzmix, maxdim, maxpar,
                                met_files, mgmin, ncycl, ndump, ninit, numpar,
                                nturb, n_hours, outdt, outfrac, output, p10f,
                                qcycle, random, splitf, tkerd, tkern, rm_dat,
                                timeout, tlfrac, tratio, tvmix, veght, vscale,
                                0, w_option, zicontroltf = F, ziscale = NA,
                                z_top, rundir)

    saveRDS(particle, output$file)   # store traj
    cat(paste(basename(output$file), "created...\n"))

  } else {
    cat('get.ground.hgt(): trajec found...\n')
    particle <- readRDS(output$file)
  }
print(str(particle))
  # select the min timestep, which is the most closed to the receptor
  sel <- abs(particle$time) == min(unique(abs(particle$time)))
  sel.part <- particle[sel, ]
  nsec <- abs(sel.part$time) * 60  # in second

  # grab instantaneous variables
  # in p1 or p2 in distCosine(), first one is longitude, second is latitude
  # distance in x- y- and z- directions [m]
  delx <- distCosine(p1 = c(sel.part$long, receptor$lati),
                     p2 = c(receptor$long, receptor$lati))
  dely <- distCosine(p1 = c(receptor$long, sel.part$lati),
                     p2 = c(receptor$long, receptor$lati))
  delz <- abs(sel.part$zagl - r_zagl)

  # U-, V- and W- velocities [m/s]
  # if n_hours is negative (ie., -1), switch wind direction
  ubar <- sign(sel.part$long - receptor$long) * delx / nsec * sign(n_hours)
  vbar <- sign(sel.part$lati - receptor$lati) * dely / nsec * sign(n_hours)
  wbar <- sign(sel.part$zagl - r_zagl) * delz / nsec * sign(n_hours)

  # ground height
  zsfc <- sel.part$zsfc

  data.frame(ubar, vbar, wbar, zsfc)
}
