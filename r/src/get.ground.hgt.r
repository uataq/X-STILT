# subroutine to get ground height using trajwind()
# DW, 10/20/2017

# 'zagl': release height for interpolating ground hgt at receptor lat/lon
# must-have variables for 'varsiwant' -- time, indx, long, lati, zagl, zsfc, temp

# Updates --
# add 0.5 deg GDAS, DW, 01/25/2018
# interpolate ground heights from multiple receptors,
# add vector forms of lat/lon/agl, DW, 05/02/2018
# version 2 for matching Ben's STILT-R version 2, call calc_trajectory DW, 05/29/2018
# allows for multiple r_zagl, when interpolating winds at receptors, DW, 08/29/2018
# allow for prescribing met files, e.g., customized WRF, DW, 08/31/2018
# simplify the code, DW, 01/18/2019

get.ground.hgt <- function(receptor = NULL, r_zagl = 5, run_trajec = T, 
                           varsiwant, conage, cpack, dxf, dyf, dzf, emisshrs, 
                           frhmax, frhs, frme, frmr, frts, frvs, hnf_plume = F, 
                           hscale, ichem, iconvect, initd, isot, kbls, kblt, 
                           kdef, khmax, kmix0, kmixd, kmsl, kpuff, krnd, kspl, 
                           kzmix, maxdim, maxpar, met_file_format, met_loc,
                           mgmin, ncycl, ndump, ninit, n_hours, outdt, outfrac, 
                           p10f, qcycle, random, splitf, tkerd, tkern, 
                           rm_dat = T, rundir, timeout, tlfrac, tratio, 
                           tvmix, veght, vscale, w_option, zicontroltf, z_top){

  # before run trajec, create new 'output'
  output  <- list()

  ### Do not change the following input variables for this simulation!!
  # only release one particle, but turn on turbulance 'nturb'
  delt    <- 1
  nturb   <- T
  numpar  <- 1 * length(r_zagl)
  n_hours <- sign(n_hours)  # only allow for one hours back or forward

  # write output
  timestr <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  output$file <- file.path(rundir, paste0(timestr, '_', receptor$long, '_', 
                                          receptor$lati, '_',
                                          ifelse(length(r_zagl) > 1, 'X', r_zagl), 
                                          '_traj.rds'))

  # replace 'receptor$zagl' with 'r_zagl'
  output$receptor <- list(run_time = receptor$run_time, lati = receptor$lati, 
                          long = receptor$long, zagl = r_zagl)

  # get met files for + or - 1 hour
  met_files <- find_met_files(receptor$run_time, met_file_format, n_hours,
                              met_loc)

  if (run_trajec) {
    cat('get.ground.hgt(): NO trajec found...\n')
    # for simply getting ground heights, no need to add ziscale or other wind errors
    # modify the code based on Ben's code, 01/20/2019 
    particle <- calc_trajectory(varsiwant, conage, cpack, delt, dxf, dyf, dzf,
                                emisshrs, frhmax, frhs, frme, frmr, frts, frvs,
                                hnf_plume, hscale, ichem, iconvect, initd, isot,
                                ivmax, kbls, kblt, kdef, khmax, kmix0, kmixd,
                                kmsl, kpuff, krnd, kspl, kzmix, maxdim, maxpar,
                                met_files, mgmin, ncycl, ndump, ninit, numpar,
                                nturb, n_hours, outdt, outfrac, output, p10f,
                                qcycle, random, splitf, tkerd, tkern, rm_dat,
                                timeout, tlfrac, tratio, tvmix, veght, vscale,
                                winderrtf = 0, w_option, zicontroltf, ziscale, 
                                z_top, rundir)

    if (length(particle) == 0) {
      cat('get.ground.hgt(): no particle can be generated, check hymodelc.out...\n')
      return()
    } else {
      saveRDS(particle, output$file)
      cat(paste(basename(output$file), 'created\n'))
    } # end if length

  } else {
    particle <- readRDS(output$file)
    cat('get.ground.hgt(): trajec found...\n')
  } # end if run_trajec

  # select the min timestep, which is the most closed to the receptor
  sel <- abs(particle$time) == min(unique(abs(particle$time)))
  sel.part <- particle[sel, ]
  nsec <- abs(sel.part$time) * 60  # in second
  
  # grab instantaneous variables
  # in p1 or p2 in distCosine(), first one is longitude, second is latitude
  # distance in x- y- and z- directions [m]
  ubar <- NULL; vbar <- NULL; wbar <- NULL; zsfc <- NULL; temp <- NULL

  # if interpolate met variables for 1+ AGL levels
  for (s in 1:length(r_zagl)) {

    tmp.part <- sel.part[s, ]
    delx <- distCosine(p1 = c(tmp.part$long, receptor$lati),
                       p2 = c(receptor$long, receptor$lati)) # in meter
    dely <- distCosine(p1 = c(receptor$long, tmp.part$lati),
                       p2 = c(receptor$long, receptor$lati))
    delz <- abs(tmp.part$zagl - r_zagl[s])

    # U-, V- and W- velocities [m/s]
    # if n_hours is negative (ie., -1), switch wind direction
    # delx always > 0, no direction, sign(long - long) provides direction on delx
    ubar <- c(ubar, sign(tmp.part$long - receptor$long) * delx / nsec[s] * sign(n_hours))
    vbar <- c(vbar, sign(tmp.part$lati - receptor$lati) * dely / nsec[s] * sign(n_hours))
    wbar <- c(wbar, sign(tmp.part$zagl - r_zagl[s]) * delz / nsec[s] * sign(n_hours))

    # ground height
    zsfc <- c(zsfc, tmp.part$zsfc)
    temp <- c(temp, tmp.part$temp - 273.15)
  } # end for s

  data.frame(ubar, vbar, wbar, zsfc, temp)

}  # end of function
