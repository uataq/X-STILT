# subroutine to get ground height using trajwind()
# DW, 10/20/2017

# 'zagl': release height for interpolating ground hgt at receptor lat/lon
# must-have variables for 'varsiwant' -- time, indx, long, lati, zagl, zsfc, temp

# Updates --
# add 0.5 deg GDAS, DW, 01/25/2018
# interpolate ground heights from multiple receptors,
# add vector forms of lat/lon/agl, DW, 05/02/2018
# version 2 for matching Ben's STILT-R version 2, call calc_trajectory DW, 05/29/2018
# allows for multiple agl, when interpolating winds at receptors, DW, 08/29/2018
# allow for prescribing met files, e.g., customized WRF, DW, 08/31/2018

get.ground.hgt <- function(receptor = NULL, agl = 5, run_trajec = T, 
                           varsiwant, conage = 48, cpack = 1, dxf = 1, dyf = 1, 
                           dzf = 0.01, emisshrs = 0.01, frhmax = 3, frhs = 1, 
                           frme = 0.1, frmr = 0, frts = 0.1, frvs = 0.1, 
                           hnf_plume = F, hscale = 10800, ichem = 0, iconvect = 0,
                           initd = 0, isot = 0, kbls = 1, kblt = 1, kdef = 1, 
                           khmax = 9999, kmix0 = 250, kmixd = 3, kmsl = 0, 
                           kpuff = 0, krnd = 6, kspl = 1, kzmix = 1, maxdim = 1, 
                           maxpar = 10000, met_file_format, met_loc,
                           mgmin = 2000, ncycl = 0, ndump = 0, ninit = 1, 
                           n_hours, outdt = 0, outfrac = 0.9, p10f = 1, 
                           qcycle = 0, random = 1, splitf = 1, tkerd = 0.18, 
                           tkern = 0.18, rm_dat = T, rundir, timeout, 
                           tlfrac = 0.1, tratio = 0.9, tvmix = 1, veght = 0.5, 
                           vscale = 200, w_option = 0, zicontroltf = 0, 
                           z_top = 25000){

  # before run trajec, create new 'output'
  output  <- list()

  ### Do not change the following input variables for this simulation!!
  # only release one particle, but turn on turbulance 'nturb'
  delt    <- 1
  nturb   <- T
  numpar  <- 1 * length(agl)
  n_hours <- sign(n_hours)  # only allow for one hours back or forward

  # write output
  timestr <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  output$file <- file.path(rundir, paste0(timestr, '_', receptor$long, '_', 
                                          receptor$lati, '_',
                                          ifelse(length(agl) > 1, 'X', agl), 
                                          '_traj.rds'))

  # replace 'receptor$zagl' with 'agl'
  output$receptor <- list(run_time = receptor$run_time, lati = receptor$lati, 
                          long = receptor$long, zagl = agl)

  # get met files for + or - 1 hour
  met_files <- find_met_files(receptor$run_time, met_file_format, n_hours,
                              met_loc)
  
  particle <- NULL
  if (file.exists(output$file)) particle <- readRDS(output$file)
  
  if (length(particle) == 0 | run_trajec) {
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
      cat('get.ground.hgt(): no particle generated, check hymodelc.out...\n')
      return()
    } else {
      saveRDS(particle, output$file)
      cat(paste('get.ground.hgt():', basename(output$file), 'created\n'))
    } # end if length  
  } 

  # select the min timestep, which is the most closed to the receptor
  sel <- abs(particle$time) == min(unique(abs(particle$time)))
  sel.part <- particle[sel, ]
  nsec <- abs(sel.part$time) * 60  # in second
  
  # grab instantaneous variables
  # in p1 or p2 in distCosine(), first one is longitude, second is latitude
  # distance in x- y- and z- directions [m]
  ubar <- NULL; vbar <- NULL; wbar <- NULL; zsfc <- NULL; temp <- NULL

  # if interpolate met variables for 1+ AGL levels
  for (s in 1:length(agl)) {

    tmp.part <- sel.part[s, ]
    delx <- distCosine(p1 = c(tmp.part$long, receptor$lati),
                       p2 = c(receptor$long, receptor$lati)) # in meter
    dely <- distCosine(p1 = c(receptor$long, tmp.part$lati),
                       p2 = c(receptor$long, receptor$lati))
    delz <- abs(tmp.part$zagl - agl[s])

    # U-, V- and W- velocities [m/s]
    # if n_hours is negative (ie., -1), switch wind direction
    # delx always > 0, no direction, sign(long - long) provides direction on delx
    ubar <- c(ubar, sign(tmp.part$long - receptor$long) * delx / nsec[s] * sign(n_hours))
    vbar <- c(vbar, sign(tmp.part$lati - receptor$lati) * dely / nsec[s] * sign(n_hours))
    wbar <- c(wbar, sign(tmp.part$zagl - agl[s]) * delz / nsec[s] * sign(n_hours))

    # ground height
    zsfc <- c(zsfc, tmp.part$zsfc)
    temp <- c(temp, tmp.part$temp - 273.15)
  } # end for s

  data.frame(ubar, vbar, wbar, zsfc, temp)

}  # end of function
