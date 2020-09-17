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
# modify the script to the latest STILT-HYSPLIT, DW, 07/03/2020 

get.ground.hgt <- function(receptor = NULL, agl = 0.1, run_trajec = T, 
                           namelist, rundir, emisshrs = 0.01, hnf_plume, 
                           met_file_format, met_path, n_hours, rm_dat = T, 
                           timeout, w_option = 0, z_top = 25000){

  # before running trajec, create new 'tmp.output', 
  # not to use 'output' to avoid conflict
  tmp.output <- list()

  ### Do not change the following input variables for this simulation!!
  # only release one particle, but turn on turbulance 'nturb'
  tmp.namelist <- namelist
  tmp.namelist$delt    <- 1
  tmp.namelist$nturb   <- T
  tmp.namelist$numpar  <- 1 * length(agl)
  n_hours <- sign(n_hours)  # only allow for one hours back or forward

  # write tmp.output
  timestr <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  tmp.output$file <- file.path(rundir, paste0(timestr, '_', receptor$long, '_', 
                                              receptor$lati, '_',
                                              ifelse(length(agl) > 1, 'X', agl), 
                                              '_traj.rds'))

  # replace 'receptor$zagl' with 'agl'
  tmp.output$receptor <- list(run_time = receptor$run_time, lati = receptor$lati, 
                              long = receptor$long, zagl = agl)

  # get met files for + or - 1 hour
  met_files <- find_met_files(receptor$run_time, met_file_format, n_hours, met_path)
  
  particle <- NULL
  if (file.exists(tmp.output$file)) particle <- readRDS(tmp.output$file)

  if (length(particle) == 0 | run_trajec) {
    # modify the code based on Ben's code, DW, 07/03/2020
    particle <- calc_trajectory(namelist = tmp.namelist, rundir, emisshrs, 
                                hnf_plume, met_files, n_hours, output = tmp.output, 
                                rm_dat, timeout, w_option, z_top)

    if (length(particle) == 0) {
      cat('get.ground.hgt(): no particle generated, check by-id...\n')
      return()
    } else {
      saveRDS(particle, tmp.output$file)
      cat(paste('get.ground.hgt():', basename(tmp.output$file), 'created\n'))
    } # end if length  
  } 

  # select the min timestep, which is the most closed to the receptor
  sel <- abs(particle$time) == min(unique(abs(particle$time)))
  sel.part <- particle[sel, ]
  nsec <- abs(sel.part$time) * 60  # in second
  
  # grab instantaneous variables
  # in p1 or p2 in distCosine(), first one is longitude, second is latitude
  # distance in x- y- and z- directions [m]
  ubar <- vbar <- wbar <- zsfc <- psfc <- temp <- NULL

  # if interpolate met variables for 1+ AGL levels
  for (s in 1 : length(agl)) {

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
    psfc <- c(psfc, tmp.part$pres)
    temp <- c(temp, tmp.part$temp - 273.15)
  } # end for s

  data.frame(ubar, vbar, wbar, zsfc, psfc, temp)

}  # end of function
