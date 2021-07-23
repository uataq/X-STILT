# subroutine to get met temp/wind/hgt for a certain hgt, DW, 10/20/2017

# 'zagl': release height for interpolating ground hgt at receptor lat/lon
# must-have variables for 'varsiwant' -- time, indx, long, lati, zagl, zsfc, temz
# !!! temz for temp at certain height; temp for temp at the lowest height 
#' @param namelist here is the one used in simulation_step in STILT

# Updates --
# add 0.5 deg GDAS, DW, 01/25/2018
# interpolate ground heights from multiple receptors,
# add vector forms of lat/lon/agl, DW, 05/02/2018
# version 2 for matching Ben's STILT-R version 2, call calc_trajectory DW, 05/29/2018
# allows for multiple agl, when interpolating winds at receptors, DW, 08/29/2018
# allow for prescribing met files, e.g., customized WRF, DW, 08/31/2018
# modify the script to the latest STILT-HYSPLIT, DW, 07/03/2020 
# update the script since no discrete level allowed for column release, DW, 02/27/2021
# **** no allowance for multiple AGL now...

get.met.varh <- function(receptor = NULL, agl = 5, run_trajec = T, namelist, 
                         rundir, hnf_plume = F, met_file_format, met_path, 
                         n_hours, timeout, z_top = 25000){
  
  if (length(agl) > 1) stop('get.met.varh(): only one AGL allowed...\n')

  # before running trajec, create new 'tmp.output', 
  # not to use 'output' to avoid conflict
  tmp.output <- list()
  n_hours <- sign(n_hours)  # only allow for one hours back or forward

  # write tmp.output
  timestr <- strftime(receptor$run_time, tz = 'UTC', format = '%Y%m%d%H')
  tmp.output$file <- file.path(rundir, paste0(timestr, '_', receptor$long, '_', 
                                              receptor$lati, '_', agl, '_traj.rds'))

  # replace 'receptor$zagl' with 'agl'
  tmp.output$receptor <- list(run_time = receptor$run_time, 
                              lati = receptor$lati, 
                              long = receptor$long, 
                              zagl = agl)

  # get met files for + or - 1 hour
  met_files <- find_met_files(t_start = receptor$run_time, met_file_format, 
                              n_hours, met_path)  
  
  if (length(met_files) == 0) 
    stop('get.met.varh(): NO meteo fields found; please check\n')

  particle <- NULL
  if (file.exists(tmp.output$file)) particle <- readRDS(tmp.output$file)
  
  # modify the code based on Ben's code, DW, 02/26/2021
  if (length(particle) == 0 | run_trajec) {
    namelist$numpar = 1
    particle <- calc_trajectory(namelist, rundir, emisshrs = 0.01, hnf_plume, 
                                met_files, n_hours, output = tmp.output, 
                                rm_dat = T, timeout, w_option = 0, z_top)

    if (length(particle) == 0) {
      cat('get.met.varh(): no particle generated, check by-id...\n')
      return()
    } else {
      saveRDS(particle, tmp.output$file)
      #cat(paste('get.met.varh():', basename(tmp.output$file), 'created\n'))
    } # end if length  
  }   # end if

  # select the min timestep, which is the most closed to the receptor
  sel <- abs(particle$time) == min(unique(abs(particle$time)))
  sel.part <- particle[sel, ]
  nsec <- abs(sel.part$time) * 60  # in second
  
  # grab instantaneous variables
  # in p1 or p2 in distCosine(), first one is longitude, second is latitude
  # distance in x- y- and z- directions [m]
  delx <- distCosine(p1 = c(sel.part$long, receptor$lati),
                     p2 = c(receptor$long, receptor$lati)) # in meter
  dely <- distCosine(p1 = c(receptor$long, sel.part$lati),
                     p2 = c(receptor$long, receptor$lati))
  delz <- abs(sel.part$zagl - agl)

  # U-, V- and W- velocities [m/s]
  # if n_hours is negative (ie., -1), switch wind direction
  # delx always > 0, no direction, sign(long - long) provides direction on delx
  ubar <- sign(sel.part$long - receptor$long) * delx / nsec * sign(n_hours)
  vbar <- sign(sel.part$lati - receptor$lati) * dely / nsec * sign(n_hours)
  wbar <- sign(sel.part$zagl - agl) * delz / nsec * sign(n_hours)

  # ground height
  zsfc <- sel.part$zsfc
  pres <- sel.part$pres
  temz <- sel.part$temz - 273.15
  data.frame(ubar, vbar, wbar, zsfc, pres, temz)

}  # end of function
