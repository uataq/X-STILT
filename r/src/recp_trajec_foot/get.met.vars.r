# function to obtain specific humidity and temp profiles from adopted met fields 
# DW, 08/26/2020 

# get specific humidity and temp profiles that will be used to calculate 
#   dry-air column density in mol m-2 for further calculating PWF and 
#   converting XCO from mol m-2 to ppb


if (F) {
    met_file_format = '%Y%m%d'
    met_path = file.path(homedir, 'gfs0p25')
    z_top = 25000
    output$receptor <- list(run_time = r_run_time, lati = r_lati,
                            long = r_long, zagl = r_zagl)
}

# still need specific humidity profiles to calculate c.dry per layer 
# returns 'output'
get.met.vars <- function(namelist, output, met_file_format, met_path, z_top = 25000) {

    tmp.namelist <- namelist
    tmp.namelist$delt <- 1
    tmp.namelist$nturb <- T
    tmp.namelist$numpar <- z_top / 10
    tmp.namelist$varsiwant <- c('time', 'indx', 'long', 'lati', 'zagl', 'zsfc', 
                                'foot', 'mlht', 'dens', 'samt', 'sigw', 'tlgr', 
                                'temp', 'pres', 'sphu', 'rhfr', 'temz')

    # write tmp.output
    tmp.nhrs <- -1
    tmp.output <- list()
    tmp.output$file <- gsub('X_traj.rds', 'met_vars.rds', output$file)
    tmp.output$receptor <- data.frame(run_time = output$receptor$run_time, 
                                      lati = output$receptor$lati, 
                                      long = output$receptor$long, 
                                      zagl = c(0, z_top))

    tmp_met_files <- find_met_files(output$receptor$run_time, met_file_format, 
                                    n_hours = tmp.nhrs, met_path)
    cat(paste('get.met.vars(): found meteo files of', tmp_met_files, '\n\n'))

    # ----------------------------------------------------------------------- #
    # generate trajectories to get specific humidity, temp, surface hgts/pressure, 
    # and u-/v-/w- wind speeds at the receptor location
    tmp.p <- calc_trajectory(namelist = tmp.namelist, 
                             rundir = dirname(output$file), emisshrs = 0.01, 
                             hnf_plume = T, met_files = tmp_met_files, 
                             n_hours = tmp.nhrs, output = tmp.output, rm_dat = T,
                             timeout = 600, w_option = 0, z_top = z_top) 
    if (is.null(tmp.p)) stop(paste('get.met.vars(): no trajec generated\n',
                             '*** Likely the current @param z_top of', z_top, 
                             'm is too high for the met field you adopted\n', 
                             '*** Lower z_top in config_xstilt.r'))

    # ----------------------------------------------------------------------- #
    # subset trajectories at min time step
    tmp.pmin <- tmp.p %>% filter(abs(time) == min(abs(time))) 

    # if specific humidity is all zero, meaning no SH is modeled in adopted met fields, 
    # use RH and saturation vapor pressure to calculate specific humidity 
    if (length(unique(tmp.pmin$sphu)) == 1) {           

        # Clausius–Clapeyron equation for saturation vapor pressure
        tmp.pmin <- tmp.pmin %>% mutate(
                    es = 6.112 * exp( (17.67 * (temz - 273.15)) / (temz - 273.15 + 243.5) ), 
                    e = rhfr * es,   # vapor pressure in air
                    w = 0.622 * e / (pres - e),   # water vapor mixing ratio
                    sphu = w / (1 + w)   # finally specific humidity in kg/kg
        ) %>% dplyr::select(-c(foot, mlht, dens, samt, sigw, tlgr, foot_no_hnf_dilution))       
    }   # end if


    # ----------------------------------------------------------------------- #
    # subset particles based on min release hgt, the most closed to the receptor
    tmp.sfc <- tmp.pmin %>% filter(xhgt == min(xhgt)) %>% rename(psfc = pres) %>% 
               dplyr::select(lati, long, time, zagl, zsfc, psfc, temp, xhgt)
    tmp.nsec <- abs(tmp.sfc$time) * 60      # in second
    
    # grab instantaneous variables
    # p1 or p2 in distCosine(), first one is longitude, second is latitude
    # distance in x- y- and z- directions, all in meters
    tmp.recp <- tmp.output$receptor[1, ]    # get receptor lat/lon

    library(geosphere)
    dx <- distCosine(p1 = c(tmp.sfc$long,  tmp.recp$lati),
                     p2 = c(tmp.recp$long, tmp.recp$lati)) # in meter
    dy <- distCosine(p1 = c(tmp.recp$long, tmp.sfc$lati),
                     p2 = c(tmp.recp$long, tmp.recp$lati))
    dz <- abs(tmp.sfc$zagl - tmp.sfc$xhgt)

    # calculate U-, V- and W- velocities [m/s]
    # if tmp.nhrs is negative (ie., -1), switch wind direction
    # dx always > 0, no direction, sign(long - long) provides direction on dx
    tmp.sfc <- tmp.sfc %>% mutate(
        ubar = sign(long - tmp.recp$long) * dx / tmp.nsec * sign(tmp.nhrs), 
        vbar = sign(lati - tmp.recp$lati) * dy / tmp.nsec * sign(tmp.nhrs), 
        wbar = sign(zagl - xhgt) * dz / tmp.nsec * sign(tmp.nhrs), 
        tsfc = temp - 273.15
    ) %>% dplyr::select(-c('lati', 'long', 'zagl', 'temp', 'xhgt', 'time'))
    

    # ----------------------------------------------------------------------- #
    # add surface wind/temp, pressure, height to 'receptor'
    # add temp and q profiles to 'output'
    output$receptor <- c(output$receptor, tmp.sfc)
    output$qt_prof <- tmp.pmin %>% dplyr::select(c('sphu', 'temz', 'pres')) %>% arrange(pres) 
    
    return(output)
}
