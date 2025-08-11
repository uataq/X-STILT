# function to obtain specific humidity and temp profiles from adopted met fields 
# DW, 08/26/2020 

# get specific humidity and temp profiles that will be used to calculate 
#   dry-air column density in mol m-2 for further calculating PWF and 
#   converting XCO from mol m-2 to ppb


if (F) {
    met_file_format = '%Y%m%d'
    met_path = file.path(homedir, 'gfs0p25')
    output$receptor = list(run_time = r_run_time, lati = r_lati,
                            long = r_long, zagl = r_zagl)
}

# still need specific humidity profiles to calculate c.dry per layer 
# returns 'output'
get.met.vars = function(namelist, output, met_file_format, met_file_tres, 
                        met_path, z_top = 25000, run_slant = FALSE) {

    tmp_namelist = namelist
    tmp_namelist$delt = 1
    
    #' @param nturb is NotTURBulence flag that turns turbulence on (FALSE or 0) or off (TRUE or 1), DW, 2022/06/02
    #' because we would like to interpolate the wind vectors from met fields, 
    #' we turn off the turbulence for this interpolation, but not for STILT runs
    tmp_namelist$nturb = T
    tmp_namelist$numpar = z_top / 10
    tmp_namelist$varsiwant = c('time', 'indx', 'long', 'lati', 'zagl', 'zsfc', 
                               'foot', 'mlht', 'dens', 'samt', 'sigw', 'tlgr', 
                               'temp', 'pres', 'sphu', 'rhfr', 'temz')

    # write tmp_output
    tmp_nhrs = -1
    tmp_output = list()
    tmp_output$file = gsub('X_traj.rds', 'met_vars.rds', output$file)

    # if nadir column
    if (!run_slant) {
        tmp_output$receptor = data.frame(run_time = output$receptor$run_time, 
                                         lati = output$receptor$lati, 
                                         long = output$receptor$long, 
                                         zagl = c(0, z_top))
    } else {
        # if slant column, use mean receptor long/lati for simplifications 
        tmp_output$receptor = data.frame(
            run_time = output$receptor$run_time, 
            lati = mean(unlist(output$receptor$lati)),
            long = mean(unlist(output$receptor$long)),
            zagl = c(0, z_top)
        )
    }
    tmp_met_files = find_met_files(output$receptor$run_time, n_hours = tmp_nhrs,
                                   met_path, met_file_format, met_file_tres)
    cat(paste('get.met.vars(): found meteo files of', tmp_met_files, '\n\n'))

    # ----------------------------------------------------------------------- #
    # generate trajectories to get specific humidity, temp, surface hgts/pressure, 
    # and u-/v-/w- wind speeds at the receptor location
    tmp_p = calc_trajectory(namelist = tmp_namelist, 
                            rundir = dirname(output$file), 
                            emisshrs = 0.01, 
                            hnf_plume = F, 
                            met_files = tmp_met_files, 
                            n_hours = tmp_nhrs, 
                            output = tmp_output, 
                            rm_dat = T, 
                            timeout = 600, 
                            w_option = 0, 
                            z_top = z_top) 
                             
    if (is.null(tmp_p)) stop(paste('get.met.vars(): no trajec generated...check error message in <output_path>/by-id/<receptor>\npossible reasons: met fields not compatible (e.g., missing met files, met_path too long, etc)'))

    # ----------------------------------------------------------------------- #
    # subset trajectories at min time step
    tmp_pmin = tmp_p %>% filter(abs(time) == min(abs(time))) 

    # if specific humidity is all zero, meaning no SH is modeled in adopted met fields, 
    # use RH and saturation vapor pressure to calculate specific humidity 
    if (length(unique(tmp_pmin$sphu)) == 1) {           

        # Clausius–Clapeyron equation for saturation vapor pressure
        tmp_pmin = tmp_pmin %>% mutate(
                    es = 6.112 * exp( (17.67 * (temz - 273.15)) / (temz - 273.15 + 243.5) ), 
                    e = rhfr * es,   # vapor pressure in air
                    w = 0.622 * e / (pres - e),   # water vapor mixing ratio
                    sphu = w / (1 + w)   # finally specific humidity in kg/kg
        ) %>% dplyr::select(-c(foot, mlht, dens, samt, sigw, tlgr))       
    }   # end if


    # ----------------------------------------------------------------------- #
    # subset particles based on min release hgt, the most closed to the receptor
    tmp_sfc = tmp_pmin %>% filter(xhgt == min(xhgt)) %>% rename(psfc = pres) %>%
              dplyr::select(lati, long, time, zagl, zsfc, psfc, temp, xhgt)
    tmp_nsec = abs(tmp_sfc$time) * 60      # in second
    
    # grab instantaneous variables
    # p1 or p2 in distCosine(), first one is longitude, second is latitude
    # distance in x- y- and z- directions, all in meters
    tmp_rlong = tmp_output$receptor[1, 'long']    # get receptor lat/lon
    tmp_rlati = tmp_output$receptor[1, 'lati']

    library(geosphere)
    dx = distCosine(p1 = c(tmp_sfc$long, tmp_rlati),
                    p2 = c(tmp_rlong, tmp_rlati))       # in meter
    dy = distCosine(p1 = c(tmp_rlong, tmp_sfc$lati),
                    p2 = c(tmp_rlong, tmp_rlati))
    dz = abs(tmp_sfc$zagl - tmp_sfc$xhgt)

    # calculate U-, V- and W- velocities [m/s]
    # if tmp_nhrs is negative (ie., -1), switch wind direction
    # dx always > 0, no direction, sign(long - long) provides direction on dx
    tmp_sfc = tmp_sfc %>% mutate(
        ubar = sign(long - tmp_rlong) * dx / tmp_nsec * sign(tmp_nhrs), 
        vbar = sign(lati - tmp_rlati) * dy / tmp_nsec * sign(tmp_nhrs), 
        wbar = sign(zagl - xhgt) * dz / tmp_nsec * sign(tmp_nhrs), 
        tsfc = temp - 273.15
    ) %>% dplyr::select(-c('lati', 'long', 'zagl', 'temp', 'xhgt', 'time'))
    

    # ----------------------------------------------------------------------- #
    # add surface wind/temp, pressure, height to 'receptor'
    # add temp and q profiles to 'output'
    output$receptor = c(output$receptor, tmp_sfc)
    output$qt_prof = tmp_pmin %>% dplyr::select(c('sphu', 'temz', 'pres')) %>% 
                     arrange(pres) 
    
    return(output)
}
