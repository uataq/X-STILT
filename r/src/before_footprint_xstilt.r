#' function for weight trajec before calculating footprint, Dien Wu, 01/18/2019
#' @param output will be grabbed from simulation_step() automatically
#                with lists of '$file', '$receptor' and '$particle'
#' @param rundir where to store the weighted trajec, in 'by-id' directory
#' @param args a list of all the customized variables 
#'             (args is defined in simulation_step), no need to define it here
#' @param ak_wgt logical flag for AK weighting, T or F 
#' @param pwf_wgt logical flag for PW weighting, T or F 

# allow for generating footprint with various horizontal resolutions, DW, 02/11/2019 
# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# add vertical AK weighting of TROPOMI column CO, DW, 08/25/2020 
# if ak_wgt == FALSE, remove dependence of any satellite, DW, 09/15/2020 
# add weighting of trajec-level footprints for NO2 modling by NOx lifetime, DW, 10/06/2020

before_footprint_xstilt = function() {

    # since xres2, yres2, and time_intergrate2 can contain multiple variables, 
    # we've already keep them in a list, now unlist 
    xres2 = unlist(args$xres2)
    yres2 = unlist(args$yres2)
    time_integrate2 = unlist(args$time_integrate2)
    obs_species = args$obs_species 
    obs_sensor = args$obs_sensor 
    obs_fn  = args$obs_fn 
    ak_wgt  = args$ak_wgt 
    pwf_wgt = args$pwf_wgt 

    # -------------------------------------------------------------------------
    # if foot_nhrs differs from trajec_nhrs, subset trajec
    if ( args$foot_nhrs < n_hours ) {

        cat('\n\nbefore_footprint_xstilt(): subset particles for calculating footprint\n')
        output$particle = output$particle %>% arrange(abs(time)) %>% 
                          filter(abs(time) <= abs(args$foot_nhrs) * 60) 
    }   # end if


    # ---------------------------------------------------------------------
    # weight trajec-level footprints using sensor-specific profiles, DW, 09/05/2020
    # ---------------------------------------------------------------------
    cat('before_footprint_xstilt(): footprint weighting ')
    if ( grepl('OCO', obs_sensor) | is.na(obs_sensor)) {

        # ---------------------------------------------------------------------
        # weight trajec-level footprints using OCO weighting profiles (ak_wgt = T)
        # or set AK = 1 for ideal case (ak_wgt = F)
        # ---------------------------------------------------------------------
        if (ak_wgt) cat('based on OCO profiles at the nearest sounding\n')
        if (!ak_wgt) cat('assuming AK = 1 (no satellite dependencies)\n')
        output = wgt.trajec.foot(output = output, oco.fn = obs_fn, 
                                 ak.wgt = ak_wgt, pwf.wgt = pwf_wgt) 

    } else if (obs_sensor == 'TROPOMI') {
        
        # ---------------------------------------------------------------------
        # weight trajec-level footprints using TROPOMI weighting profiles (ak_wgt = T)
        # ---------------------------------------------------------------------
        cat('based on TROPOMI profiles at the nearest sounding\n')
        output = wgt.trajec.foot.tropomi(output = output, tropomi.fn = obs_fn, 
                                         tropomi.species = obs_species, 
                                         ak.wgt = ak_wgt, pwf.wgt = pwf_wgt) 
                
        # ---------------------------------------------------------------------
        # weight trajec-level footprints by NOx lifetime, DW, TBD/....
        # ---------------------------------------------------------------------
    } 
    
    cat('before_footprint_xstilt(): DONE with footprint weighting...\n')



    # ---------------------------------------------------------------------
    # calculate footprint with different configurations if there is 
    # ---------------------------------------------------------------------
    if ( !(NA %in% xres2) & !(NA %in% yres2) ) {

        cat('\n\nbefore_footprint_xstilt(): generate footprint with different config...\n')
        for (f in 1 : length(xres2)) {
            
            tmp_fn = file.path(rundir, paste0(basename(rundir), '_', 
                                              signif(xres2[f], 3), 'x', 
                                              signif(yres2[f], 3), '_foot.nc'))
            
            if (!time_integrate2[f])
                tmp_fn = file.path(rundir, paste0(basename(rundir), '_', 
                                                  signif(xres2[f], 3), 'x', 
                                                  signif(yres2[f], 3), 
                                                  '_hrly_foot.nc'))

            # use weighted output for footprint
            calc_footprint(p = output$particle, output = tmp_fn, 
                           r_run_time = r_run_time, smooth_factor = smooth_factor,
                           time_integrate = time_integrate2[f],
                           xmn = xmn, xmx = xmx, xres = as.numeric(xres2[f]),
                           ymn = ymn, ymx = ymx, yres = as.numeric(yres2[f]))
        } # end for f
    }   # end if is na

    cat('\nbefore_footprint_xstilt(): DONE with optional footprint calculation..\n')


    # -------------------------------------------------------------------------
    # calculate horizontal trans error if needed, DW, 01/29/2019
    # -------------------------------------------------------------------------
    # calculate the trajec-level value, e.g., CO2 (FxEmiss), or others (FxGRIDs)
    if ( args$run_hor_err ) {
        cat('\n\nbefore_footprint_xstilt(): run_hor_err = T; estimating trajec-level value...\n')

        ## calculate trajec-level CO2 (emiss x foot) concentration using trajec 
        # with/without wind errors
        stat.file = cal.trajfoot.stat(workdir = xstilt_wd, output = output, 
                                      outdir = output_wd, met = args$met, 
                                      emiss.file = args$emiss_fn, 
                                      combine.prof = output$wgt, 
                                      ct.ver = args$ct_ver, 
                                      ctflux.path = args$ctflux_path, 
                                      ctmole.path = args$ctmole_path, 
                                      r_run_time, r_lati, r_long, r_zagl)
    }   # end if run_hor_err

    cat('before_footprint_xstilt(): move onto footprint calculation with the default config...\n')
    return(output) # return weighted output
}
