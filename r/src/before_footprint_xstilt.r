#' function for weight trajec before calculating footprint, Dien Wu, 01/18/2019
#' @param output will be grabbed from simulation_step() automatically
#                with lists of '$file', '$receptor' and '$particle'
#' @param rundir where to store the weighted trajec, in 'by-id' directory
#' @param args a list of all the customized variables 
#'             (args is defined in simulation_step), no need to define it here
#' @param oco.path path for OCO-2/3 data, needed to frab AK, PW vertical profiles
#' @param ak.wgt logical flag for AK weighting, T or F 
#' @param pwf.wgt logical flag for PW weighting, T or F 

# allow for generating footprint with various horizontal resolutions, DW, 02/11/2019 
# minor update for using OCO-3 data, i.e., change variable names, DW, 06/28/2020
# add vertical AK weighting of TROPOMI column CO, DW, 08/25/2020 
# if ak.wgt == FALSE, remove dependence of any satellite, DW, 09/15/2020 
# add weighting of trajec-level footprints for NO2 modling by NOx lifetime, DW, 10/06/2020

before_footprint_xstilt <- function() {

    # after weighting the foot column, new foot with AK PW weighting is referred to as 'foot', 
    # the foot before weighting will be changed to 'foot_no_wgt', DW, 07/03/2020
    #if ('foot_before_weight' %in% colnames(output$particle) & args$overwrite.wgt == F) {    
    #    cat("before_footprint_xstilt(): 'foot_before_weight' column exists, no need to weight foot...\n")
    
    #} else {

    # ---------------------------------------------------------------------
    # weight trajec-level footprints using TROPOMI profiles, DW, 09/05/2020
    # ---------------------------------------------------------------------
    # whether to calc footprint for TROPOMI
    if ( !NA %in% unlist(args$tropomi.speci) ) {   

        cat('before_footprint_xstilt(): need to generate footprints based on TROPOMI profiles\n')
        p.tropomi <- wgt.trajec.foot.tropomi(output = output, 
                                             tropomi.fn = unlist(args$tropomi.fn), 
                                             tropomi.speci = unlist(args$tropomi.speci), 
                                             ak.wgt = args$ak.wgt, 
                                             pwf.wgt = args$pwf.wgt)

        # generate spatial footprint based on TROPOMI profiles 
        for ( tmp.speci in unlist(args$tropomi.speci) ) {
            
            cat(paste('\n\nbefore_footprint_xstilt(): calculating spatial footprint for TROPOMI', 
                        tmp.speci, '\n'))

            if (tmp.speci == 'CO') 
                p.tmp <- p.tropomi %>% rename(foot_before_wgt = foot, foot = foot_wgt_co)

            if (tmp.speci == 'CH4') 
                p.tmp <- p.tropomi %>% rename(foot_before_wgt = foot, foot = foot_wgt_ch4)

            # ---------------------------------------------------------------------
            # weight trajec-level footprints by NOx lifetime, DW, 10/06/2020
            # ---------------------------------------------------------------------
            if (tmp.speci == 'NO2') {
                
                p.prep <- p.tropomi %>% rename(foot_before_wgt = foot, 
                                               foot = foot_wgt_no2)
                
                # look for corresponding CTM file 
                ctm.fn <- find_ctm_files(r_run_time, n_hours, 
                                         ctm_file_format = args$ctm_file_format, 
                                         ctm_path = args$ctm_path)
                                         
                if (length(ctm.fn) == 0) {
                    cat('NO CTM files found, will skip NO2 modeling'); next }

                # performing footprint weighting using estimated K and modeled [OH]
                p.tmp <- wgt_foot_oh(p = p.prep, r_run_time, ctm = args$ctm_name, 
                                     ctm.fn = ctm.fn, xstilt_wd) 
                
                #cat('before_footprint_xstilt(): Saving chem-adjusted trajec...\n')
                #saveRDS(p.tmp, gsub('X_traj.rds', 'X_chem_traj.rds', output$file))
            }   # end if

            # use weighted output for footprint
            tmp.file <- file.path(rundir, paste0(basename(rundir), '_foot_TROPOMI_', tmp.speci, '.nc'))
            calc_footprint(p = p.tmp, output = tmp.file, r_run_time = r_run_time,
                           smooth_factor = smooth_factor, 
                           time_integrate = time_integrate, 
                           xmn = xmn, xmx = xmx, xres = xres, 
                           ymn = ymn, ymx = ymx, yres = yres)

            gc(); cat('\n')
        }   # end for tropomi speci
        
        cat('DONE with foot for TROPOMI..\n')
        if ( args$tropomiTF ) 
            stop('***THIS IS NOT AN ERROR, but no need to calculate foot for OCO, so stop the model...\n')
    }   # end if is.na()

    # ---------------------------------------------------------------------
    # weight trajec-level footprints by OCO-2 weighting profiles
    # this will be carried out by default
    # ---------------------------------------------------------------------
    # Weight footprint: call wgt.trajec.footv3() to weight trajec-level
    # footprint before calculating gridded footprint, DW, 06/01/2018
    cat("\n\nbefore_footprint_xstilt(): weight trajec-level foot using OCO's vertical profiles...\n")
    
    # for ideal simulation (ak.wgt == FALSE) then args$oco.path == NA
    # for real OCO simulation (ak.wgt == TRUE), will use AK from OCO files
    output <- wgt.trajec.foot(output, oco.path = args$oco.path, 
                              ak.wgt = args$ak.wgt, pwf.wgt = args$pwf.wgt)
    
   
    # -------------------------------------------------------------------------
    # calculate horizontal trans error if needed, DW, 01/29/2019
    # -------------------------------------------------------------------------
    # calculate the trajec-level value, e.g., CO2 (FxEmiss), or others (FxGRIDs)
    if ( args$run_hor_err ) {
        cat('\n\nbefore_footprint_xstilt(): run_hor_err = T; estimating trajec-level value...\n')

        ## calculate trajec-level CO2 (emiss x foot) concentration using trajec 
        # with/without wind errors
        stat.file <- cal.trajfoot.stat(workdir = stilt_wd, output = output, 
                                       outdir = output_wd, met = args$met, 
                                       emiss.file = args$emiss.file, 
                                       combine.prof = output$wgt.prof, 
                                       ct.ver = args$ct.ver, 
                                       ctflux.path = args$ctflux.path, 
                                       ctmole.path = args$ctmole.path, 
                                       r_run_time, r_lati, r_long, r_zagl)
    } # end if run_hor_err


    # -------------------------------------------------------------------------
    # if foot_nhrs is different from trajec_nhrs, subset trajec
    if ( args$foot_nhrs < n_hours ) {
        cat('\n\nbefore_footprint_xstilt(): subset particles for calculating footprint\n')
        output$particle <- output$particle %>% arrange(abs(time)) %>% 
                           filter(abs(time) <= abs(args$foot_nhrs) * 60) 
    }   # end if


    # -------------------------------------------------------------------------
    # we also need to generate footprint with different resolutions, if needed
    # -------------------------------------------------------------------------
    # DW, 02/11/2019 
    xres2 <- unlist(args$xres2)
    yres2 <- unlist(args$yres2)
    if ( !(NA %in% xres2) & !(NA %in% yres2) ) {

        cat('\n\nbefore_footprint_xstilt(): generate footprint with diff res...\n')
        for (f in 1 : length(xres2)) {
            foot_file <- file.path(rundir, paste0(basename(rundir), '_', 
                                                  signif(xres2[f], 3), 'x', 
                                                  signif(yres2[f], 3), 
                                                  '_foot.nc'))
            
            # use weighted output for footprint
            foot <- calc_footprint(p = output$particle, output = foot_file,
                                   r_run_time = r_run_time,
                                   smooth_factor = smooth_factor,
                                   time_integrate = time_integrate,
                                   xmn = xmn, xmx = xmx, 
                                   xres = as.numeric(xres2[f]),
                                   ymn = ymn, ymx = ymx, 
                                   yres = as.numeric(yres2[f]))
        } # end for f
    } # end if is na

    cat('END of before_footprint_xstilt(), start calculating gridded foot...\n')
    return(output) # return weighted output
}
