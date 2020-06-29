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

before_footprint_xstilt <- function() {

    # check whether weighted trajec exists already, DW, 07/13/2018
    #rundir <- dirname(output$file)
    wgt.file <- file.path(rundir, paste0(basename(rundir), '_wgttraj.rds'))

    if (file.exists(wgt.file) & !args$overwrite_wgttraj) {    
        wgt.output <- readRDS(wgt.file)     # grab from by-id directory
        
    } else {
       
        # get OCO-2/3 profile first according to lat/lon of receptor, 
        # return a list
        oco.info <- get.oco.info(oco.path = args$oco.path, receptor = output$receptor)

        if (is.null(oco.info)) {
            warning('before_footprint_xstilt(): 
                     NO OCO info found for this receptor lat/lon\n'); return()
        } # end if is.null

        # Weight footprint: call wgt.trajec.footv3() to weight trajec-level
        # footprint before calculating gridded footprint, DW, 06/01/2018
        cat("before_footprint_xstilt(): 
             weight trajec-level foot using OCO's profiles for X-STILT...\n")
        output$file <- wgt.file         # overwrite filename
        wgt.output  <- wgt.trajec.footv3(output = output, oco.info = oco.info,
                                         ak.wgt = args$ak.wgt, 
                                         pwf.wgt = args$pwf.wgt)
    }  # end if file.exists()


    ### if horizontal trans error is turned on, DW, 01/29/2019
    # calculate the trajec-level value, e.g., CO2 (FxEmiss), or others (FxGRIDs)
    if (args$run_hor_err) {
        cat('before_footprint_xstilt(): run_hor_err = T; estimating trajec-level value...\n')

        ## calculate trajec-level CO2 (emiss x foot) concentration using trajec 
        # with/without wind errors
        stat.file <- cal.trajfoot.stat(workdir = stilt_wd, output = output, 
                                       outdir = output_wd, met = args$met, 
                                       emiss.file = args$emiss.file, 
                                       combine.prof = wgt.output$wgt.prof, 
                                       ct.ver = args$ct.ver, 
                                       ctflux.path = args$ctflux.path, 
                                       ctmole.path = args$ctmole.path, 
                                       r_run_time, r_lati, r_long, r_zagl)
    } # end if run_hor_err

    # if foot_nhrs is different from trajec_nhrs, subset trajec
    if ( args$foot.nhrs < min(wgt.output$particle$time) / 60 ) {
        cat('before_footprint_xstilt(): subset particles for calculating footprint\n')
        wgt.output$particle <- wgt.output$particle %>% arrange(abs(time)) %>% 
                               filter(abs(time) <= abs(args$foot.nhrs) * 60) 
    }   # end if

    # we would like to generate footprint with different resolutions, if needed
    # DW, 02/11/2019 
    xres2 <- unlist(args$xres2)
    yres2 <- unlist(args$yres2)
    if ( !(NA %in% xres2) & !(NA %in% yres2) ) {

        cat('before_footprint_xstilt(): generate footprint with diff res...\n')
        for (f in 1: length(xres2)) {
            foot_file <- file.path(rundir, paste0(basename(rundir), '_', 
                                                  signif(xres2[f], 3), 'x', 
                                                  signif(yres2[f], 3), 
                                                  '_foot.nc'))
            
            # use weighted output for footprint
            foot <- calc_footprint(wgt.output$particle, output = foot_file,
                                   r_run_time = r_run_time,
                                   smooth_factor = smooth_factor,
                                   time_integrate = time_integrate,
                                   xmn = xmn, xmx = xmx, 
                                   xres = as.numeric(xres2[f]),
                                   ymn = ymn, ymx = ymx, 
                                   yres = as.numeric(yres2[f]))
        } # end for f
    } # end if is na
    #stop()
    cat('before_footprint_xstilt(): END OF FUNCTION, start to calculate foot...\n')
    return(wgt.output) # return weighted output
}
