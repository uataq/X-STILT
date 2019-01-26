#' subroutines for weight trajec before calculating footprint
#' @author: Dien Wu, 01/18/2019

#' @param output will be grabbed from simulation_step() automatically
#                with lists of '$file', '$receptor' and '$particle'
#' @param rundir where to store the weighted trajec, in 'by-id' directory

#' @param args a list of all the customized variables 
#'             (args is defined in simulation_step), no need to define it here
#'      @param oco2.path path for OCO-2 data, needed to frab AK, PW vertical profiles
#'      @param ak.wgt logical flag for AK weighting, T or F 
#'      @param pwf.wgt logical flag for PW weighting, T or F 

before_footprint_xstilt <- function() {

    # check whether weighted trajec exists already, DW, 07/13/2018
    #rundir <- dirname(output$file)
    wgt.file <- file.path(rundir, paste0(basename(rundir), '_wgttraj.rds'))

    if (file.exists(wgt.file)) {    # if true, grab from by-id directory
        wgt.output <- readRDS(wgt.file)
        
    } else {
        # get OCO-2 profile first according to lat/lon of receptor, 
        # return a list
        oco2.info <- get.oco2.info(oco2.path = args$oco2.path, 
                                   receptor = output$receptor)

        if (is.null(oco2.info)) {
            warning('before_footprint_xstilt(): 
                     NO OCO-2 info found for this receptor lat/lon\n')
            return()
        } # end if is.null

        # Weight footprint: call wgt.trajec.footv3() to weight trajec-level
        # footprint before calculating gridded footprint, DW, 06/01/2018
        cat('before_footprint_xstilt(): 
             weight trajec-level foot using OCO-2 profiles for X-STILT...\n')
        output$file <- wgt.file         # overwrite filename
        wgt.output  <- wgt.trajec.footv3(output = output, 
                                         oco2.info = oco2.info,
                                         ak.wgt = args$ak.wgt, 
                                         pwf.wgt = args$pwf.wgt)
    }  # end if file.exists()

    cat('before_footprint_xstilt(): 
         return weighted trajec and start to calculate foot if possible...\n')
    return(wgt.output) # return weighted output
}
