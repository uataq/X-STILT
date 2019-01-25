#' subroutines for weight trajec before calculating footprint
#' @author: Dien Wu, 01/18/2019

#' @param output will be grabbed from simulation_step() automatically
#                with lists of '$file', '$receptor' and '$particle'
#' @param rundir where to store the weighted trajec, in 'by-id' directory
#' @param oco2.path path for OCO-2 data, needed to frab AK, PW vertical profiles
#' @param ak.wgt logical flag for AK weighting, T or F 
#' @param pwf.wgt logical flag for PW weighting, T or F 
#' environment of this script (including arguments, ak.wgt, pwf.wgt, oco2.path) 
#'     will be passed to simulation_step. 

generate_before_footprint_xstilt <- function(ak.wgt, pwf.wgt, oco2.path) {

    before_footprint_xstilt <- function() {

        # check whether weighted trajec exists already, DW, 07/13/2018
        rundir <- dirname(output$file)
        wgt.file <- file.path(rundir, paste0(basename(rundir), '_wgttraj.rds'))

        if (file.exists(wgt.file)) {    # if true, grab from by-id directory
            wgt.output <- readRDS(wgt.file)
            
        } else {
            # get OCO-2 profile first according to lat/lon of receptor, 
            # return a list
            oco2.info <- get.oco2.info(oco2.path = oco2.path, 
                                       receptor = output$receptor)

            if (is.null(oco2.info)) {
                warning('before_footprint_xstilt(): 
                NO OCO-2 info found for this receptor lat/lon\n')
                return()
            } # end if is.null

            # Weight footprint: call wgt.trajec.footv3() to weight trajec-level
            # footprint before calculating gridded footprint, DW, 06/01/2018
            cat('before_footprint_xstilt(): weight foot column from trajec for X-STILT...\n')
            output$file <- wgt.file         # overwrite filename
            wgt.output  <- wgt.trajec.footv3(output = output, 
                                             oco2.info = oco2.info,
                                             ak.wgt = ak.wgt, 
                                             pwf.wgt = pwf.wgt)
        }  # end if file.exists()

        return(wgt.output) # return weighted output
    }

    return(before_footprint_xstilt)  # return before_foot function
}  # end of function
