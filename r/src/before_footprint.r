#' subroutines for weight trajec before calculating footprint
#' @author: Dien Wu, 01/18/2019

#' @param output created from simulation_step() 
#                with lists of '$file', '$receptor' and '$particle'
#' @param rundir where to store the weighted trajec, in 'by-id' directory
#' @param oco2.path path for OCO-2 data, needed to frab AK, PW vertical profiles
#' @param ak.wgt logical flag for AK weighting, T or F 
#' @param pwf.wgt logical flag for PW weighting, T or F 

before_footprint <- function(output, rundir, oco2.path, ak.wgt = T, pwf.wgt = T) {

    # check whether weighted trajec exists already, DW, 07/13/2018
    wgt.file <- file.path(rundir, paste0(basename(rundir), '_wgttraj.rds'))
    cat(paste0('wgttraj file:', wgt.file, '\n'))

    if (file.exists(wgt.file)) {    # if true, grab from by-id directory
        wgt.output <- readRDS(wgt.file)
        
    } else {
        # get OCO-2 profile first according to lat/lon of receptor, return a list
        oco2.info <- get.oco2.info(oco2.path, receptor = output$receptor)

        if (is.null(oco2.info)) {
            warning('before_footprint(): NO OCO-2 info found for this receptor lat/lon')
            return()
        } # end if is.null

        # Weight footprint: call wgt.trajec.footv3() to weight trajec-level
        # footprint before calculating gridded footprint, DW, 06/01/2018
        cat('before_footprint(): weight foot column from trajec for X-STILT...\n')
        output$file <- wgt.file         # overwrite filename
        wgt.output  <- wgt.trajec.footv3(output = output, 
                                         oco2.info = oco2.info,
                                         ak.wgt = ak.wgt, 
                                         pwf.wgt = pwf.wgt)
    }  # end if file.exists()

    # overwrite 'particle' with weighted trajec
    particle <- wgt.output$particle
    return(particle)
}  # end of function
