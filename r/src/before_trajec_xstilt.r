#' subroutines for estimating ground height before generating trajectories
#' @author: Dien Wu, 01/18/2019

#' "output" created from simulation_step() with lists of '$file' and '$receptor';
#' since before_trajec() will have the same environment as the env of simulation_step(), 
#' no need to add variables

before_trajec_xstilt <- function() {

    # need modeled ground height to interpolate pres-hgt relation to
    # interpolate satellite weighting profiles from OCO-2 20 levels to model
    # levels, added by Dien Wu, 05/26/2018
    cat('before_trajec_xstilt(): estimating modeled ground heights for X-STILT ...\n')

    # store trajec from 5mAGL in the same copy dir 'rundir' by calling
    # get.ground.height() that calls calc_trajectory() to estimate ground
    # height [m] along w. u-, v- and w- component instantaneous wind
    # given receptor lat/lon/time/agl = 5 (near ground)
    recp.var <- get.ground.hgt(receptor = output$receptor, agl = 0, run_trajec, 
                               namelist, rundir, emisshrs, hnf_plume, 
                               met_file_format, met_loc, n_hours, rm_dat, 
                               timeout, w_option, z_top)
                                
    # paste interpolated info to output$receptor
    output$receptor <- c(output$receptor, recp.var)

    cat('END of before_trajec_xstilt(), start calculating X-trajectories...\n')
    return(output)
}  # end of function
