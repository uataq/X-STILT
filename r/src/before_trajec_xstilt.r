#' subroutines for estimating ground height before generating trajectories
#' @author: Dien Wu, 01/18/2019

#' "output" created from simulation_step() with lists of '$file' and '$receptor';
#' 

before_trajec_xstilt <- function(output, r_zagl = 5, run_trajec, varsiwant, 
                                 conage, cpack, dxf, dyf, dzf, emisshrs, frhmax,
                                 frhs, frme, frmr, frts, frvs, hnf_plume, hscale,
                                 ichem, iconvect, initd, isot, kbls, kblt, kdef,
                                 khmax, kmix0, kmixd, kmsl, kpuff, krnd, kspl, 
                                 kzmix, maxdim, maxpar, met_file_format, met_loc,
                                 mgmin, ncycl, ndump, ninit, n_hours, outdt, 
                                 outfrac, p10f, qcycle, random, splitf,
                                 tkerd, tkern, rm_dat = T, rundir, timeout, 
                                 tlfrac, tratio, tvmix, veght, vscale, w_option,
                                 zicontroltf, z_top) {

    # need modeled ground height to interpolate pres-hgt relation to
    # interpolate satellite weighting profiles from OCO-2 20 levels to model
    # levels, added by Dien Wu, 05/26/2018
    cat('before_trajec(): estimating modeled ground heights for X-STILT ...\n')

    # store trajec from 5mAGL in the same copy dir 'rundir' by calling
    # get.ground.height() that calls calc_trajectory() to estimate ground
    # height [m] along w. u-, v- and w- component instantaneous wind
    # given receptor lat/lon/time/agl=5 (near ground)
    # remove ziscale and zicontroltf from get.ground.hgt(), DW
    recp.var <- get.ground.hgt(receptor = output$receptor, r_zagl, run_trajec, 
                               varsiwant, conage, cpack, dxf, dyf, dzf, emisshrs, 
                               frhmax, frhs, frme, frmr, frts, frvs, hnf_plume, 
                               hscale, ichem, iconvect, initd, isot, kbls, kblt, 
                               kdef, khmax, kmix0, kmixd, kmsl, kpuff, krnd, kspl, 
                               kzmix, maxdim, maxpar, met_file_format, met_loc,
                               mgmin, ncycl, ndump, ninit, n_hours, outdt, outfrac, 
                               p10f, qcycle, random, splitf, tkerd, tkern, 
                               rm_dat, rundir, timeout, tlfrac, tratio, tvmix, 
                               veght, vscale, w_option, zicontroltf, z_top)
                                
    # paste interpolated info to output$receptor
    output$receptor <- c(output$receptor, recp.var)
    return(output)

}  # end of function
