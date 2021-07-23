#' subroutines for estimating ground height before generating trajectories
#' @author: Dien Wu, 01/18/2019

#' "output" created from simulation_step() with lists of '$file' and '$receptor';
#' since before_trajec() will have the same environment as the env of simulation_step(), 
#' no need to add variables

before_trajec_xstilt = function() {

    # need modeled ground height to interpolate pres-hgt relation to
    # interpolate satellite weighting profiles from OCO-2 20 levels to model
    # levels, added by Dien Wu, 05/26/2018
    cat('before_trajec_xstilt(): extracting meteo vars from meteo fields...\n')

    # obtain multiple variables from met fields, DW, 09/14/2020
    # e.g., get specific humidity and temp profiles that will be used to 
    #   calculate dry-air column density in mol m-2 for further calculating PWF &
    #   convert XCO from mol m-2 to ppb
    output = get.met.vars(namelist, output, met_file_format, met_path, z_top)

    cat('END of before_trajec_xstilt(), start calculating X-trajectories...\n')
    return(output)
}  # end of function
