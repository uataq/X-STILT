# read satellite data and craete output dir for TROPOMI or OCO-3, DW, 07/02/2021


# obs_* can be NA, for column simulations without satellite data
create.outwd = function(timestr, obs_species, obs_sensor, obs_path, lon_lat, 
                        store_path, met, run_hor_err) {

  if ( is.na(obs_sensor) ) {    # ideal run without satellite data
    obs_info = data.frame(timestr = timestr, fn = NA, stringsAsFactors = F)
    output_wd = file.path(store_path, paste('out', timestr, met, 'ideal', sep = '_'))

  } else if ( obs_sensor == 'TROPOMI' ) {
    
    # get satellite info (need to look for the correct nc file for TROPOMI) ----
    cat(paste('\n\n --- Obtaining TROPOMI', obs_species, 'info --- \n'))
    obs_info = find.tropomi(tropomi_path = obs_path, timestr, lon_lat) 
    if (is.null(obs_info)) stop('NO TROPOMI overpass found...terminating...\n')
    obs_info = obs_info %>% mutate(species = obs_species)
    
    # get TROPOMI overpass hour near the site of interest ----
    timestr = unique(substr(obs_info$overpass.start.time, 1, 10))

    # path for storing trajec, foot for TROPOMI, DW, 07/05/2021
    output_wd = file.path(store_path, 
                          paste('out', timestr, met, obs_sensor, sep = '_'), 
                          obs_species)

  } else if ( grepl('OCO', obs_sensor) ) {
    
    cat(paste('\n\n --- Obtaining', obs_sensor, obs_species, 'info --- \n'))
    obs_fn = list.files(obs_path, paste0('_', substr(timestr, 3, 8), '_'), 
                        full.names = T, recursive = T)
    obs_info = data.frame(timestr = timestr, fn = obs_fn, stringsAsFactors = F)

    # path for storing trajec, foot for OCO-2/3, DW, 07/31/2018
    output_wd = file.path(store_path, paste('out', timestr, met, obs_sensor, sep = '_'))

  } else stop("ONLY TROPOMI, OCO-2, OCO-3, NA are currently implemented, please check @param obs_sensor...\n")

  if (run_hor_err)  # for error output 
    output_wd = gsub(paste0('out_', timestr), paste0('outerr_', timestr), output_wd)


  cat(paste('\n\nOutput dir - ', output_wd, '\n'))
  return(list(obs_info = obs_info, output_wd = output_wd))
}
