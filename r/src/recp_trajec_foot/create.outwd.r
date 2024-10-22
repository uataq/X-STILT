# read satellite data and craete output dir for TROPOMI or OCO-3, DW, 07/02/2021
# add TCCON, DW, 04/21/2023

# obs_* can be NA, for column simulations without satellite data
create.outwd = function(timestr, obs_species, obs_sensor, obs_path, obs_fn = NA,
                        lon_lat, store_path, met, run_hor_err) {
  
  if ( is.na(obs_sensor) ) {    # ideal run without satellite data
    obs_info = data.frame(timestr = timestr, fn = NA, stringsAsFactors = F)

    if (length(timestr) > 1) {
      dir_wd = paste0(min(timestr), '-', max(timestr))
    } else {
      dir_wd = ifelse(is.na(timestr), '',  timestr)
    }
  
    output_wd = file.path(store_path, 
                          paste0('out_', dir_wd, '_', met, '_ideal'))
    
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
    if (length(obs_fn) == 0) 
      stop('NO OCO-2/3 files found for this time, please check...\n')

    obs_info = data.frame(timestr = timestr, fn = obs_fn, stringsAsFactors = F)

    # path for storing trajec, foot for OCO-2/3, DW, 07/31/2018
    output_wd = file.path(store_path, paste('out', timestr, met, obs_sensor, sep = '_'))

  } else if ( grepl ('TCCON', obs_sensor) ) {
    
    if ( is.na(obs_fn) ) stop('create.outwd(): NO TCCON file assigned...\n')
    obs_info = data.frame(timestr = timestr, fn = obs_fn, stringsAsFactors = F)

    # path for storing trajec, foot for OCO-2/3, DW, 07/31/2018
    if ( length(timestr) > 1) {
      output_wd = file.path(store_path, 
                            paste('out', min(timestr), max(timestr), met, 
                                  obs_sensor, sep = '_'))
    } else output_wd = file.path(store_path, paste('out', timestr, met, obs_sensor, sep = '_'))

  } else stop("ONLY TROPOMI, OCO-2, OCO-3, NA are currently implemented, please check @param obs_sensor...\n")

  if (run_hor_err)  # for error output 
    output_wd = gsub(paste0('out_', timestr), 
                     paste0('outerr_', timestr), output_wd)

  cat(paste('\n\nOutput dir -', output_wd, '\n'))
  return(list(obs_info = obs_info, output_wd = output_wd))
}
