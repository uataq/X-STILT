# subroutine to determine receptor locations based on obs soundings (centered lat/lon)
# DW, 07/02/2021 

# possible obs_filter --- 
# c('QF', 0): OCO data with QF of 0 
# c('QA', #): TROPOMI data with QA > #, e.g., 0.4/0.5/0.7/1 (see TROPOMI ATBD)
# NULL: NO data filtering, deal with all soundings

# randomly create receptors within the satellite polygon if jitterTF = T
# DW, 07/03/2021

get.recp.sensor = function(timestr, obs_filter = NULL, obs_fn, obs_sensor, 
                           obs_ver, obs_path, lon_lat, selTF = FALSE, 
                           jitterTF = FALSE, num_jitter, peak_lat, num_bg, 
                           num_peak, agl, run_trajec, output_wd = NULL){

  # grab OCO or TROPOMI data
  qfTF = FALSE
  if (!is.null(obs_filter)) { qfTF = TRUE; tropomi_qa = as.numeric(obs_filter[2]) } 
  obs_dat = get.sensor.obs(site, timestr, sensor = obs_sensor, 
                           sensor_gas = obs_species, sensor_fn = obs_fn, 
                           sensor_path = obs_path, sensor_ver = obs_ver, 
                           qfTF, tropomi_qa, lon_lat) 
  
  if (obs_sensor == 'TROPOMI') 
    obs_dat = obs_dat %>% rename(lon = center_lon, lat = center_lat, vertices = corner)

  # remove vertices cordinate 
  uni_dat = obs_dat %>% dplyr::select(-c('lons', 'lats', 'vertices')) %>% unique()


  # ------------------- Step 2. SET UP the STILT receptors ----------------- #
  output_fns = list.files(file.path(output_wd, 'by-id'), 'X_traj.rds', 
                          recursive = T, full.names = T)

  if ( length(output_fns) > 0 & !run_trajec ) {       # if trajec data exists

    cat('Found existing trajectories...\n')           
    recp_info = ident.to.info(ident = basename(output_fns), stilt.ver = 2) %>% 
                mutate(run_time = as.POSIXct(substr(timestr, 1, 12), 'UTC', 
                                             format = '%Y%m%d%H%M')) %>% 
                dplyr::select('lati' = 'recp.lat', 'long' = 'recp.lon', 'run_time')

  } else {    

    # -------------------------------------------------------------------------
    # if no existing trajec found, grab obs soundings as receptors
    if ( grepl('OCO', obs_sensor) & selTF ) {  # only for OCO mission 
      
      # place denser receptors within lat range with high XCO2
      recp_indx = c(seq(lon_lat$minlat, peak_lat[1],    1 / num_bg),
                    seq(peak_lat[1],    peak_lat[2],    1 / num_peak),
                    seq(peak_lat[1],    lon_lat$maxlat, 1 / num_bg))

      # select lat, lon based on OCO soundings
      uni.lat  = uni_dat$lat
      sort.lat = sort(uni.lat)
      recp.lat = sort.lat[findInterval(recp_indx, sort.lat)]
      match.index = unique(match(recp.lat, uni.lat))
      recp_info = uni_dat[match.index, ]

    } else if ( jitterTF & !is.na(num_jitter) ) {

      # call inner subroutine to create more recptors around centered lat/lon 
      # for a given satellite sounding
      recp_info = jitter.obs4recp(obs_dat, num_jitter)

    } else recp_info = obs_dat        # simply use all obs around site
    
    # round lat, lon for each sounding, fix bug, DW, 07/31/2018
    recp_info = recp_info %>% dplyr::select(run_time = datestr, lati = lat, 
                                            long = lon) %>% arrange(lati) %>% 
                mutate(lati = signif(lati, 6), long = signif(long, 7))
  } # end if trajec file existed


  ## add release height
  recp_info$zagl = list(agl)

  # return receptor info
  recp_info = recp_info %>% arrange(lati)
  return(recp_info)
} # end of subroutine
