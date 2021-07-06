### subroutine to draw STILT particle locations on 2D or 3D plot
# need joinPolys() and point.in.polygon() functions from PBSmapping
# td for the threshold for cutting polygons
# by Dien Wu, 06/21/2017

# convert kernel density to percentile, DW, 11/15/2017
# return latitude range, DW, 11/16/2017
# use readRDS instead of getr (have changed in Trajecmulti), DW, 07/29/2018
# add customized data filtering, DW, 08/20/2018

# add dlat range for background range, DW, 08/24.2018
# add numbers of soundings used for background, DW, 09/05/2018
# add background uncertainty (including spread sd + retrieval err), DW, 09/07/2018
# bug fixed if there's no intersection between overpass and forward plume, DW, 10/03/2018
# minor update to the script in using OCO-3 data, DW, 06/30/2020
# add TROPOMI data so that we can define background for CO and NO2 
#   and grant update to the background code, DW, 03/04/2021

#' @param timestr can be in format of YYYYMMDD, but satellite data required!!
calc.bg.forward.trajec = function(site, timestr, sensor, sensor_path, sensor_gas, 
                                  sensor_ver, sensor_qa = 0.5, qfTF = T, 
                                  store_path, met, td = 0.1, bg_dlat = 0.5, 
                                  bg_dlon = 0.5, zoom = 8, api.key, 
                                  lon_lat = NULL, font.size = rel(1.0), 
                                  pp_fn = NULL){
  
  library(ggpubr); register_google(key = api.key)

  # figure name --
  traj_path = file.path(store_path, 'out_forward')
  plot_path = file.path(traj_path, 'plot'); dir.create(plot_path, showWarnings = F)

  cat(paste('\n\n# --------------------------------- #\n',
      'calc.bg.forward.trajec(): loading satellite obs and forward trajec for', 
      timestr, ', it takes time...\n'))

  # -------------------------------------------------------------------------- #
  # 1. call get.sensor.obs() to get satellite data, currently available include
  # OCO-2, OCO-3, TROPOMI CO, CH4, and NO2, DW, 07/06/2021
  if (is.null(lon_lat)) lon_lat = get.lon.lat(site, dlat = 1.5, dlon = 1.5)
  obs_df = get.sensor.obs(site, timestr, sensor, sensor_gas, sensor_fn = NULL,
                          sensor_path, sensor_ver, qfTF, tropomi_qa = sensor_qa, lon_lat)
  if ( is.null(obs_df) ) return()

  # revise column names for OCO data ------------------------------------------
  if ( grepl('OCO', sensor) ) 
    obs_df = obs_df %>% rename(val = xco2, val_uncert = xco2.uncert, time_utc = timestr) 
  
  # rename columne names for TROPOMI data -------------------------------------
  if ( grepl('TROPOMI', sensor)) {
    obs_df = obs_df %>% rename(lon = center_lon, lat = center_lat) %>% 
                        mutate(time_utc = substr(time_utc, 1, 10))
    if (sensor_gas == 'CO') obs_df = obs_df %>% rename(val = xco, val_uncert = xco_uncert)
    if (sensor_gas == 'NO2') obs_df = obs_df %>% rename(val = tno2, val_uncert = tno2_uncert)
    if (sensor_gas == 'CH4') obs_df = obs_df %>% rename(val = ch4_bc, val_uncert = ch4_uncert)

    # change to use orbit as ref, instead of timestr, DW, 05/25/2021
    # since TROPOMI does not have orbit like OCO, create fake 'orbit' based on timestr
    uni_time = obs_df %>% distinct(time_utc) %>% mutate(orbit = row_number())
    obs_df = obs_df %>% left_join(uni_time, by = 'time_utc')
  }  # end if

  # check to see IF multiple orbits are found for a single day, DW, 04/30/2021
  orbit = unique(obs_df$orbit)


  # ------------------------------------------------------------------------- #
  bg_df = NULL
  for ( tt in 1 : length(orbit) ) {

    if (length(orbit) > 1) 
      cat(paste('\n\nFound multiple orbits in a day => working on orbit #', orbit[tt], '\n'))

    # change to use orbit as ref, instead of timestr, DW, 05/25/2021
    #obs_tmp = obs_df %>% filter(time_utc == timestr[tt])
    obs_tmp = obs_df[obs_df$orbit == orbit[tt], ]
    timestr_tmp = min(obs_tmp$time_utc)

    # ----------------------------------------------------------------------- #
    # 2. call fit.kde.plume() to calculate 2D kernel density from particle distributions
    plm_list = fit.kde.plume(site, timestr = timestr_tmp, traj_path, 
                             obs_df = obs_tmp, sensor, td)
    if (is.null(plm_list) & is.null(bg_df)) return()
    if (is.null(plm_list) & !is.null(bg_df)) return(bg_df)

    # plm_list: 
    # recp_box  - df of four corner of receptor box
    # recp_info - df of receptor info (e.g., lat, lon, time, numpar, etc.)
    # sel_traj - df of trajectories during overpass duration
    # densf  - df of fitted 2D kernel density of the particles distribution
    # obs_tmp - df of modified obs with @param plmTF for those within urban plume 
    # plm_df - df of lat/lon coordinates of urban plume
    recp_box  = plm_list$recp_box
    recp_info = plm_list$recp_info
    sel_traj  = plm_list$sel_traj
    densf     = plm_list$densf
    obs_tmp   = plm_list$obs_df
    plm_df    = plm_list$plm_df
    td        = plm_list$td            # td may be updated
    intersectTF = plm_list$intersectTF
    site_lon = lon_lat$citylon
    site_lat = lon_lat$citylat


    # -------------------------------------------------------------------------- #
    # 3. if there is an intersect, calculate the background over upwind region
    if (intersectTF) {
      rds_fn = paste0('obs_plume_', site, '_', timestr_tmp, '_', sensor, 
                      '_', sensor_gas, '_qf', qfTF, '.rds')
      saveRDS(obs_tmp, file = file.path(traj_path, rds_fn)) 

      # get background values
      bg_tmp = calc.bg.upwind(site_lon, site_lat, obs_df = obs_tmp, sensor, 
                              bg.dlat = bg_dlat, bg.dlon = bg_dlon, perc = 0.1) 
      if (!is.null(bg_tmp)) bg_tmp = bg_tmp %>% mutate(timestr = timestr_tmp) %>% 
                                                relocate(timestr, .before = swath)                         
    } else bg_tmp = NULL


    # -------------------------------------------------------------------------- #
    # 4. call plot.urban.plume to generate figure 
    picname = file.path(plot_path, paste0('forward_plume_', site, '_', timestr_tmp, '_', 
                                          met, '_', sensor, '_', sensor_gas, '.png'))
    
    p1 = plot.bg(site, site_lon, site_lat, sensor, sensor_gas, recp_box, 
                 recp_info, sel_traj, densf, obs_df = obs_tmp, plm_df, 
                 intersectTF, bg_df = bg_tmp, bg_dlat, bg_dlon, zoom, td, 
                 picname, font.size, pp_fn) 
    
    # store all background info
    bg_df = rbind(bg_df, bg_tmp)
  } # end for 

  return(bg_df)
}  # end of subroutine

