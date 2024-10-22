# subroutine to determine receptor locations based on obs soundings (centered lat/lon)
# DW, 07/02/2021 

# possible obs_filter --- 
# c('QF', 0): OCO data with QF of 0 
# c('QA', #): TROPOMI data with QA > #, e.g., 0.4/0.5/0.7/1 (see TROPOMI ATBD)
# NULL: NO data filtering, deal with all soundings

# randomly create receptors within the satellite polygon if jitterTF = T
# DW, 07/03/2021

# if combine_locTF is TRUE, requires OCO path, DW, 07/23/2021
# use OCO-3 soundings as receptors loc for TROPOMI, 
# but preserve TROPOMI overpass time for receptor time

# attempt to evenly select satellite soundings, DW, 06/17/2022 
# implement slant column release (recalc recp lon/lat based on SZA and SAA), DW, 04/21/2023

get.recp.sensorv2 = function(timestr, obs_filter = NULL, obs_fn, obs_sensor, 
                             obs_path, obs_spescies, lon_lat, jitterTF = FALSE, 
                             num_jitter, numpar, nf_dlat, nf_dlon, num_bg_lat, 
                             num_bg_lon, num_nf_lat, num_nf_lon, agl, run_trajec, run_slant = FALSE, output_wd, 
                             combine_locTF = FALSE, oco_path = NULL){
    
    # -------------------------------------------------------------------------
    # 1 - grab all available OCO or TROPOMI data first 
    # -------------------------------------------------------------------------
    qfTF = FALSE
    if (!is.null(obs_filter)) { 
        qfTF = T; tropomi_qa = as.numeric(obs_filter[2]) }

    # make sure selTF is correct based on num_*
    selTF = TRUE; if (is.na(num_bg_lat) | is.na(num_nf_lat)) selTF = FALSE
    obs_df = get.sensor.obs(site, timestr, sensor = obs_sensor, 
                            sensor_gas = obs_species, sensor_fn = obs_fn, 
                            sensor_path = obs_path, qfTF, tropomi_qa, lon_lat) 
    
    rename_lookup = c(lon = 'center_lon', lat = 'center_lat', 
                      vertices = 'corner')
    obs_df = obs_df %>% rename(any_of(rename_lookup))
    
    # remove vertices coordinate from obs
    if ( 'lons' %in% colnames(obs_df) ) {
        uni_df = unique(obs_df %>% dplyr::select(-c(lons, lats, vertices)))
    } else uni_df = obs_df %>% unique()
    

    # -------------------------------------------------------------------------
    # 2 - select observations
    # -------------------------------------------------------------------------
    output_fns = list.files(file.path(output_wd, 'by-id'), 'X_traj.rds', 
                            recursive = T, full.names = T)

    if ( length(output_fns) > 0 & !run_trajec ) {       # if trajec data exists

        cat('get.recp.sensorv2(): Found existing trajectories...\n')           
        recp_info = ident.to.info(ident = basename(output_fns),stilt.ver = 2)%>%
                    mutate(run_time = as.POSIXct(substr(timestr, 1, 12), 'UTC', 
                                                 format = '%Y%m%d%H%M')) %>% 
                    dplyr::select(lati = recp.lat, long = recp.lon, run_time)

    } else if ( combine_locTF == TRUE & obs_sensor == 'TROPOMI' ) {     

        # ----------------------------------------------------------------------
        # use OCO soundings as receptor locations and TROPOMI overpass time
        # for setting up receptors
        cat('get.recp.sensorv2(): since combine_locTF == TRUE, use OCO loc and TROPOMI time for receptors...\n')
        oco_df = get.sensor.obs(site, timestr, sensor = 'OCO-3', 
                                sensor_gas = 'CO2', sensor_fn = NULL, 
                                sensor_path = oco_path, qfTF, lon_lat = lon_lat)
        
        if (nrow(oco_df) == 0) 
            stop('get.recp.sensor(): NO OCO loc found for placing receptors with TROPOMI time')
        oco_df = unique(oco_df %>% dplyr::select(-c(lons, lats, vertices)))
        
        # place denser receptors within lat range with high XCO2
        if (selTF) oco_df = sel.obs4recpv2(df = oco_df, nf_dlat, nf_dlon, 
                                           lon_lat, num_bg_lat, num_bg_lon, 
                                           num_nf_lat, num_nf_lon)
        
        # replace OCO receptor time with TROPOMI time 
        trp_sel = uni_df %>% filter(lat >= min(oco_df$lat), 
                                    lat <= max(oco_df$lat), 
                                    lon >= min(oco_df$lon), 
                                    lon <= max(oco_df$lon))
        trp_time = seq(min(trp_sel$datestr), max(trp_sel$datestr), 
                       length = nrow(oco_df))
        recp_info = oco_df %>% arrange(lat) %>% mutate(datestr = trp_time)

    } else if ( selTF == TRUE ) {   
        
        # ----------------------------------------------------------------------
        # spatially select soundings as receptors; works for OCO and TROPOMI
        if (obs_sensor %in% c('OCO-2', 'OCO-3', 'TROPOMI')) {
            recp_info = sel.obs4recpv2(df = uni_df, nf_dlat, nf_dlon, lon_lat, 
                                       num_bg_lat, num_bg_lon, 
                                       num_nf_lat, num_nf_lon)

        } #else if (obs_sensor == 'TCCON') {
            # for TCCON, since it's continuous, if selTF is true, 
            # let's perform an hourly integration 
            
        #}

    } else if ( jitterTF & !is.na(num_jitter) ) {
        
        # ----------------------------------------------------------------------
        # create more recptors around centered lat/lon for a given sounding
        # works for either OCO and TROPOMI
        cat('get.recp.sensorv2(): jittering option turned on...\n')
        recp_info = jitter.obs4recp(obs_df, num_jitter)

    } else recp_info = uni_df        # use all obs around site
    
    
    # ----------------------------------------------------------------------
    # round lat, lon for each sounding, fix bug, DW, 07/31/2018
    if ( !'run_time' %in% colnames(recp_info) ) {

        # requires SZA and SAA info from observation or users
        if ( run_slant ) {
            
            if ( !'sza' %in% colnames(recp_info) | 
                 !'saa' %in% colnames(recp_info) ) 
            stop('You turned on slant column release @param run_slant, but no SZA/SAA found in the observation data for calc slant x-receptors...\n')
            
            recp_info = recp_info %>% dplyr::select(run_time = datestr, 
                                                    lati = lat, long = lon, 
                                                    sza, saa)

        } else recp_info = recp_info %>% dplyr::select(run_time = datestr, 
                                                       lati = lat, long = lon)

        recp_info = recp_info %>% 
                    mutate(lati = signif(lati, 6), long = signif(long, 7)) %>% 
                    arrange(lati)
    }   # end if

    ## add release height
    recp_info$zagl = list(agl)


    # ----------------------------------------------------------------------
    # now tilt the column receptor if slant column release is enabled 
    if (run_slant) {
       
        R = 6378137     # earth radius in m
        
        # work on each row and calc the lon/lat along slant column 
        # make lati/long a list   
        recp_slant = NULL
        for ( r in 1 : nrow(recp_info) ) {
          
          tmp = recp_info[r, ]
          xhgt_min = min(unlist(tmp$zagl))
          xhgt_max = max(unlist(tmp$zagl))
          xhgt_rng = xhgt_max - xhgt_min
          xhgt_step = xhgt_rng / numpar
          xhgt = (1:numpar) * xhgt_step + xhgt_min
          
          # expand fixed lon/lat at Nadir to slant lon/lat
          # create new data frame where lati, long, zagl are all list()
          dl_m = xhgt * tan(tmp$sza * pi / 180)

          # deviations in lon/lat per particle
          dlon_m = dl_m * sin(tmp$saa * pi / 180)
          dlat_m = dl_m * cos(tmp$saa * pi / 180)

          # calc new lon/lat along slant column 
          long_slant = tmp$long + dlon_m / R * 180 / pi
          lati_slant = tmp$lati + dlat_m / (R * cos(tmp$lati * pi / 180)) * 180 / pi

          tmp_slant = data.frame(run_time = tmp$run_time, lati0 = tmp$lati, 
                                 long0 = tmp$long)
          tmp_slant$lati = list(lati_slant)
          tmp_slant$long = list(long_slant)
          tmp_slant$zagl = list(xhgt)

          recp_slant = rbind(recp_slant, tmp_slant)
        }   # end for

        recp_info = recp_slant
    }        

    return(recp_info)
} # end of subroutine



# -------------------------------------------------------------------------
# inner function to select soundings within NEAR-FIELD 
# 
#           2*dlat by 2*dlon for ENTIRE region
# |----------------------------------------------------|
# |           x              x              x          |
# |           x              x              x          |
# |           x              x              x          |
# |         2*nf_dlat by 2*nf_dlon for NEAR-FIELD      |
# |           x     |---------------|       x          |
# |           x     | x x x x x x x |       x          |
# |           x     | x x x x x x x |       x          |
# |           x     | x x x x x x x |       x          |
# |           x     | x x x x x x x |       x          |
# |           x     | x x x x x x x |       x          |
# |           x     |---------------|       x          |
# |           x              x              x          |
# |           x              x              x          |
# |           x              x              x          |
# |           x              x              x          |
# |----------------------------------------------------| 

sel.obs4recpv2 = function(df, nf_dlat, nf_dlon, lon_lat, num_bg_lat, 
                          num_bg_lon, num_nf_lat, num_nf_lon) {

    # figure out the # of soundings per long or per lati based on num_nf
    nf_y = seq(lon_lat$site_lat - nf_dlat, lon_lat$site_lat + nf_dlat, 
               length = num_nf_lat)
    nf_x = seq(lon_lat$site_lon - nf_dlon, lon_lat$site_lon + nf_dlon, 
               length = num_nf_lon)
    bg_y = seq(lon_lat$minlat, lon_lat$maxlat, length = num_bg_lat)
    bg_x = seq(lon_lat$minlon, lon_lat$maxlon, length = num_bg_lon)

    nf_xy = expand.grid(x = nf_x, y = nf_y, fac = 'NF', stringsAsFactors = F)
    bg_xy = expand.grid(x = bg_x, y = bg_y, fac = 'BG', stringsAsFactors = F)
    xy = rbind(nf_xy, bg_xy) 

    # find the closest satellite soundings (based on their centered lat/lon)
    # loop over each grid and find out the closest satellite soundings
    cat('sel.obs4recpv2(): selecting satellite soundings...\n')
    for (r in 1: nrow(xy)) {

        tmp_xy = xy[r, ]
        tmp_dist = do.call(c, lapply(1 : nrow(df), function(rr){
            geosphere::distCosine(p1 = c(tmp_xy$x, tmp_xy$y), 
                                  p2 = c(df$lon[rr], df$lat[rr]))
        }))
        tmp_indx = which(tmp_dist == min(tmp_dist))
        tmp_polygon = df$polygon[tmp_indx]
        xy[r, 'polygon'] = tmp_polygon
    }

    recp_info = df %>% filter(polygon %in% unique(xy$polygon))
    return(recp_info)
}




if (F) {

    nf_dlat = 0.3; nf_dlon = 0.5; num_nf_lat = 7; num_nf_lon = 14
    num_bg_lon = num_bg_lat = 3

    tlt = paste('BACKGROUND: num_bg_lat =', num_bg_lat, 
                '; num_bg_lon =', num_bg_lon, 
                '\nNEARFIELD: num_nf_lat =', num_nf_lat, 
                'num_nf_lon =', num_nf_lon)
    t1 = ggplot(data = obs_df %>% 
                        mutate(selTF = polygon %in% unique(xy$polygon))) + 
         geom_polygon(aes(lons, lats, group = polygon, fill = tno2, 
                            alpha = selTF), color = 'black', linewidth = 0.1) +
         scale_fill_gradientn(colors = def.col()) + theme_dark() + 
         scale_alpha_manual(values = c(0.2, 0.9)) + 
         geom_point(data = xy, aes(x, y), color = 'white', shape = 8) +
         labs(title = tlt, x = 'LONGITUDE', y = 'LATITUDE')


}