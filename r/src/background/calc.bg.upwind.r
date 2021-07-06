# script to calculate the background based on observations and wind direction 
# DW, 03/05/2021

if (F) obs_df = obs_tmp     # for debugging

#' @param perc % of buffer adding to the edge of the urban plume 
calc.bg.upwind = function(site_lon, site_lat, obs_df, sensor, bg.dlat = 0.5, 
                          bg.dlon = 0.5, perc = 0.1) {
    
    # obs within and outside the urban plume as indicated by the convex hull
    obs_out = obs_df %>% filter(plmTF == FALSE)     # outside the plume
    if (nrow(obs_out) == 0) { cat('NO obs within background region...\n'); return() } 
    
    ### calculate multiple background based on 'swath' column in obs
    # multiple swaths for SAMs; only one swath for non-SAM, DW, 08/04/2020
    # ------------------------------------------------------------------------ #
    # select soundings to the background side BUT exclude those in the obs_plm
    nn = select.obs.side(obs_df, sensor, 'north', bg.dlat, bg.dlon, perc)$obs_bg
    if ( nrow(nn) > 0 ) { nn_df = get.bg.inner(nn, 'north') } else nn_df = NULL

    # ------------------------------------------------------------------------ #
    ss = select.obs.side(obs_df, sensor, 'south', bg.dlat, bg.dlon, perc)$obs_bg
    if ( nrow(ss) > 0 ) { ss_df = get.bg.inner(ss, 'south') } else ss_df = NULL

    # ------------------------------------------------------------------------ #
    ee = select.obs.side(obs_df, sensor, 'east', bg.dlat, bg.dlon, perc)$obs_bg
    if ( nrow(ee) > 0 ) { ee_df = get.bg.inner(ee, 'east') } else ee_df = NULL 

    # ------------------------------------------------------------------------ #
    ww = select.obs.side(obs_df, sensor, 'west', bg.dlat, bg.dlon, perc)$obs_bg
    if ( nrow(ww) > 0 ) { ww_df = get.bg.inner(ww, 'west') } else ww_df = NULL 
    
    if (is.null(nn_df) & is.null(ww_df) & is.null(ee_df) & is.null(ss_df)) {
        cat('NO background obs found...\n'); return() }
    
    # otherwise return background info
    bg_df = rbind(nn_df, ss_df, ee_df, ww_df) %>% 
            mutate(bg.sd = sqrt(bg.sd.spread ^ 2 + bg.sd.retrv ^ 2)) 
            #left_join(plm_stat, by = 'swath')

    return(bg_df)
}



# --------------------------------------------------------------------------- #
# inner function to get plume statistics, e.g., min, max lat/lon
get.plm.stat = function(obs_df, perc) {
    
    # get urban enhanced lat range for each swath for SAM, DW, 08/04/2020
    # now get the boundary lat/lon of urban plume 
    plm_stat = obs_df %>% filter(plmTF == TRUE) %>% group_by(swath) %>% 
               dplyr::summarize(plm.xmn = min(lon), plm.xmx = max(lon),
                                plm.ymn = min(lat), plm.ymx = max(lat),
                                dy = plm.ymx - plm.ymn, dx = plm.xmx - plm.xmn, 

            # adjust enhanced lat range if buffer exists (i.e., dx * perc)
                                plm.ymn = plm.ymn - dy * perc,
                                plm.ymx = plm.ymx + dy * perc, 
                                plm.xmn = plm.xmn - dx * perc,
                                plm.xmx = plm.xmx + dx * perc)

    obs_df = obs_df %>% group_by(swath) %>% left_join(plm_stat, by = 'swath') %>% ungroup()
    return(obs_df)
}


# --------------------------------------------------------------------------- #
# inner function to subset observations
# if for TROPOMI (no swath available), create mini lat/lon bins
select.obs.side = function(obs_df, sensor, bg.side, bg.dlat, bg.dlon, perc = 0.1) {
    
    # for northern or southern bg, bin up longitudes
    if (bg.side %in% c('north', 'south') & sensor == 'TROPOMI') {
        bin_x = seq(floor(min(obs_df$lon)), ceiling(max(obs_df$lon)), 0.3)
        obs_df = obs_df %>% arrange(lon) %>% 
                 mutate(swath = findInterval(lon, bin_x), bin = bin_x[swath])

    } else if (grepl('OCO', sensor)) {

        # if OCO, bin as the min lon for each swath
        obs_df = obs_df %>% group_by(swath) %>% mutate(bin = min(lon))
    }  # end if binning long

    # for eastern or western bg, bin up latitudes
    if (bg.side %in% c('east', 'west') & sensor == 'TROPOMI') {
        bin_y = seq(floor(min(obs_df$lat)), ceiling(max(obs_df$lat)), 0.3)
        obs_df = obs_df %>% arrange(lat) %>% 
                 mutate(swath = findInterval(lat, bin_y), bin = bin_y[swath])

    }  else if (grepl('OCO', sensor)) { 

        # if OCO, bin as the min lat for each swath
        obs_df = obs_df %>% group_by(swath) %>% mutate(bin = min(lat))
    }  # end if binning lat


    # ------------------------------------------------------------------------ #
    # calculate the plume statistics and paste to obs_df
    obs_st  = get.plm.stat(obs_df, perc)
    obs_out = obs_st %>% filter(plmTF == FALSE)

    # select soundings to the background side BUT exclude those in the obs_plm
    # plm.xmn, plm.xmx, plm.ymn, plm.ymx can be NA 
    #   because no plume is found within certain swaths -> meaning all 'clean' obs
    #   thus, should also select obs with NA plm.* 
    if (bg.side == 'north') 
        bb = obs_out %>% group_by(swath, bin) %>% filter(lat >= plm.ymx & lat <= plm.ymx + bg.dlat)

    if (bg.side == 'south') 
        bb = obs_out %>% group_by(swath, bin) %>% filter(lat >= plm.ymn - bg.dlat & lat <= plm.ymn)

    if (bg.side == 'east') 
        bb = obs_out %>% group_by(swath, bin) %>% filter(lon >= plm.xmx & lon <= plm.xmx + bg.dlon) 

    if (bg.side == 'west') 
        bb = obs_out %>% group_by(swath, bin) %>% filter(lon >= plm.xmn - bg.dlon & lon <= plm.xmn)

    bb = bb %>% ungroup()

    return(list(obs_bg = bb, obs_df = obs_df))
}


# --------------------------------------------------------------------------- #
# include retrieval error in background uncert, DW, 09/06/2018
get.bg.inner = function(bb, bg.side) {
    
    if (nrow(bb) == 0) return(NULL)
    bg_df = bb %>% group_by(swath, bin) %>% 
            dplyr::summarize(bg.xmn = min(lon), bg.xmx = max(lon), 
                             bg.ymn = min(lat), bg.ymx = max(lat), 
                             bg.mean = mean(val), bg.median = median(val), 
                             bg.sd.spread = sd(val), 
                             bg.sd.retrv = sqrt(mean(val_uncert ^ 2)), 
                             .groups = 'drop') %>% 
            left_join(bb %>% group_by(swath, bin) %>% tally(), by = c('swath', 'bin')) %>% 
            mutate(bg.side = bg.side)
    
    return(bg_df)
}

