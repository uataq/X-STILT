# script to calculate the background based on observations and wind direction 
# DW, 03/05/2021

if (F) {    # for debugging
    obs_df = obs_tmp     
    bg_deg = 0.2
    perc = 0.1
    bin_deg = 0.3
    plotTF = T
}

#' @param perc % of buffer adding outside the edge of the urban plume 
#' @param rm_outlierTF for removing potential FF signals
calc.bg.upwind = function(site_lon, site_lat, obs_df, sensor, sensor_gas, 
                          bg_deg = 0.3, perc = 0.1, bin_deg = 0.3, 
                          rm_outlierTF = T, plotTF = T, plot_path, map) {
    
    # check if there are any obs outside the plume, if not return NULL
    obs_out = obs_df %>% filter(plmTF == FALSE)
    if (nrow(obs_out) == 0) { cat('NO obs within background region...\n'); return()} 

    # will remove bg observations outside 2-sigma or only within [2.3%tile, 97.7%tile]
    if (rm_outlierTF) {
        library(cowplot)
        cat('calc.bg.upwind(): removing outlier observations from background calc...\n')

        # calculate the empirical cumulative density function of background obs
        Fn = ecdf(obs_out$val); ecdf_df = NULL

        # obs spacing that mod can be resolved, e.g., 0.1 ppm for XCO2, 1 ppb for XCO
        if (sensor_gas == 'CO2') by_obs = 0.1   
        if (sensor_gas != 'CO2') by_obs = 1
        for (test in seq(min(obs_out$val), max(obs_out$val), by_obs)) 
            ecdf_df = rbind(ecdf_df, data.frame(x = test, y = Fn(test)))

        edge = min(ecdf_df[ecdf_df$y >= 0.9, 'x'])
        #edge = quantile(obs_out$val, c(2.3, 97.7)/100)
        #edge = mean(bb$val) + sd(bb$val) * 2 * c(-1, 1)
        #edge = median(bb$val) + sd(bb$val) * 2 * c(-1, 1)

        # sanity check on the distribution of background obs, DW, 07/12/2021 
        if (plotTF) {
            o1 = map + scale_fill_gradientn(colours = def.col()) + 
                 scale_alpha_manual(values = c(0.5, 0.9)) +
                 geom_polygon(data = obs_df, aes(lons, lats, group = polygon, 
                                                 fill = val, alpha = plmTF)) 
                 
            if (sensor_gas != 'CO2') 
                o1 = o1 + geom_text(data = obs_df, aes(lon, lat, label = signif(val, 2)), 
                                    size = 1.5, color = 'white') 

            h1 = ggplot() + theme_classic() + labs(x = 'Xobs', y = 'Count') + 
                 geom_histogram(data = obs_df, aes(val, fill = plmTF), color = 'gray30')
            d1 = ggplot() + theme_classic() + labs(x = 'Xobs_bg', y = NULL, title = 'Cum Emp Dens') + 
                 geom_point(data = ecdf_df, aes(x, y), size = 0.4) + 
                 geom_vline(xintercept = edge, color = 'blue')

            xloc = 0.5; yloc = 0.6; if (sensor_gas == 'CO2') {xloc = 0.6; yloc = 0.7}
            hd = ggdraw() + draw_plot(h1) + draw_plot(d1, x = xloc, y = yloc, width = .3, height = .3)
            picname = file.path(plot_path, paste0('background_demo_x', tolower(sensor_gas), 
                                                  '_', min(unique(obs_df$time_utc)), '.png'))
            ggsave(ggarrange(o1, hd), filename = picname, width = 11, height = 5)
        }

        obs_out = obs_out %>% filter(val <= edge)
        obs_df = rbind(obs_df %>% filter(plmTF), obs_out)
    }   # end if


    ### calculate multiple background based on 'swath' column in obs
    # multiple swaths for SAMs; only one swath for non-SAM, DW, 08/04/2020
    # ------------------------------------------------------------------------ #
    # select soundings to the background side BUT exclude those in the obs_plm
    nn = select.obs.side(obs_df, sensor, 'north', bg_deg, perc, bin_deg)$obs_bg
    if ( nrow(nn) > 0 ) nn_df = get.bg.inner(nn, 'north') else nn_df = NULL

    # ------------------------------------------------------------------------ #
    ss = select.obs.side(obs_df, sensor, 'south', bg_deg, perc, bin_deg)$obs_bg
    if ( nrow(ss) > 0 ) ss_df = get.bg.inner(ss, 'south') else ss_df = NULL

    # ------------------------------------------------------------------------ #
    ee = select.obs.side(obs_df, sensor, 'east', bg_deg, perc, bin_deg)$obs_bg
    if ( nrow(ee) > 0 ) ee_df = get.bg.inner(ee, 'east') else ee_df = NULL 

    # ------------------------------------------------------------------------ #
    ww = select.obs.side(obs_df, sensor, 'west', bg_deg, perc, bin_deg)$obs_bg
    if ( nrow(ww) > 0 ) ww_df = get.bg.inner(ww, 'west') else ww_df = NULL 
    
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
select.obs.side = function(obs_df, sensor, bg.side, bg_deg, perc = 0.1, 
                           bin_deg = 0.3) {
    
    # for northern or southern bg, bin up longitudes --------------------------
    if ( bg.side %in% c('north', 'south') & sensor == 'TROPOMI' ) {
        bin_x = seq(floor(min(obs_df$lon)), ceiling(max(obs_df$lon)), bin_deg)
        obs_df = obs_df %>% arrange(lon) %>% 
                 mutate(swath = findInterval(lon, bin_x), bin = bin_x[swath])

    } else if ( bg.side %in% c('north', 'south') & grepl('OCO', sensor) ) {

        # if OCO, bin as the min lon for each swath
        obs_df = obs_df %>% group_by(swath) %>% mutate(bin = min(lon))
    }  # end if binning long


    # for eastern or western bg, bin up latitudes ---------------------------- 
    if ( bg.side %in% c('east', 'west') & sensor == 'TROPOMI' ) {
        bin_y = seq(floor(min(obs_df$lat)), ceiling(max(obs_df$lat)), bin_deg)
        obs_df = obs_df %>% arrange(lat) %>% 
                 mutate(swath = findInterval(lat, bin_y), bin = bin_y[swath])

    }  else if ( bg.side %in% c('east', 'west') & grepl('OCO', sensor) ) { 

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
        bb = obs_out %>% group_by(swath, bin) %>% filter(lat >= plm.ymx & lat <= plm.ymx + bg_deg)

    if (bg.side == 'south') 
        bb = obs_out %>% group_by(swath, bin) %>% filter(lat >= plm.ymn - bg_deg & lat <= plm.ymn)

    if (bg.side == 'east') 
        bb = obs_out %>% group_by(swath, bin) %>% filter(lon >= plm.xmx & lon <= plm.xmx + bg_deg) 

    if (bg.side == 'west') 
        bb = obs_out %>% group_by(swath, bin) %>% filter(lon >= plm.xmn - bg_deg & lon <= plm.xmn)

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
                             bg.mean = mean(val), 
                             bg.median = median(val), 
                             bg.sd.spread = sd(val), 
                             bg.sd.retrv = sqrt(mean(val_uncert ^ 2)), 
                             .groups = 'drop') %>% 
            left_join(bb %>% group_by(swath, bin) %>% 
            tally(), by = c('swath', 'bin')) %>% 
            mutate(bg.side = bg.side)
    
    return(bg_df)
}

# end of functions 
