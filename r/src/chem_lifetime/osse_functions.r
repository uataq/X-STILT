# -----------------------------------------------------------------------------
grab_wrfchem_ver = function(fn, max_level, ext, maxagl = 2000) {

    df = NULL
    
    pblh = load_wrfout_var(fn, 'PBLH', level = 1, ext)  # PBLH in meter
    hgt = load_wrfout_var(fn, 'HGT', level = 1, ext)   # Terrain Height in m

    for (level in 1 : max_level) {
        
        print(level)
        no2 = load_wrfout_var(wrfout_fns = fn, name = 'no2', level, ext)
        no  = load_wrfout_var(wrfout_fns = fn, name = 'no', level, ext)

        # grab total pressure = perturbation pres + state pres in Pa -----------
        pp = load_wrfout_var(fn, 'P', level, ext)
        pb = load_wrfout_var(fn, 'PB', level, ext)
        pres = (pp + pb) / 100        # convert to mb

        # compute altitudes, perturbation + base-state geopotential in m2 s-2
        ph = load_wrfout_var(fn, 'PH', level, ext)
        phb = load_wrfout_var(fn, 'PHB', level, ext)
        zasl = (ph + phb) / 9.81
        zagl = zasl - hgt   
        zagl[zagl < 0] = 0      # negative due to rounding error
        
        # merge all rasterLayers
        stk = stack(no2, no, pres, zasl, zagl, pblh, hgt)
        names(stk) = c('no2', 'no', 'pres', 'zasl', 'zagl', 'pblh', 'zsfc')

        tmp_utc = unique(names(no2))
        tmp_df = as.data.frame(stk, xy = T) %>% 
                 mutate(nox = no + no2, level = level, 
                        datestr = gsub('X', '',  tmp_utc), 
                        date = as.POSIXct(datestr, 'UTC', 
                                          format ='%Y.%m.%d.%H.%M.%S'))
        
        df = rbind(df, tmp_df)
        if (maxValue(zagl) > maxagl) break
    }
    
    return(df)
}

# -----------------------------------------------------------------------------
calc_tno2_wrfchem = function(obs_df, wrf_df, dd) {

    obs_sel = obs_df %>% filter(substr(datestr, 1, 13) == dd)
    wrf_sel = wrf_df %>% filter(substr(date, 1, 13) == dd, level <= 11,
                                x >= min(obs_sel$lons), 
                                x <= max(obs_sel$lons), 
                                y >= min(obs_sel$lats), 
                                y <= max(obs_sel$lats)) %>% 
              group_by(x, y) %>% mutate(min_pres = min(pres)) %>% ungroup() %>% 
              dplyr::select(x, y, date, no2, no, nox, pres, level, 
                            datestr, min_pres)

    # calculate the tropospheric column of NO2
    # for each WRF grid volume, find the nearest obs for tropospheric column
    wrf_add = NULL
    for (rr in 1 : nrow(wrf_sel)) {

        tmp_rr = wrf_sel[rr, ]
        if (rr %% 30 == 0) 
            print(paste('# ----', signif(rr / nrow(wrf_sel) * 1E2, 3), 
                        '% ---- #'))

        for (p in unique(obs_sel$polygon)) {
            obs_p = obs_sel %>% filter(polygon == p)
            pipTF = point.in.polygon(tmp_rr$x, tmp_rr$y, obs_p$lons, obs_p$lats)
            
            if (pipTF == 0) next 
            if (pipTF > 0) {    # if find the nearest TROPOMI sounding ...
                obs_u = obs_p %>% dplyr::select(polygon, tno2, psfc, 
                                                tropo_lower_pres, 
                                                tropo_upper_pres, air_vcd,
                                                air_vcd_tropo) %>% unique()
                tmp_rr = cbind(tmp_rr, obs_u)
                break
            }   # end if
        }   # end for p 

        if (rr > 1) { if (ncol(tmp_rr) != ncol(wrf_add)) next }
        wrf_add = rbind(wrf_add, tmp_rr)
    }   # end for rr


    # for modeled tropospheric column 
    wrf_avg = wrf_add %>% group_by(x, y) %>% summarise_if(is.numeric, mean) %>% 
              mutate(date = unique(wrf_sel$date), 
                     datestr = unique(wrf_sel$datestr), 
                     sf = (psfc - min_pres) / (psfc - tropo_lower_pres), 
                     level = 'X', no2 = no2 * sf, nox = nox * sf) %>% 
              ungroup() %>% dplyr::select(-sf)

    wrf_add = wrf_add %>% mutate(level = formatC(level, width = 2, flag = 0))
    wrf_com = full_join(wrf_add, wrf_avg, by = colnames(wrf_add))
    rds_list = list(obs_df = obs_sel, wrf_df = wrf_com)

    return(rds_list)
}



# -----------------------------------------------------------------------------
# try to use sf::st_join
calc_tno2_wrfchemv2 = function(obs_df, wrf_df, dd) {

    library(sf); library(sfheaders)

    obs_sel = obs_df %>% filter(substr(datestr, 1, 13) == dd)
    wrf_sel = wrf_df %>% filter(substr(date, 1, 13) == dd, level <= 11,
                                x >= min(obs_sel$lons) + 0.05, 
                                x <= max(obs_sel$lons) - 0.05, 
                                y >= min(obs_sel$lats) + 0.05, 
                                y <= max(obs_sel$lats) - 0.05) %>% 
              group_by(x, y) %>% mutate(min_pres = min(pres)) %>% ungroup() %>% 
              dplyr::select(x, y, date, no2, no, nox, pres, level, 
                            datestr, min_pres)

    # convert data.frame to sf object 
    wrf_sf = st_as_sf(x = wrf_sel, coords = c('x', 'y'), crs = 26918)
    obs_sf = sf_polygon(obj = obs_sel, x = 'lons', y = 'lats',           
                        polygon_id = 'polygon')     # sf object with polygon
    sf::st_crs(obs_sf) = 26918
    
    # label wrf grids with their corresponding TROPOMI polygon id
    # polygon may be NA, because it cannot find a TROPOMI grid for a WRF point
    wrf_add = st_join(wrf_sf, obs_sf) %>% arrange(polygon) %>% 
              left_join(obs_sel %>% dplyr::select(polygon, tno2, psfc, 
                                                  tropo_lower_pres, 
                                                  tropo_upper_pres, 
                                                  air_vcd, air_vcd_tropo), 
                        by = 'polygon') %>% sf_to_df(fill = T) %>% 
               dplyr::select(-c(sfg_id, point_id)) %>% unique() %>% na.omit()

    # for modeled tropospheric column 
    #xy_df = wrf_add %>% dplyr::select(x, y) %>% group_by(x, y) %>% tally()
    wrf_x = wrf_add %>% group_by(x, y) %>%    # x,y are wrf centered lat/lon
            summarise_if(is.numeric, mean) %>% 
            mutate(date = unique(wrf_sel$date), 
                   datestr = unique(wrf_sel$datestr), 
                   sf = (psfc - min_pres) / (psfc - tropo_lower_pres), 
                   no = no * sf, no2 = no2 * sf, nox = nox * sf) %>% 
            ungroup() %>% 
            left_join(obs_sel %>% dplyr::select(lons, lats, polygon), 
                      by = 'polygon')

    return(wrf_x)
}


# -----------------------------------------------------------------------------
plot_mod_obs = function(obs_x, wrf_x, stilt_x = NULL, stilt_x2 = NULL, 
                        emiss_df = NULL, m0, dd, site, 
                        title2 = 'STILT regrid') {

    # select data
    fn = paste0('wrfchem_tropomi_no2_', site,'_', timestr, '.png')
    obs_p = obs_x %>% filter(lons >= min(m0$data$lon), 
                             lons <= max(m0$data$lon), 
                             lats >= min(m0$data$lat), 
                             lats <= max(m0$data$lat))

    wrf_p = wrf_x %>% filter(x >= min(m0$data$lon), 
                             x <= max(m0$data$lon), 
                             y >= min(m0$data$lat), 
                             y <= max(m0$data$lat))
    mx = max(c(wrf_p$no2, obs_p$tno2))
    if (!is.null(stilt_x)) {
        stilt_p = stilt_x %>% filter(x >= min(m0$data$lon), 
                                     x <= max(m0$data$lon), 
                                     y >= min(m0$data$lat), 
                                     y <= max(m0$data$lat))
        mx = max(c(wrf_p$no2, obs_p$tno2, stilt_p$tno2))
        fn = paste0('wrfchem_stilt_tropomi_no2_', site,'_', timestr, '.png')
    }

    if (!is.null(stilt_x2)) {
        stilt_p2 = stilt_x2 %>% filter(x >= min(m0$data$lon), 
                                       x <= max(m0$data$lon), 
                                       y >= min(m0$data$lat), 
                                       y <= max(m0$data$lat))
        mx = max(c(wrf_p$no2, obs_p$tno2, stilt_p$tno2, stilt_p2$tno2))
    }


    dz = 0.5; if (mx < 2) dz = 0.2; if (mx > 5) dz = 1
    labs = seq(0, 10, dz)
    col = def.col()[-c(1, length(def.col()))]

    title = paste('Modeled vs. observed [NO2] in ppb on', dd, 'over', site)
    m1 = m0 + coord_equal(1.2) + labs(x = 'LONGITUDE', y = 'LATITUDE') + 
         scale_fill_gradientn(colors = col, name = NULL) + 
         theme(legend.key.height = unit(0.3, 'cm'), 
               legend.key.width = unit(1.5, 'cm'))
        
    o1 = m1 + ggtitle('TROPOMI') + theme_bw() + 
         geom_polygon(data = obs_p, 
                      aes(lons, lats, group = polygon, fill = tno2), 
                      color = 'gray70', alpha = 0.7, size = 0.2) + 
         geom_point(data = emiss_df, aes(lon, lat), color = 'white',
                    size = 3, shape = 4) 


    n1 = m1 + ggtitle('WRF-chem') + theme_bw() + 
         geom_tile(data = wrf_p, aes(x, y, fill = no2), 
                   color = 'gray70', alpha = 0.7) + 
         geom_point(data = emiss_df, aes(lon, lat), color = 'white',
                    size = 3, shape = 4) 

    width = 7; height = 4
    on = ggarrange(o1, n1, ncol = 2) 
    if (!is.null(stilt_x)) {
        s1 = m1 + ggtitle('STILT') + theme_bw() + 
             geom_polygon(data = stilt_p, aes(x, y, group = polygon, 
                                              fill = tno2), 
                          color = 'gray70', alpha = 0.7, size = 0.2) + 
             geom_point(data = emiss_df, aes(lon, lat), 
                        color = 'white', size = 3, shape = 4) 
        
        on = ggarrange(o1, n1, s1, ncol = 3); width = 10

        if (!is.null(stilt_x2)) {
            
            s2 = m1 + ggtitle(title2) + theme_bw() + 
                 geom_tile(data = stilt_p2, aes(x, y, fill = layer), 
                           color = 'gray70', alpha = 0.7)
            on = ggarrange(o1, n1, s1, s2, ncol = 2, nrow = 2)
            height = 9
        }       
    }

    on = annotate_figure(on, top = title)
    timestr = gsub('-', '', substr(dd, 1, 10))
    ggsave(on, filename = fn, width = width, height = height)

}


plot_wrf_stilt_obs = function(stilt_x, lon_lat, wrf_df1, wrf_df2, avg_x, 
                              w0, mx = NA) {

    # for plotting WRF-chem and TROPOMI ---------------------------------------
    library(RColorBrewer)
    w0 = w0 + scale_fill_gradientn(colors = hcl.colors(10), name = NULL) + 
              xlim(c(min(stilt_x$lons), max(stilt_x$lons))) + 
              ylim(c(min(stilt_x$lats), max(stilt_x$lats)))

    if (!is.na(mx)) 
        w0 = w0 + scale_fill_gradientn(colors = hcl.colors(10), 
                                        name = NULL, limits = c(0, mx))

    w2 = w0 + ggtitle('WRF-chem 12km tNO2 [ppb]') + 
        geom_tile(data = wrf_df2, aes(x, y, fill = no2, group = polygon), 
                  color = 'gray70', alpha = 0.8) + 
        geom_point(data = lon_lat, aes(citylon, citylat), 
                    color = 'white', shape = 4, size = 3) + 
        annotate('text', lon_lat$citylon - 0.25, lon_lat$citylat + 0.4, 
                 label = signif(avg_x$sim_tno2_wrf_12km, 3), size = 6, 
                 color = 'white')

    w1 = w0 + ggtitle('TROPOMI tNO2 [ppb]') + guides(alpha = 'none') + 
        scale_alpha_manual(name = NULL, values = c(0.6, 0.95)) + 
        #scale_fill_gradientn(colors = hcl.colors(10), name = NULL, 
        #                    limits = c(0, max(wrf_df1$tno2))) +
        geom_polygon(data = stilt_x, aes(lons, lats,fill = obs_tno2, 
        #geom_polygon(data = wrf_df1, aes(lons, lats, fill = tno2, 
                     group = indx, 
                     alpha = sim_tno2_ppb > avg_x$obs_tno2_uncert),
                     color = 'gray70', size = 0.2) +
        geom_point(data = lon_lat, aes(citylon, citylat), 
                    color = 'white', shape = 4, size = 3) + 
        annotate('text', lon_lat$citylon - 0.25, lon_lat$citylat + 0.4, 
                label = signif(avg_x$obs_tno2, 3), size = 6, color = 'white')

    w3 = w0 + ggtitle('WRF-chem 3km tNO2 [ppb]') + 
        geom_tile(data = wrf_df1, aes(x, y, fill = no2), 
                  color = 'gray70', alpha = 0.9, size = 0.2) + 
        geom_point(data = lon_lat, aes(citylon, citylat), 
                    color = 'white', shape = 4, size = 3) + 
        annotate('text', lon_lat$citylon - 0.25, lon_lat$citylat + 0.4, 
                  label = signif(avg_x$sim_tno2_wrf_3km, 3), size = 6, 
                color = 'white')
    ww = ggarrange(w1, w2, w3, ncol = 3)

    # for plotting STILT runs
    s0 = w0 + scale_alpha_manual(name = NULL, values = c(0.6, 0.95))
                        #limits = c(0, max(stilt_x$sim_tno2_nochem_ppb))) 

    s1 = s0 + ggtitle('STILT 1km mixing tNO2 [ppb]') + guides(alpha = 'none') + 
         geom_polygon(data = stilt_x, aes(lons, lats, fill = sim_tno2_ppb, 
                      group = indx, 
                      alpha = sim_tno2_ppb > avg_x$obs_tno2_uncert), 
                      color = 'gray70', size = 0.2) +
         geom_point(data = lon_lat, aes(citylon, citylat), color = 'white', 
                    shape = 4, size = 3) + 
         annotate('text',lon_lat$citylon - 0.25, lon_lat$citylat + 0.4, 
                  label = signif(avg_x$sim_tno2_ppb, 3), 
                  size = 6, color = 'white')

    s2 = s0 + ggtitle('STILT non-mixing 10m tNO2 [ppb]') + 
         geom_polygon(data = stilt_x, aes(lons, lats, fill = sim_tno2_nomix_ppb,
                      group = indx, 
                      alpha = sim_tno2_ppb > avg_x$obs_tno2_uncert), 
                      color = 'gray70', size = 0.2) +
         geom_point(data = lon_lat, aes(citylon, citylat), color = 'white', 
                    shape = 4, size = 3) + guides(alpha = 'none') + 
         annotate('text', lon_lat$citylon - 0.25, lon_lat$citylat + 0.4, 
                  label = signif(avg_x$sim_tno2_nomix_ppb, 3), size = 6, 
                  color = 'white')

    s3 = s0 + guides(alpha = 'none') + 
         ggtitle('STILT non-chem tNO2\nassumed NO2-NOx ratio of 0.74 [ppb]') +
         geom_polygon(data = stilt_x, aes(lons, lats, 
                            fill = sim_tno2_nochem_ppb, group = indx, 
                            alpha = sim_tno2_ppb > avg_x$obs_tno2_uncert), 
                            color = 'gray70', size = 0.2) +
        geom_point(data = lon_lat, aes(citylon, citylat), color = 'white', 
                    shape = 4, size = 3) + 
        annotate('text', lon_lat$citylon - 0.25, lon_lat$citylat + 0.4, 
                label = signif(avg_x$sim_tno2_nochem_ppb, 3), size = 6, 
                color = 'white')

    ss = ggarrange(s1, s2, s3, nrow = 1)
    ws = ggarrange(ww, ss, nrow = 2)

    return(ws)
}