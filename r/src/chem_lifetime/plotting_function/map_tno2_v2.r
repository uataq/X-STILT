
map_tno2_v2 = function(site, pr_df, emiss_df, zoom, labTF, ncol, 
                    mnz = NULL, mxz = NULL, title = NULL, 
                    xmn = NULL, xmx = NULL, ymn = NULL, ymx = NULL) {
    
    library(RColorBrewer); library(viridis)
    co2TF = FALSE; if ( 'obs_xco2_ff' %in% colnames(pr_df) ) co2TF = T 
    colnms = c('corner', 'lons', 'lats', 'lon', 'lat')
    pm_df = pr_df %>% reshape2::melt(id.vars = c('site', 'timestr', 'emiss', 
                                              'met', 'corner', 'lons', 'lats', 
                                              'lon', 'lat', 'polygon', 
                                              'obs_tno2_v1', 'obs_tno2_v2', 
                                              'obs_tno2_uncert_v1', 
                                              'obs_tno2_uncert_v2'))

    if (is.null(mnz)) mnz = max(min(pm_df$value, na.rm = T), -0.1)
    if (is.null(mxz)) {
        mxz = max(c(pm_df$value), na.rm = T)
        if (mxz < 0.4) mxz = 0.4; if (mxz > 8) mxz = 8
    }

    if (mxz < 1) dx = 0.1; if (mxz >= 1) dx = 0.2; if (mxz > 1.5) dx = 0.3
    if (mxz > 2) dx = 0.5; if (mxz > 5) dx = 1

    m0 = ggplot.map(map = 'ggmap', zoom = zoom, center.lat = mean(pr_df$lat),
                    center.lon = mean(pr_df$lon))[[1]] + #coord_equal(1.2) +
        labs(x = NULL, y = NULL, title = title) + coord_map() +
        xlim(c(min(pm_df$lons) - 0.05, max(pm_df$lons) + 0.05)) + 
        ylim(c(min(pm_df$lats) - 0.05, max(pm_df$lats) + 0.05))

    # 1. spatial map of modeled and observed columns -------------------------
    fll = def.col()[-12] #rev(brewer.pal(11, 'RdYlBu')
    fil_name = 'tNO2 [ppb]'
    new_var = init_var = sort(unique(pm_df$variable))
    names(new_var) = init_var

    #pm_df = pm_df %>% filter(variable == 'HRRR_EPA_TNO2_MIX')
    if (co2TF) {
        pm_df2 = pm_df %>% filter(variable == 'OBS_XCO2_FF')
        pm_df = pm_df %>% filter(variable != 'OBS_XCO2_FF')
    }

    m0 = m0 + theme_bw() + guides(color = guide_legend(ncol = 3)) + 
        facet_wrap(~variable, ncol = ncol, 
                    labeller = labeller(variable = new_var)) + 
        theme(legend.position = 'bottom', legend.key.height = unit(0.4, 'cm'),
            legend.key.width = unit(1.1, 'cm')) 
        
    m1 = m0 + scale_fill_gradientn(colors = fll, limits = c(mnz, mxz), 
                            breaks = seq(0, 10, dx), labels = seq(0, 10, dx), 
                            name = fil_name, na.value = '#2d0e4a') +
        geom_polygon(data = pm_df, aes(lons, lats, group = polygon, fill = value),
                    color = 'white', linewidth = 0.1, alpha = 0.9) +
        scale_size_continuous(range = c(4, 6)) #name = expr_dfession('ENOx [umol m'^-2~'s'^-1~']'))
    
    if (!is.null(xmn) & !is.null(xmx) & !is.null(ymn) & !is.null(ymx)) 
        m1 = m1 + annotate('rect', xmin = xmn, xmax = xmx, ymin = ymn, 
                            ymax = ymx, color = 'gray80', size = 0.4, 
                            fill = 'gray80', alpha = 0.1, linetype = 2)

    if (labTF) {
        lab_df = pm_df %>% group_by(variable) %>% 
                        dplyr::summarise(mean = mean(value, na.rm = T)) 
        m1 = m1 + geom_label(data = lab_df, aes(min(pm_df$lons) + 0.15, 
                            max(pm_df$lats) - 0.05, label = signif(mean, 2)), alpha = 0.5)
    }


     return(m12)
}
