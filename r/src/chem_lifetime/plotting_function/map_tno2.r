
map_tno2 = function(site, pr, emiss_df, zoom, labTF, ncol, 
                    mnz = NULL, mxz = NULL, title = NULL, 
                    xmn = NULL, xmx = NULL, ymn = NULL, ymx = NULL, 
                    facet_dir = c('v', 'h')[1]) {
     
     library(RColorBrewer); library(viridis)
     #if (is.null(title))
     #  title = paste('a) Spatial map of modeled vs. retrieved tNO2 over', site)
     
     co2TF = FALSE; if ( 'obs_xco2_ff' %in% colnames(pr) ) co2TF = T 
     colnms = c('corner', 'lons', 'lats', 'lon', 'lat')
     pm = pr %>% dplyr::select(-ends_with('_ratio')) %>% 
          reshape2::melt(id.vars = c(colnms, 'polygon')) %>% 
          mutate(variable = toupper(variable))

     if (is.null(mnz)) mnz = max(min(pm$value, na.rm = T), -0.1)
     if (is.null(mxz)) {
          mxz = max(c(pm$value), na.rm = T)
          if (mxz < 0.4) mxz = 0.4; if (mxz > 8) mxz = 8
     }
     
     if (mxz < 1) dx = 0.1; if (mxz >= 1) dx = 0.2; if (mxz > 1.5) dx = 0.3
     if (mxz > 2) dx = 0.5; if (mxz > 5) dx = 1

     m0 = ggplot.map(map = 'ggmap', zoom = zoom, center.lat = mean(pr$lat),
                     center.lon = mean(pr$lon))[[1]] + #coord_equal(1.2) +
          #ggplot() + 
          labs(x = NULL, y = NULL, title = title) + coord_map() +
          xlim(c(min(pm$lons) - 0.05, max(pm$lons) + 0.05)) + 
          ylim(c(min(pm$lats) - 0.05, max(pm$lats) + 0.05))
     
     emiss_df = emiss_df %>% filter(lon >= min(pm$lons) - 0.05, 
                                    lon <= max(pm$lons) + 0.05, 
                                    lat >= min(pm$lats) - 0.05, 
                                    lat <= max(pm$lats) + 0.05) %>% 
                mutate(enox = signif(enox, 3))
     print(emiss_df)

     # 1. spatial map of modeled and observed columns -------------------------
     fll = def.col()[-12] #rev(brewer.pal(11, 'RdYlBu')
     fil_name = 'tNO2 [ppb]'
     new_var = init_var = sort(unique(pm$variable))
     names(new_var) = init_var

     #pm = pm %>% filter(variable == 'HRRR_EPA_TNO2_MIX')
     if (co2TF) {
          pm2 = pm %>% filter(variable == 'OBS_XCO2_FF')
          pm = pm %>% filter(variable != 'OBS_XCO2_FF')
     }
     
     m0 = m0 + theme_bw() + guides(color = guide_legend(ncol = 3)) + 
          facet_wrap(~variable, ncol = ncol, dir = facet_dir, 
                     labeller = labeller(variable = new_var)) + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.4, 'cm'),
                legend.key.width = unit(1.4, 'cm'), 
                text = element_text(family = 'Arial')) 
     
     library(viridis)
     m1 = m0 + 
          scale_fill_gradientn(colors = fll, limits = c(mnz, mxz), 
          #scale_fill_viridis(option = 'magma', limits = c(mnz, mxz), 
                              breaks = seq(0, 10, dx), labels = seq(0, 10, dx), 
                              name = fil_name, na.value = 'black') +
          geom_polygon(data = pm, aes(lons, lats, group = polygon,fill = value),
                       color = 'white', linewidth = 0.04, alpha = 0.9) +
          geom_point(data = emiss_df, aes(lon, lat, size = enox), alpha = 0.4, 
                     color = 'white', shape = 19, show.legend = F) +
          scale_size_continuous(range = c(4, 6)) #name = expression('ENOx [umol m'^-2~'s'^-1~']'))
     
     if (!is.null(xmn) & !is.null(xmx) & !is.null(ymn) & !is.null(ymx)) 
          m1 = m1 + annotate('rect', xmin = xmn, xmax = xmx, ymin = ymn, 
                              ymax = ymx, color = 'gray80', size = 0.4, 
                              fill = 'gray80', alpha = 0.1, linetype = 2)

     if (labTF) {
          lab_df = pm %>% group_by(variable) %>% 
                          dplyr::summarise(mean = mean(value, na.rm = T)) 
          m1 = m1 + geom_label(data = lab_df, aes(min(pm$lons) + 0.15, 
                               max(pm$lats) - 0.05, label = signif(mean, 2)), alpha = 0.5)
     }

     if ( co2TF) {
          m2 = m0 + labs(title = '  ') + 
               theme(legend.position = 'bottom', 
                     legend.key.height = unit(0.4, 'cm'),
                     legend.key.width = unit(0.9, 'cm')) +
               scale_fill_gradientn(colors = fll, 
                                    limits = c(-1, max(pm2$value)), 
                                    breaks = seq(-2, 10, 0.5), 
                                    labels = seq(-2, 10, 0.5), 
                                    name = 'XCO2_ff\n[ppm]', 
                                    na.value = '#2d0e4a') +
               geom_polygon(data = pm2, aes(lons, lats, group = polygon,
                              fill = value), color = 'white', 
                              linewidth = 0.03, alpha = 0.9) +
               geom_point(data = emiss_df, aes(lon, lat, size = enox), 
                          alpha = 0.4, color = 'white', shape = 19, 
                          show.legend = F) +
               scale_size_continuous(range = c(4, 6)) 

          library(cowplot)
          m12 = plot_grid(m1, m2, nrow = 1, align = 'h',
                          rel_widths = c(ncol - 1, 1.1))
     } else m12 = m1

     return(m12)
}
