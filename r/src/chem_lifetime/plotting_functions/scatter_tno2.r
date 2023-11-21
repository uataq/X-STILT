
scatter_tno2 = function(site, timestr, pr, lm_sma, mns = 0.05, mxs = NULL, 
                        title = NULL) {

     library(RColorBrewer)
     facs = c('01 - EDGAR (no chem, mix)', '02 - EDGAR (chem, mix)', 
              '03 - EDGAR (chem, non-mix)', '04 - EDGAR_FF (chem, mix)',  

               '05 - EPA hr (no chem, mix)', 
               '06 - EPA hr (chem, mix)', 
               '07 - EPA hr (chem, non-mix)',    
               '08 - EPA hr_FF (chem, mix)', 

               '09 - ODIAC (no chem, mix)',  
               '10 - ODIAC (chem, mix)', 
               '11 - ODIAC (chem, non-mix)')

     # 2. scatter plot of model-data comparisons ------------------------------
     if (is.null(title))
       title = paste('b) Scatter plot over', site, 'on', timestr)

     if (is.null(mxs)) {
          colnms = c('corner', 'lons', 'lats', 'lon', 'lat')
          pm = pr %>% dplyr::select(-ends_with('_ratio')) %>% 
                      reshape2::melt(id.vars = c(colnms, 'polygon')) %>% 
                      mutate(variable = toupper(variable)) %>% 
                      filter(!variable %in% c('OBS_XCO2_FF'))
          mxs = max(pm$value, na.rm = T)
     }

     if (mxs < 1) dx = 0.1; if (mxs >= 1) dx = 0.2; if (mxs > 1.5) dx = 0.3
     if (mxs > 2) dx = 0.5; if (mxs > 5) dx = 1
     
     sz = 1.8; cols = c(ggdef.col(nrow(lm_sma)), 'black'); brks = seq(0, 10, dx)
     s1 = ggplot(data = pr) + theme_bw() + theme(legend.position = 'bottom') + 
          geom_abline(slope = 1, intercept = 0) + 
          geom_abline(data = lm_sma, aes(slope = Slope, intercept = Intercept, 
                      color = fac), linetype = 2) + 
          geom_text(data = lm_sma, 
                    aes(mxs * 0.75, mxs * 0.35 - row * (mxs / 20),
                    label = text, color = fac)) + 
          guides(size = 'none', shape = 'none') 
     
     if ('epa_tno2_mix' %in% colnames(pr)) {
          s1 = s1 + geom_errorbarh(aes(xmin = obs_tno2_tot - obs_tno2_uncert, 
                                       xmax = obs_tno2_tot + obs_tno2_uncert, 
                                       y = epa_tno2_mix), color = cols[6], 
                                   size = 0.1, height = mxs / 50, 
                                   linetype = 2, alpha = 0.8) +
                    geom_errorbar(aes(ymin = epa_tno2_mix * 0.9, 
                                      ymax = epa_tno2_mix * 1.1, 
                                      x = obs_tno2_tot), color = cols[6], 
                                      size = 0.1, alpha = 0.8, linetype = 2)
          
     } else s1 = s1 + geom_errorbarh(aes(xmin = obs_tno2_tot - obs_tno2_uncert, 
                                         xmax = obs_tno2_tot + obs_tno2_uncert, 
                                         y = edgar_tno2_mix), color = cols[2],
                                    size = 0.1, height = mxs / 50, alpha = 0.8)+
                      geom_errorbar(aes(ymin = edgar_tno2_mix * 0.9, 
                                        ymax = edgar_tno2_mix * 1.1, 
                                        x = obs_tno2_tot), color = cols[2], 
                                        size = 0.1, alpha = 0.8, linetype = 2)
     
     shp = 1
     if ('edgar_tno2_mix' %in% colnames(pr)) 
     s1 = s1 + geom_point(aes(obs_tno2_tot, edgar_tno2_nochem, 
                              color = facs[1]), shape = shp)+
          geom_point(aes(obs_tno2_tot, edgar_tno2_mix, 
                         color = facs[2]), shape = shp) 

     s1 = s1 + scale_size_manual(name = NULL, values = c(sz - 0.3, sz + 0.5)) + 
          scale_shape_manual(name = NULL, values = c(1, 20)) + 
          scale_color_manual(name = NULL, values = cols) + 
          scale_x_continuous(limits = c(mns, mxs), 
                             breaks = brks, labels = brks) +
          scale_y_continuous(limits = c(mns, mxs), 
                             breaks = brks, labels = brks) +
          labs(x = 'TROPOMI-observed tNO2 [ppb]', y = 'XSTILT-NOx tNO2 [ppb]', 
               title = title) 

     if ('epa_tno2_mix' %in% colnames(pr)) 
          s1 = s1 + geom_point(aes(obs_tno2_tot, epa_tno2_mix, 
                                   color = facs[6]), shape = shp) 

     return(s1)
}
