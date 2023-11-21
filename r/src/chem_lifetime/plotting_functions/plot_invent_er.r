
plot_invent_er = function(site, lon_lat, timestr, eno_fn, eco2_fn, 
                          min_eno = 0.05, epa_fn = NA, plotTF = FALSE) {

    # load emission for specific domain, unit now in umol m-2 s-1
    xmn = lon_lat$minlon
    xmx = lon_lat$maxlon 
    ymn = lon_lat$minlat 
    ymx = lon_lat$maxlat 

    eno_rt = load_eno(timestr, invent = 'edgar', emiss_fn = eno_fn, 
                      xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx)$emiss_rt  

    eco2_rt = load_eco2(timestr, invent = 'edgar', eco2_edgar_fn = eco2_fn, 
                        xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx)$emiss_rt
    
    emiss_df = as.data.frame(stack(eno_rt, eco2_rt), xy = T) %>% 
               rename(lon = x, lat = y, ENO = layer.1, ECO2 = layer.2) %>%
               mutate(ER = ENO / ECO2 * 1E3)
    fr = mean(emiss_df[emiss_df$ENO > min_eno, 'ER'])
       

    if (!is.na(epa_fn)) {

        hr_df = read.csv(file = epa_fn, row.names = NULL) 
        colnames(hr_df) = colnames(hr_df)[-1]; hr_df = hr_df[, -ncol(hr_df)]

        # convert hourly total pounds to mean flux on 0.1deg in umol m-2 s-1 
        epa_df = hr_df %>% 
                 mutate(NOx = NOx..pounds. * 0.453592 * 1E3 / 46 * 1E6 / 3600 / (0.1 * 111)^2 / 1E6) %>% 
                 filter(as.Date(Date) == 
                        as.Date(substr(timestr, 1, 8), format = '%Y%m%d')) %>% 
                 mutate(date = as.POSIXct(paste0(as.character(Date), ' ', 
                                                 formatC(Hour, width = 2, 
                                                         flag = 0)), 
                                          format = '%Y-%m-%d %H'), 
                        sf_pp = NOx / emiss_pp) %>% 
                 dplyr::select(date, sf_pp, epa_pp = NOx)
    }


    if (plotTF) {
        
        # plot emissions from inventories
        library(viridis)
        e0 = ggplot(data = emiss_df, 
                    aes(x = lon, y = lat, alpha = ENO > min_eno)) +
             scale_fill_viridis() + theme_bw() + 
             theme(legend.position = 'bottom') +
             scale_alpha_manual(values = c(0.2, 1)) + guides(alpha = 'none') 
        e1 = e0 + geom_raster(aes(fill = ENO)) + ggtitle('ENO [umol m-2 s-1]')
        e2 = e0 + geom_raster(aes(fill = ECO2)) + ggtitle('ECO2 [umol m-2 s-1]')
        e3 = e0 + geom_raster(aes(fill = ER)) + ggtitle('ER [mmol mol-1]')

        ee = ggarrange(e1, e2, e3, ncol = 3)
        ee = annotate_figure(ee, top = paste('Mean EDGAR-based ER for ENO >',
                                             min_eno,' =', signif(fr, 3), 
                                             'mmol mol-1 for', site))
        ggsave(ee, filename = paste0('EDGAR_ER_', site, '.png'), 
                   width = 10, height = 5)
    }
    
    return(emiss_df)
}