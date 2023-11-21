
plot_nox_curve = function(bin_df, mxz = 72) {

    brk.z = c(-mxz, seq(-60, 60, 10), mxz)
    lab.z = c('< -3 days', seq(-60, 60, 10), '> 3 days')
    brk.x = c(1E-3, 0.01, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000)

    pltx = bin_df$nox %>% mutate(ts_nox = ifelse(ts_nox >= mxz, mxz, ts_nox), 
                                 ts_nox = ifelse(ts_nox <= -mxz, -mxz, ts_nox))
    plt2 = bin_df$no2 %>% mutate(ts_no2 = ifelse(ts_no2 >= mxz, mxz, ts_no2), 
                                 ts_no2 = ifelse(ts_no2 <= -mxz, -mxz, ts_no2)) 
    pz = max(c(abs(pltx$ts_nox), abs(plt2$ts_no2)))

    ### lifetime contour ------------------------------------------------------
    library(RColorBrewer); cols = brewer.pal(11, 'Spectral')#'PuOr')
    s0 = ggplot() + theme_bw() + labs(y = 'SZA') + 
         theme(legend.position = 'bottom', 
               panel.grid.minor = element_blank(), 
               legend.key.width = unit(2,'cm'),
               legend.key.height = unit(0.5,'cm'))+ 
         scale_fill_gradientn(name = 'Net loss timescale\n[hr]', colors = cols, 
                              breaks = brk.z, limits = c(-pz, pz), 
                              labels = lab.z) +
         scale_x_continuous(breaks = brk.x, labels = brk.x, trans = 'log') + 
         scale_y_continuous(breaks = seq(0, 180, 10), labels = seq(0, 180, 10)) 

    s1 = s0 + labs(x = 'NOx [ppb]') + 
         geom_raster(data = pltx, aes(bin_nox, bin_sza, fill = ts_nox)) + 
         geom_contour(data = pltx, aes(bin_nox, bin_sza, z = ts_nox), 
                      linetype = 2, bins = 6, size = 0.3, color = 'gray20') +
         ggtitle(expression('[NO'[x]*'] / (L[NO'[x]*'] - P[NO'[x]*'])')) + 
         geom_hline(yintercept = 80, linetype = 2, color = 'white') 

    s2 = s0 + labs(x = 'NO2 [ppb]') + 
         geom_raster(data = plt2, aes(bin_no2, bin_sza, fill = ts_no2)) + 
         geom_contour(data = plt2, aes(bin_no2, bin_sza, z = ts_no2), 
                      linetype = 2, bins = 6, size = 0.3, color = 'gray20') +
         ggtitle(expression('[NO'[2]*'] / (L[NO'[2]*'] - P[NO'[2]*'])')) + 
         geom_hline(yintercept = 80, linetype = 2, color = 'white') 

    #s3 = s0 + geom_raster(data = pltx, aes(bin_nox, bin_sza, fill = n / 1E4)) +
    #     theme(legend.key.width = unit(0.8, 'cm'), legend.position = 'top') + 
    #     geom_hline(yintercept = 80, linetype = 2, color = 'white') +
    #     scale_fill_gradientn(name = 'Grid # (1e4)', 
    #                          colors = brewer.pal(9, 'YlGnBu'))


    # plot tau and NO2-NOx fractions against SZA -------------------------------
    title_string = '+ SZA from monthly WRF-chem runs'
    title1 = bquote(NO[x]~'net loss timescale [hr] vs.'~NO[x]~.(title_string))
    max_nox = NULL; if (site == 'Shanghai' | site == 'ALL') max_nox = 500
    s4 = plot_tau_sza(df = bin_df$nox %>% rename(tau = ts_nox), title = title1, 
                      min_tau = 0.1, max_tau = mxz, min_x = 0.003, 
                      max_x = max_nox, logTF = T, revTF = F, x = 'nox') 

    title2 = bquote(NO[2]~'net loss timescale [hr] vs.'~NO[2]~'+ SZA')
    s5 = plot_tau_sza(df = bin_df$no2 %>% rename(tau = ts_no2), title = title2, 
                      min_tau = 0.1, max_tau = mxz, min_x = 0.003, 
                      max_x = max_nox, logTF = T, revTF = F, x = 'no2') 

    s12 = ggarrange(s1, s2, ncol = 2, common.legend = T)
    #s123 = ggarrange(s12, s3, ncol = 2, widths = c(2, 1))
    s45 = ggarrange(s4, s5, ncol = 2, common.legend = T)

    ss = ggarrange(s12, s45, nrow = 2)
    ggsave(ss, filename = paste0('tau_raster_', site, '_v2022.png'), 
                width = 10, height = 10)

}