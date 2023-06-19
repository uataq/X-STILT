
# -----------------------------------------------------------------------------
plot_tau_sza = function(df, title, min_ts = 0, max_ts = 50, min_x = 0.003, 
                        max_x = NULL, max_sza = 80, logTF = T, revTF = F, 
                        x = c('nox', 'no2', 'vocr'), cols = NULL) {
    extrafont::loadfonts(device = "postscript")

    # plot ts vs. [NOx] for each SZA bin
    brk.y = sort(unique(c(0.1, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, max_ts)))
    brk.y = brk.y[brk.y <= max_ts]
    if (revTF) brk.y = seq(0, 0.5, 0.1)
    lab.y = brk.y; lab.y[lab.y == max_ts] = expression(infinity)
    
    if (is.null(max_x)) {
        if (x == 'nox') max_x = max(ceiling(df$bin_nox))
        if (x == 'no2') max_x = max(ceiling(df$bin_no2))
        if (x == 'vocr') max_x = max(ceilling(df$bin_vocr))
    }

    brk.x = c(1E-3, 3E-3, 0.01, 0.03, 0.1, 0.3, 1, 2, 3, 5, 10, 30, 50, 100, 
              200, 500, 1000)

    brk.z = seq(0, max_sza, 10)
    #cols = RColorBrewer::brewer.pal(11, 'Spectral')
    if (is.null(cols)) cols = RColorBrewer::brewer.pal(11, 'RdYlBu')
    df = df %>% filter(ts_nox >= min_ts, ts_nox <= max_ts)

    if (x == 'nox') { xlab = bquote(NO[x]~' bins [ppb]')
                      s1 = ggplot(data = df, aes(x = bin_nox)) 
                      df = df %>% filter(bin_nox >= min_x, bin_nox <= max_x) }
    if (x == 'no2') { xlab = bquote(NO[2]~'bins [ppb]')
                      s1 = ggplot(data = df, aes(x = bin_no2))
                      df = df %>% filter(bin_no2 >= min_x, bin_no2 <= max_x) }
    if (x == 'vocr') { xlab = bquote(VOC[R]~'bins [s-1]')
                      s1 = ggplot(data = df, aes(x = bin_vocr))
                      df = df %>% filter(bin_vocr >= min_x, bin_vocr <= max_x) }
    
    s1 = s1 + theme_bw() + theme(legend.position = 'bottom', 
                                 #text = element_text(family = 'Arial'),
                                 panel.grid.minor = element_blank(), 
                                 legend.key.width = unit(1.4, 'cm'), 
                                 legend.key.height = unit(0.4, 'cm')) + 
              scale_x_continuous(breaks = brk.x, labels = brk.x, trans = 'log', 
                                 limits = c(min_x, max_x)) + 
              scale_color_gradientn(name = 'Binned SZA', colors = cols, 
                                    limits = c(0, max_sza),
                                    breaks = brk.z, labels = brk.z) 
    
    if (logTF) { 
        s1 = s1 + scale_y_continuous(breaks = brk.y, labels = lab.y, 
                                     trans = 'log10', 
                                     limits = c(min_ts, max_ts))
    } else s1 = s1 + scale_y_continuous(breaks = brk.y, labels = lab.y, 
                                        limits = c(min_ts, max_ts))
    
    if ('ts_nox_sd' %in% colnames(df) & !logTF) 
        s1 = s1 + geom_ribbon(aes(ymin = ts_nox - ts_nox_sd, 
                                  ymax = ts_nox + ts_nox_sd, 
                                  fill = bin_sza, group = bin_sza), 
                                  alpha = 0.01, size = 0.2) + 
                  scale_fill_gradientn(colors = cols, limits = c(0, max_sza),
                                       breaks = brk.z, labels = brk.z) + 
                  guides(fill = 'none')

    if (revTF) {    # plot frequencies
        s1 = s1 +geom_hline(yintercept = 0.1, linetype = 2, color = 'gray50')+
             geom_point(aes(y = 1/ts_nox, color = bin_sza), 
                        shape = 42, size = 7) +
             scale_linetype_manual(name = NULL, values = c(2, 2)) + 
             labs(x = xlab, y = bquote(freq[chem]~'[hr-1]'), title = title)
        
    } else s1 = s1 + 
                geom_hline(yintercept = 1, linetype = 2, color = 'gray50') +
                geom_point(aes(y = ts_nox, color = bin_sza), 
                           shape = 19, size = 0.7) +
                geom_smooth(aes(y = ts_nox, color = bin_sza, group = bin_sza),
                            se = F, linewidth = 0.1, alpha = 0.3, span = 0.6) + 
                scale_linetype_manual(name = NULL, values = c(2, 2)) + 
                labs(x = xlab, y = bquote(ts[chem]~'[hr]'), title = title)

    if ('fac' %in% colnames(df)) 
        s1 = s1 + facet_wrap(~fac, ncol = 1)

    return(s1)
}
 


# -----------------------------------------------------------------------------
# NO2-NOx fractions with SZA and temp
plot_ratio_sza = function(df, title, min_sza = 0, max_sza = 180) {

    # plot fraction vs. SZA for each SZA bin
    brk.y = seq(0.1, 1.0, 0.1)
    brk.x = seq(-40, 60, 5)
 
    sza_cols = RColorBrewer::brewer.pal(11, 'Spectral')
    f1 = ggplot(data = df, aes(x = bin_tc)) + theme_bw() + 
         theme(legend.position = 'bottom', 
               panel.grid.minor = element_blank(), 
               legend.key.width = unit(2.2, 'cm'), 
               legend.key.height = unit(0.5, 'cm')) + 
         geom_hline(yintercept = c(0.5, 0.8), linetype = 2, color = 'gray50') +
         geom_point(aes(y = ratio, color = bin_sza), shape = 21, size = 1.2) + 
         scale_linetype_manual(name = NULL, values = c(2, 2)) + 
         scale_y_continuous(breaks = brk.y, labels = brk.y) +
         scale_x_continuous(breaks = brk.x, labels = brk.x) + 
         geom_smooth(aes(y = ratio, color = bin_sza, group = bin_sza),
                         se = F, size = 0.1, alpha = 0.4) + 
         scale_color_gradientn(name = 'Binned SZA', colors = sza_cols, 
                               limits = c(min_sza, max_sza),
                               breaks = seq(0, max_sza, 10), 
                               labels = seq(0, max_sza, 10)) + 
         labs(x = 'Temp [degC]', y = 'NO2:NOx ratio', title = title)
    return(f1)
}
 