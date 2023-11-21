
# -----------------------------------------------------------------------------
plot_tau_tsfc = function(df, title, max_tau = 60, storeTF = F, tf = NULL) {

    # add solar zenith angle
    bin_tsfc = seq(-30, 50, 1)
    temp_df = df %>% mutate(bin_tsfc = bin_tsfc[findInterval(tsfc, bin_tsfc)]) %>% 
              group_by(bin_nox, bin_tsfc) %>% 
              dplyr::summarise(n = mean(n), 
                               mean_tau = mean(tau), 
                               sd_tau = sd(tau), 
                               median_tau = median(tau), 
                               meanlog_tau = exp(mean(log(tau))) ) %>% 
              ungroup() %>% filter(n > 1)
    if (storeTF) write.table(temp_df, file = tf, sep = ',', row.names = F)

    # plot tau vs. [NOx] for each temp bin
    brk.y = c(0.5, 1, 2, 5, 10, 15, 20, 30, 40, 60, 80)
    brk.x = c(0.2, 0.5, 1, 2, 5, 10, 20, 30, 40, 60, 80, 100, 120)

    tsfc_cols = rev(RColorBrewer::brewer.pal(11, 'Spectral'))
    t1 = ggplot(data = temp_df, aes(x = bin_nox)) + theme_bw() + 
         theme(legend.position = 'bottom', panel.grid.minor = element_blank(), 
               legend.key.width = unit(2, 'cm'), 
               legend.key.height = unit(0.5, 'cm')) + 
         geom_hline(yintercept = c(0.5, 1, 5), linetype = 2, color = 'gray50') +
         geom_point(aes(y = meanlog_tau, color = bin_tsfc), 
                    shape = 42, size = 7) + 
         scale_linetype_manual(name = NULL, values = c(2, 2)) + 
         scale_color_gradientn(name = 'Binned Tsfc', colors = tsfc_cols, 
                               breaks = bin_tsfc, labels = bin_tsfc) + 
         scale_y_continuous(breaks = brk.y, labels = brk.y, 
                            limits = c(0, max_tau)) +
         scale_x_continuous(breaks = brk.x, labels = brk.x, trans = 'log') + 
         labs(x = bquote(NO[x]~'[ppb]'), y = bquote(tau[chem]~'[hr]'), 
              title = title)

    return(t1)
}
 
