
# -----------------------------------------------------------------------------
plot_tau_nox = function(df, seq_df, title, xints, storeTF = F, fn = NULL, 
                        min_tau = 10, max_tau = 50, levelTF = T, hrTF = F, 
                        revTF = F) {
    
    library(viridis)
    
    if (revTF) {    # plot net loss frequency instead of net loss lifetime
        r1 = ggplot(data = df) + theme_bw() + ylim(c(min_tau, max_tau)) +
            theme(plot.title = element_text(face = 1), 
                  #text = element_text(family = 'Arial'), 
                  legend.position = 'bottom') + 
            labs(x = bquote(NO[x]~'[ppb]'), 
                 y = bquote(freq[chem]~'[hr-1]'), title = title) +
            geom_vline(xintercept = xints, color = 'gray70', linetype = 2) +
            annotate('rect', xmin = xints[1], xmax = xints[2], 
                             ymin = -Inf, ymax = Inf, 
                             fill = 'darkgreen', alpha = 0.1) + 
            annotate('rect', xmin = xints[2], xmax = xints[3], 
                             ymin = -Inf, ymax = Inf, 
                             fill = 'orange', alpha = 0.1) + 
            annotate('rect', xmin = xints[3], xmax = xints[4], 
                             ymin = -Inf, ymax = Inf, 
                             fill = 'brown', alpha = 0.1) + 
            geom_text(aes(x = indx, y = -2, label = n), size = 2, angle = 30) + 
            geom_boxplot(aes(x = indx, y = 1/tau, group = indx, fill = bin_nox),
                        size = 0.3, width = 0.7, outlier.shape = NA) + 
            scale_x_continuous(breaks = seq_df$indx[seq(1, nrow(seq_df), 2)], 
                     labels = signif(seq_df$nox[seq(1, nrow(seq_df), 2)], 2)) +
            scale_fill_viridis(option = 'inferno', name = NULL) + 
            guides(fill = 'none') +
            stat_summary(fun = median, aes(x = indx, y = 1/tau, group = indx), 
                    geom = 'point', shape = '*', color = 'white', size = 4) +
            stat_summary(fun = mean, aes(x = indx, y = 1/tau, group = indx), 
                    geom = 'point', shape = 2, color = 'white', size = 0.5) 

    } else r1 = ggplot(data = df) + theme_bw() + ylim(c(min_tau, max_tau)) +
         theme(plot.title = element_text(face = 1), 
               #text = element_text(family = 'Arial'), 
               legend.position = 'bottom')+
         labs(x = bquote(NO[x]~'[ppb]'), 
              y = bquote(tau[chem]~'[hr]'), title = title) +
         geom_vline(xintercept = xints, color = 'gray70', linetype = 2) +
         annotate('rect', xmin = xints[1], xmax = xints[2], 
                          ymin = -Inf, ymax = Inf, 
                          fill = 'darkgreen', alpha = 0.1) + 
         annotate('rect', xmin = xints[2], xmax = xints[3], 
                          ymin = -Inf, ymax = Inf, 
                          fill = 'orange', alpha = 0.1) + 
         annotate('rect', xmin = xints[3], xmax = xints[4], 
                          ymin = -Inf, ymax = Inf, 
                          fill = 'brown', alpha = 0.1) + 
         geom_text(aes(x = indx, y = -2, label = n), size = 2, angle = 30) + 
         geom_boxplot(aes(x = indx, y = tau, group = indx, fill = bin_nox), 
                      size = 0.3, width = 0.7, outlier.shape = NA) + 
         scale_x_continuous(breaks = seq_df$indx[seq(1, nrow(seq_df), 2)], 
                     labels = signif(seq_df$nox[seq(1, nrow(seq_df), 2)], 2)) +
         scale_fill_viridis(option = 'inferno', name = NULL) + 
         guides(fill = 'none') +
         stat_summary(fun = median, aes(x = indx, y = tau, group = indx), 
                      geom = 'point', shape = '*', color = 'white', size = 4) +
         stat_summary(fun = mean, aes(x = indx, y = tau, group = indx), 
                      geom = 'point', shape = 2, color = 'white', size = 0.5) 

    height = 6; width = 10
    if (levelTF) { r1 = r1 + facet_wrap(~level, ncol = 1); height = 12; width = 6 }
    if (hrTF) { r1 = r1 + facet_wrap(~date, ncol = 2); height = 7; width = 9 }
    if (storeTF) ggsave(r1, filename = fn, width = width, height = height)     
    return(r1)
}

