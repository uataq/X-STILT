

color_rmp = function(num = NULL) {
    cols1 = RColorBrewer::brewer.pal(3, 'Dark2')
    cols2 = RColorBrewer::brewer.pal(3, 'Accent')
    cols3 = RColorBrewer::brewer.pal(3, 'Paired')

    cols2 = c(cols2[1], cols2[3], cols2[2])
    cols = c(cols1[1], cols2[1], 
             cols1[2], cols2[2], 
             cols1[3], cols2[3], 
             cols3[2], cols3[1])
    cols = c(cols, 'black')
    if (!is.null(num)) cols = cols[1:num]
    return(cols)
}

plot_rmp = function(plt_df) {

    brks = seq(-180, 180, 60)
    cols = color_rmp()
    f0 = ggplot(data = plt_df, aes(x = angle)) + theme_bw() + 
         theme(legend.position = 'bottom') +
         guides(color = guide_legend(nrow = 2)) +
         labs(x = 'Rotated angles [deg]', 
              subtitle = paste('as a function of rotated angles [degrees] for', sel_timestr)) + 
         scale_x_continuous(breaks = brks, labels = brks, 
                            sec.axis = sec_axis(~.*pi / 180, name = 'Radian',
                                            breaks = brks * pi / 180, 
                                            labels = signif(brks * pi / 180, 3))) + 
         scale_color_manual(name = 'UNrotated X rotated', values = cols) +
         scale_fill_manual(name = 'UNrotated X rotated', values = cols) +
         geom_vline(data = plt_df,# %>% filter(!variable %in% nochem_vars), 
                    aes(xintercept = mu_deg, color = fac), 
                    linetype = 2, alpha = 0.8)

    library(ggrepel)
    txt_df = unique(plt_df[, c('pred_rmp_mx', 'pred_rmp_mn', 
                               'auc', 'fac', 'variable', 'timestr')])

    f1 = f0 + geom_point(aes(y = rmp, color = fac), size = 1.4) + 
         labs(y = 'SQRT(MEAN(product)) [ppb]', 
              title = expression('a) Root-mean-product (RMP) between two sets of tNO'[2]*' plumes [ppb]')) +
         scale_y_continuous(breaks = seq(0, 10, 0.1), 
                            labels = seq(0, 10, 0.1)) + 
         geom_line(aes(y = pred_rmp, color = fac, group = fac), 
                   linewidth = 0.4) +
         geom_ribbon(aes(ymin = pred_rmp_mn, ymax = pred_rmp, xmin = -180, 
                         xmax = 180, color = fac, group = fac, fill = fac), 
                         alpha = 0.1) + 
         geom_text(data = txt_df, aes(x = 5, y = pred_rmp_mx + 0.005, 
                   label = signif(auc / pi * 180, 3), color = fac)) +
         geom_text(data = txt_df, aes(x = -120, y = pred_rmp_mn -0.008, 
                   label = signif(pred_rmp_mn, 3), color = fac)) 

    f2 = f0 + geom_point(aes(y = norm, color = fac), size = 1.5) + 
         scale_y_continuous(breaks = seq(0, 1, 0.4), labels = seq(0, 1, 0.4)) + 
         labs(y = 'Normalized values', 
              title = expression('b) Normalized RMP between two sets of tNO'[2]*' plumes [unitless]')) + 
         geom_line(aes(y = pred_norm, color = fac, group = fac), 
                   linewidth = 0.4) 

    ff = ggarrange(f1, f2, ncol = 2, common.legend = T, legend = 'bottom')
    return(ff)
}