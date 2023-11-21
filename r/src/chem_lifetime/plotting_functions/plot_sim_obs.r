
plot_sim_obs = function(pr, lm_sma, site, timestr, emiss_df, zoom = 10, 
                        picname = NULL, mns = 0.01, mxs = NULL, 
                        labTF = T, ncol = 4, met = NULL) {
     
     m1 = map_tno2(site, pr, emiss_df, zoom, labTF, ncol, mns, mxs)
     s1 = scatter_tno2(site, timestr, pr, lm_sma)

     # merging
     setwd('/home/dienwu/postdoc_proj/NOx')
     if (is.null(picname)) {
          picname = paste0('sim_obs_no2_', site, '_', timestr, '_tm5.png')
          if (!is.null(met)) 
               picname = paste0('sim_obs_no2_', site, '_', timestr, 
                                '_', met, '_tm5.png')
     }
     print(picname)


     if ( !'epa_tno2_mix' %in% colnames(pr) ) {
          s1 = s1 + guides(color = guide_legend(ncol = 2))
          nn = ggarrange(m1, s1, ncol = 2, widths = c(ncol * 0.5, 1))
          ggsave(nn, filename = picname, width = 11, height = 5, bg = 'white')

     } else {
          
          s1 = s1 + guides(color = guide_legend(ncol = 2))
          nn = ggarrange(m1, s1, ncol = 2, widths = c(1, 0.9))
          ggsave(nn, filename = picname, width = 12, height = 7, bg = 'white')
     }

     return(m1)
}
