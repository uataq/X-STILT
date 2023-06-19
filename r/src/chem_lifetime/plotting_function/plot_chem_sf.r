plot_chem_sf = function(pr, fr, emiss_df, site, timestr, ppTF, zoom) {
     
     # 3. NO2-NOx ratio ppb/ppm -----------------------------------------------
     cols2 = RColorBrewer::brewer.pal(9, 'RdPu')
     pr = pr %>% mutate(edgar_ratio = edgar_tno2_mix / edgar_xco2)

     m0 = ggplot.map(map = 'ggmap', zoom = zoom, center.lat = mean(pr$lat),
                         center.lon = mean(pr$lon))[[1]] + coord_equal(1.25) +
          labs(x = NULL, y = NULL)

     r0 = m0 + theme(legend.position = 'bottom',
                         legend.key.height = unit(0.3, 'cm'), 
                         legend.key.width = unit(0.5, 'cm'), 
                         plot.title = element_text(size = 9)) +
          scale_fill_gradientn(colors = cols2, #limits = c(0, fr * 3), 
                              name = NULL) +
          scale_color_manual(values = c('white', 'gray50')) + 
          guides(color = 'none') 

          
     r1 = r0 + ggtitle('Modeled EnhR, tNO2 ppb / XCO2 ppm (EDGAR, chem)') + 
          geom_polygon(data = pr, aes(lons, lats, group = polygon, 
                                      color = enhTF, fill = edgar_ratio), 
                       alpha = 0.7, size = 0.2) +
          geom_text(data = pr, aes(lon, lat, label = signif(edgar_ratio, 2)), 
                                   size = 2) +
          geom_point(data = emiss_df, aes(lon, lat), shape = 4, size = 3) 
     

     r4 = r0 + ggtitle('Chemical scaling factor [ER_EDGAR / EnhR_STILT]') + 
               geom_polygon(data = pr, aes(lons, lats, group = polygon, 
                                        color = enhTF, fill = fr / edgar_ratio),
                         alpha = 0.7, size = 0.2) +
               geom_text(data = pr, aes(lon, lat, 
                         label = signif(fr/edgar_ratio, 2)), size = 2) + 
               scale_fill_gradientn(colors = cols2) + 
               geom_point(data = emiss_df, aes(lon, lat), shape = 4, size = 3) 


    if (F) {
        r2 = r0 + ggtitle('Modeled EnhR, tNO2 ppb / XCO2 ppm (EPA, chem)') + 
             geom_polygon(data = pr, aes(lons, lats, group = polygon, 
                                        color = enhTF, fill = epa_ratio), 
                          alpha = 0.7, size = 0.2) +
             geom_text(data = pr, aes(lon, lat, label = signif(epa_ratio, 2)),
                       size = 2) +
             geom_point(data = emiss_df, aes(lon, lat), shape = 4, size = 3) 

        r3 = r0 + ggtitle('Chemical scaling factor [ER_EPA / EnhR_STILT]') + 
             geom_polygon(data = pr, aes(lons, lats, group = polygon, 
                                        color = enhTF, fill = fr / epa_ratio), 
                         alpha = 0.7, size = 0.2) +
             geom_text(data = pr, aes(lon, lat, label = signif(fr/epa_ratio,2)),
                    size = 2) + scale_fill_gradientn(colors = cols2) + 
             geom_point(data = emiss_df, aes(lon, lat), shape = 4, size = 3) 

        rr = ggarrange(r1, r2, r3, r4, ncol = 2, nrow = 2)
        height = 8
    } else { rr = ggarrange(r1, r4, ncol = 2); height = 5 }
  

    ggsave(rr, filename = paste0('chemsf_', site, '_', timestr, '_tm5.png'),
               width = 8, height = height)

}