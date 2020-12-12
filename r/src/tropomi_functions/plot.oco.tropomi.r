

plot.oco.tropomi <- function(site, timestr, lon.lat, xco2.obs, sif.obs = NULL, 
                             xco.obs, xno2.obs, oco.sensor, xco2.qf = T, 
                             xco.qa = 0.4, xno2.qa = 0.7, zoom = 9, zsfcTF = T,
                             plot.dir = './plot') {

    library(ggpubr)
    #if (grepl('fossil_', site)) site = gsub('fossil_', '', site)
    tropomi.hr <- unique(substr(xco.obs$time_utc, 9, 10))
    oco.hr <- unique(substr(xco2.obs$timestr, 9, 10))

    # plot google map
    m1 <- ggplot.map(map = 'ggmap', zoom = zoom, center.lat = lon.lat$citylat,
                     center.lon = lon.lat$citylon)[[1]] 
    col <- def.col()[-c(1, length(def.col()))]
    #col <- rev(brewer.pal(11, 'RdYlBu'))


    ### --------------------------- plot XCO2 using vertex lat/lon 
    if (xco2.qf) xco2.obs <- xco2.obs %>% filter(qf == 0)
    c2 <- m1 + theme_bw() + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_polygon(data = xco2.obs, aes(lons, lats, fill = xco2, group = indx), 
                       alpha = 0.7, color = NA, size = 0.5) + 
          scale_fill_gradientn(name = 'XCO2 [ppm]', colours = col) +
          labs(title = paste(oco.sensor, 'XCO2 [QF = 0] on', 
                             substr(timestr, 1, 8), oco.hr, 'UTC')) + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                legend.key.width = unit(1.2, 'cm'))


    ### --------------------------- plot xCO
    df1 <- xco.obs %>% filter(qa >= xco.qa) 
    zero.indx <- which(df1$corner == 0)
    df1$group <- findInterval(as.numeric(rownames(df1)), zero.indx)

    c1 <- m1 + theme_bw() + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_polygon(data = df1, 
                       aes(lons, lats, fill = xco * 6.02214E19, group = group), 
                       alpha = 0.7, color = 'white', size = 0.5) + 
          scale_fill_gradientn(name = 'XCO\n[molec cm-2]', colours = col) +
          labs(title = paste0('TROPOMI XCO [QA >= ', xco.qa, '] on ', 
                              substr(timestr, 1, 8), ' ', tropomi.hr, ' UTC')) + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                legend.key.width = unit(1.2, 'cm'))



    ### --------------------------- plot tropo xNO2
    df3 <- xno2.obs %>% filter(qa >= xno2.qa) %>% 
           mutate(tropo_xno2 = ifelse(tropo_xno2 < 0, 0, tropo_xno2)) 
    #count.corner = table(df3$corner)
    #max.corner.indx = names(count.corner)[count.corner == max(count.corner)]
    zero.indx <- which(df3$corner == 0)
    df3$group <- findInterval(as.numeric(rownames(df3)), zero.indx)

    n1 <- m1 + theme_bw() + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_polygon(data = df3, 
                       aes(lons, lats, fill = tropo_xno2 * 6.02214E19, group = group), 
                       alpha = 0.7, color = 'white', size = 0.5) + 
          scale_fill_gradientn(name = 'Tropospheric NO2\n[molec cm-2]', colours = col) +
          labs(title = paste0('Tropospheric NO2 [QA >= ', xno2.qa, '] on ', 
                              substr(timestr, 1, 8), ' ', tropomi.hr, ' UTC')) + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                legend.key.width = unit(1.2, 'cm'))


    ### plot sfc altitude 
    width = 17
    if (zsfcTF) {
        h1 <- m1 + theme_bw() +
            geom_polygon(data = df3, aes(lons, lats, fill = hsfc, group = group), 
                        alpha = 0.7, color = 'white', size = 0.5) + 
            scale_fill_gradientn(name = 'Zsfc [m]', colours = col) +
            labs(title = 'TROPOMI-retrieved Zsfc', x = 'LONGITUDE', y = 'LATITUDE') + 
            theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                    legend.key.width = unit(1.2, 'cm'))
        ccnh <- ggarrange(c2, c1, n1, h1, ncol = 4)
        width = 13
    } else ccnh <- ggarrange(c2, c1, n1, ncol = 3)


    ### --------------------------- plot SIF and land cover
    height = 6; all = ccnh
    fn <- paste0(oco.sensor, '_XCO2_XCO_tNO2_', site, '_', substr(timestr, 1, 8), '.png')

    if (!is.null(sif.obs)) {
        if (!nrow(sif.obs) == 0) {
            melt.sif <- sif.obs %>% 
                    dplyr::select('timestr', 'lat', 'lon', 'sif757', 'sif771', 'avg.sif') %>% 
                    melt(id.var = c('timestr', 'lat', 'lon'))

            title <- paste(oco.sensor, 'SIF [W/m2/sr/µm] and IGBP on', timestr)
            s1 <- m1 + geom_point(data = melt.sif, aes(lon, lat, colour = value), size = 0.9) +
                facet_wrap(~variable, ncol = 3) + theme(legend.position = 'bottom') + 
                scale_colour_gradientn(name = paste(oco.sensor, 'SIF'), colours = col,
                                        limits = c(-1, max(2.5, max(melt.sif$value))),
                                        breaks = seq(-4, 4, 0.5), labels = seq(-4, 4, 0.5)) + 
                labs(x = NULL, y = NULL, title = title) + 
                theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                        legend.key.width = unit(1.2, 'cm'))

            l1 <- ggmap.igbp(sif.obs, m1, legend.ncol = 5) + labs(x = NULL, y = NULL) + 
                theme(legend.position = 'bottom')
            
            si  <- ggarrange(s1, l1, ncol = 2, widths = c(2.8, 1))
            all <- ggarrange(ccnh, si, nrow = 2, heights = c(1, 0.9))
            height = 10
            fn <- paste0(oco.sensor, '_XCO2_SIF_XCO_tNO2_', site, '_', substr(timestr, 1, 8), '.png')
        } 
    }
    
    cat(paste('Saving plot in', file.path(plot.dir, fn), '\n\n'))
    all = ggpubr::annotate_figure(all, top = site)
    ggsave(all, filename = file.path(plot.dir, fn), width = width, height = height)
}

