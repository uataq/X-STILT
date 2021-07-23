#' subroutine to compute M2H background based on Hakkarainen et al., 2016
#' @author: Dien Wu, 09/16/2018

#' require: 
#' continent info, oco2 path

calc.bg.M2H <- function(lon.lat, all.timestr, output.path, oco2.ver, oco2.path, 
                        txtfile, plotTF = F) {

    bg <- NULL
    for (tt in 1:length(all.timestr)) {

        timestr <- all.timestr[tt]
        
        # minlon, maxlon, minlat, maxlat, used in Hakkareinen et al., 2016
        if (lon.lat$regid == 'Eurasia'){
            reg.lon.lat <- data.frame(minlon = -15, maxlon = 60, 
                                      minlat = 0, maxlat = 60)
        } else if (lon.lat$countryid == 'China') {
            reg.lon.lat <- data.frame(minlon = 60, maxlon = 150, 
                                      minlat = 0, maxlat = 60)
        }
        if (lon.lat$regid == 'North America') {
            reg.lon.lat <- data.frame(minlon = -130, maxlon = -60, 
                                      minlat = 0, maxlat = 60)
        }

        print(reg.lon.lat)

        if (is.null(reg.lon.lat)) {
            cat('NO region defined for backgorund\n'); return()}

        # load map
        if (plotTF) mm <- ggplot.map(map = 'black',
            minlat = reg.lon.lat$minlat, minlon = reg.lon.lat$minlon, 
            maxlat = reg.lon.lat$maxlat, maxlon = reg.lon.lat$maxlon)

        # grab observations, default is to filter QF = 0
        obs <- grab.oco2(oco2.path, timestr, reg.lon.lat, oco2.ver) %>% 
               filter(qf == 0)

        if (plotTF) {
            title <- paste('Observed XCO2 (QF = 0) for overpass on', timestr, 
                           'with regional daily median of', 
                           signif(median(obs$xco2), 5))

            m1 <- mm + labs(x = 'LONGITUDE', y = 'LATITUDE', title = title) +
                       geom_point(data = obs, aes(lon, lat, colour = xco2), size = 0.4) +
                       scale_colour_gradientn(colours = def.col(), name = 'XCO2')

            ggsave(m1, width = 9, height = 7, filename = file.path(output, 
                    paste0('M2H_bg_overpass_', oco2.ver, '_', timestr, '.png')))
        } # end if plotTF

        bg <- c(bg, median(obs$xco2))
    } # end for tt
  
    bg.df <- data.frame(timestr = all.timestr, hakka.bg = bg)
    write.table(bg.df, quote = F, row.names = F, sep = ',', file = txtfile)

    return(bg.df)
}