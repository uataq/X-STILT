# subroutine to plot OCO-2 info for each selected overpass, DW, 01/25/2019

ggmap.obs.info <- function(plotTF, site, store.path = NULL, all.timestr = NULL, 
                           oco.sensor = NULL, oco.ver = NULL, oco.path = NULL, 
                           sif.path = NULL, lon.lat = NULL, dlat.urban = NULL, 
                           dlon.urban = NULL, zoom = 8, qfTF = T) {

    if (plotTF) {
        
        plotdir <- file.path(store.path, 'plot')
        dir.create(plotdir, showWarnings = F, recursive = T)

        # zoom scales 1 to 10, larger the value, more zoomed in
        # qfTF = T for only plotting data with QF = 0
        for (t in 1 : length(all.timestr)) {

            x1 <- ggmap.obs.xco2(site, all.timestr[t], oco.sensor, oco.ver, 
                                 oco.path, lon.lat, plotdir, zoom, qfTF = qfTF, 
                                 box.dlat = dlat.urban, box.dlon = dlon.urban)

            s1 <- ggmap.obs.sif(site, all.timestr[t], oco.sensor, oco.ver, 
                                sif.path, lon.lat, plotdir, zoom)
        } # end for t

    } else {
        cat(paste('ggmap.obs.info(): NO need to plot', oco.sensor, 
                  'data as plotTF = F..\n'))
    } # end if plotTF

}