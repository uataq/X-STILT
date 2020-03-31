# subroutine to plot OCO-2 info for each selected overpass, DW, 01/25/2019

ggmap.obs.info <- function(plotTF, site, store.path = NULL, all.timestr = NULL, 
                           oco2.ver = NULL, oco2.path = NULL, lon.lat = NULL, 
                           workdir = NULL, dlat.urban = NULL, dlon.urban = NULL) {

    if (plotTF) {

        plotdir <- file.path(store.path, 'plot')
        dir.create(plotdir, showWarnings = F, recursive = T)

        # zoom scales 1 to 10, larger the value, more zoomed in
        # qfTF = T for only plotting data with QF = 0
        for (t in 1 : length(all.timestr)) {

            x1 <- ggmap.obs.xco2(site, all.timestr[t], oco2.ver, oco2.path, 
                                 lon.lat, workdir, plotdir, zoom = 8, qfTF = T, 
                                 box.dlat = dlat.urban, box.dlon = dlon.urban)

            s1 <- ggmap.obs.sif(site, all.timestr[t], sif.path, lon.lat, workdir, 
                                plotdir, zoom = 8, box.dlon = dlon.urban, 
                                box.dlat = dlat.urban)
        } # end for t

    } else {
        cat('ggmap.obs.info(): NO need to plot OCO-2 data as plotTF = F..\n')
    } # end if plotTF

}