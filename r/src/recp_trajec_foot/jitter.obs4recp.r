
jitter.obs4recp = function(obs.dat, num.jitter) {

    
    # (OPTIONAL) create more receptors within each satellite sounding, 
    # useful when dealing with relatively large TROPOMI sounding, DW, 07/02/2021
    library(sp)
    tmp.dat = obs.dat %>% dplyr::select(lons, lats, polygon)
              
    # create list for each polygon, remove polygon column
    tmp.list = split(tmp.dat, tmp.dat$polygon)
    tmp.list = lapply(tmp.list, function(x) { x['polygon'] <- NULL; x })

    # concert to Polygon and add id variable for each list element
    tmp.pl = lapply(tmp.list, Polygon)    # as Polygon, but still as list
    tmp.pi = lapply(seq_along(tmp.pl), function(i) Polygons(list(tmp.pl[[i]]), 
                                                    ID = names(tmp.list)[i] ))

    # create SpatialPolygons object, class of SpatialPolygons now
    tmp.sp = SpatialPolygons(tmp.pi, proj4string = CRS("+proj=longlat +datum=WGS84") ) 
    
    # ---------------------------------
    # now sample each polygon and create new points within each polygon
    # ref from https://rstudio-pubs-static.s3.amazonaws.com/200263_e77d00d6b6b24aa8890c8c4f074bcdff.html
    # possible jitter types: hexagonal, random, regular, stratified etc. 
    tmp.spl = slot(tmp.sp, 'polygons')

    # sanity check
    if (F) {
        regular.list = sapply(tmp.spl, function(i) spsample(i, n = num.jitter, type = 'regular'))
        random.list = sapply(tmp.spl, function(i) spsample(i, n = num.jitter, type = 'random'))
        hexa.list = sapply(tmp.spl, function(i) spsample(i, n = num.jitter, type = 'hexagonal'))
        stra.list = sapply(tmp.spl, function(i) spsample(i, n = num.jitter, type = 'stratified'))
    
        plot(tmp.sp, main = "Sampling points within each sounding")
        points(random.list[[35]], col = 'red', pch = 3, cex = .5)
        points(regular.list[[33]], col = 'blue', pch = 3, cex = .5)
        points(hexa.list[[31]], col = 'darkgreen', pch = 3, cex = .5)
        points(stra.list[[29]], col = 'orange', pch = 3, cex = .5)
    }

    jitter.list = sapply(tmp.spl, function(i) spsample(i, n = num.jitter, type = 'regular'))
    jitter.df = do.call(rbind, lapply(as.list(1: length(jitter.list)), function(x) {
                        as.data.frame(unlist(jitter.list[[x]])) %>% mutate(polygon = x)
                })) %>% 
                rename(lon = x1, lat = x2) %>% 
                left_join(obs.dat %>% dplyr::select(datestr, polygon) %>% unique(), by = 'polygon') %>% 
                dplyr::select(-polygon)
    
    return(jitter.df)
}