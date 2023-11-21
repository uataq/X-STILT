
load_eco2_vulcan = function(vulc_fn, xmin = NA, xmax = NA, ymin = NA, ymax = NA) {
                                    
    library(ncdf4); library(raster)
    vulc_dat = nc_open(vulc_fn)
    lon = ncvar_get(vulc_dat, 'lon')
    lat = ncvar_get(vulc_dat, 'lat')
    x_m = ncvar_get(vulc_dat, 'x')         # x coordinate of projection in m
    y_m = ncvar_get(vulc_dat, 'y')         # y coordinate of projection in m
    #day = ncvar_get(vulc_dat, 'time')     # days since 2010-01-01 00:00:00 UTC
    #yrs = as.numeric(substr(as.Date(day, origin = '2010-01-01 00:00:00'), 1, 4))

    # convert metric tons of carbon km-2 year-1 to umol-CO2 m-2 s-1
    eco2 = ncvar_get(vulc_dat, 'carbon_emissions') * 1e6 / 12 * 1e6 / 1e6 / (365 * 24 * 3600)
    eco2 = eco2[,,6]                        # only use the latest year of 2015
    dimnames(eco2) = dimnames(lon) = dimnames(lat) = list(x_m, y_m)
    
    # ------------------------- 
    var_name = list(c('lon', 'lat', 'eco2'))
    var_list = list(lon, lat, eco2)
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions

    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] = c('x_m', 'y_m')
    nc_close(vulc_dat)
    cat('load_eco2_vulcan(): pass ncdf loading...\n')

    # ------------------------- 
    # dealing with Vulcan's LCC projection 
    crs_lcc = '+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
    crs0 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

    # ------------------------- 
    if ( !is.na(xmin) & !is.na(xmax) & !is.na(ymin) & !is.na(ymax) ) {

        # if (xmin > -180) {
        #     library(sf)
        #     cr_df = data.frame(long = c(xmin, xmax, xmax, xmin), 
        #                        lati = c(ymin, ymin, ymax, ymax))
        #     # convert xmin, xmax, ymin, ymax to LCC projection
        #     cr_sf = st_as_sf(cr_df, coords = c('long', 'lati'), crs = crs0)
        #     cr_pj = st_transform(cr_sf, crs = crs_lcc)
        #     cr_xy = as.data.frame(st_coordinates(cr_pj))
        #     xmin = min(cr_xy$X); xmax = max(cr_xy$X)
        #     ymin = min(cr_xy$Y); ymax = max(cr_xy$Y)
        # }
        # emis_df = var_df %>% filter(x_m >= xmin, x_m <= xmax,  
        #                             y_m >= ymin, y_m <= ymax)

        emis_df = var_df %>% filter(lon >= xmin, lon <= xmax,  
                                    lat >= ymin, lat <= ymax)
    } else emis_df = var_df
    cat('load_eco2_vulcan(): pass cropping...\n')
    
    emis_rt = rasterFromXYZ(emis_df[, c('x_m', 'y_m', 'eco2')])
    crs(emis_rt) = crs_lcc
    cat('load_eco2_vulcan(): pass rasterizing...\n')

    return(emis_rt)
}
