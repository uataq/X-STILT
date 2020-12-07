#' subroutine to grab TROPOMI data, given a spatial domain and a date
#' @author: Dien Wu

### ------------------------------- 
# tropomi module used in main script, config_tropomi() being called by config_xstilt()
config_tropomi <- function(timestr, tropomi.speci, tropomi.path, lon.lat) {

    # get TROPOMI info, if tropomi overpass hour differs from OCO, DW, 10/26/2020
    if (!NA %in% tropomi.speci) {
        
        oco.hr <- substr(timestr, 9, 10)
        tropomi.fn <- NULL 
        for (ss in 1 : length(tropomi.speci)) {

            cat(paste('- Obtaining TROPOMI', tropomi.speci[ss], 'info...\n'))
            tmp.path <- tropomi.path[grep(tropomi.speci[ss], tropomi.path)]
            tropomi.info <- find.tropomi(tropomi.path = tmp.path, 
                                         timestr = substr(timestr, 1, 8), 
                                         lon.lat = lon.lat)
            
            if (ss == length(tropomi.speci)) {

                # if TRUE, use TROPOMI overpass hour instead of OCO
                # while still keep the same lat/lon soundings of OCO as receptor locations
                # if FALSE, combine trajec simulation as one run, given their same overpass hour
                tropomi.hr <- substr(tropomi.info[, c(1, 2)], 10, 11)
                cat(paste('\nOCO overpass hour:', oco.hr, 
                          'UTC; TROPOMI overpass hour:', tropomi.hr[1], 
                          'to', tropomi.hr[2], 'UTC\n'))

                tropomiTF <- !as.numeric(oco.hr) %in% seq(as.numeric(tropomi.hr[1]), 
                                                          as.numeric(tropomi.hr[2]), 1)
            }   # end if ss

            tropomi.fn <- c(tropomi.fn, file.path(tmp.path, tropomi.info$fn))
        }   # end for ss
       
    } # end if NA 

    return(list(tropomiTF = tropomiTF, tropomi.fn = tropomi.fn))
}



### ------------------------------- 
# timestr needs to have a format of YYYYMMDD
find.tropomi <- function(tropomi.path, timestr, lon.lat) {

    library(ncdf4)    
    fns <- list.files(tropomi.path, paste0('____', timestr))
    if (length(fns) == 0) { cat('NO TROPOMI files found...please check..\n'); return() }
    
    tropomi_df <- NULL
    for (fn in fns) {
        
        dat <- nc_open(file.path(tropomi.path, fn))
        indx_along_track  <- ncvar_get(dat, 'PRODUCT/scanline')
        indx_across_track <- ncvar_get(dat, 'PRODUCT/ground_pixel')

        ## variables
        lat <- ncvar_get(dat, 'PRODUCT/latitude')
        lon <- ncvar_get(dat, 'PRODUCT/longitude')

        # for CO and NO2, 

        # for CH4, A continuous quality descriptor, varying between 0 (no data) 
        # and 1 (full quality data). Recommend to ignore data with qa_value < 0.5
        qa  <- ncvar_get(dat, 'PRODUCT/qa_value')  # quality assurances


        if (unique(grepl('_CO_', fn))) 
            val <- ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column')
        
        # bias corrected dry_atmosphere_mole_fraction_of_methane
        if (unique(grepl('_CH4_', fn)))     # unit of 1e-9
            val <- ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_bias_corrected')
            #val <- ncvar_get(dat, 'PRODUCT/methane_mixing_ratio')

        if (unique(grepl('_NO2_', fn))) 
            val <- ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column')

        dimnames(lat) <- dimnames(lon) <- dimnames(qa) <- dimnames(val) <- 
            list(indx_across_track, indx_along_track)

        # merge matrix 
        var_name <- list(c('lat', 'lon', 'qa', 'val'))
        var_list <- list(lat, lon, qa, val)
        
        #it is assumed all matrices in the list have equal dimensions
        var_array <- array(
            data = do.call(cbind, var_list), 
            dim = c(dim(var_list[[1]]), length(var_list)), 
            dimnames = c(dimnames(var_list[[1]]), var_name)
        )
        
        var_df <- dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
        colnames(var_df)[1:2] <- c('indx_across_track', 'indx_along_track')

        # qa_value > 0.75. For most users this is the recommended pixel filter. 
        # This removes cloud-covered scenes (cloud radiance fraction > 0.5), 
        # part of the scenes covered by snow/ice, errors and problematic retrievals.
        # qa_value > 0.50. This adds the good quality retrievals over clouds and over scenes covered by snow/ice. 
        # Errors and problematic retrievals are still filtered out. 
        # In particular this choice is useful for assimilation and model comparison
        sel_df <- var_df %>% filter(!is.na(val), 
                                    lon >= lon.lat$minlon, lon <= lon.lat$maxlon,
                                    lat >= lon.lat$minlat, lat <= lon.lat$maxlat) 
        nc_close(dat)

        if (nrow(sel_df) > 0) {
            cat(paste0('found the data that match the criteria; see nc file: ', 
                       fn, '; stop searching...\n'))
            fn_df <- strsplit.to.df(basename(fn)) 

            if (unique(grepl('_CO_', fn))) { start.time = fn_df$V10; end.time = fn_df$V11 }
            if (unique(grepl('_CH4_', fn))) { start.time = fn_df$V9; end.time = fn_df$V10 } 
            if (unique(grepl('_NO2_', fn))) { start.time = fn_df$V9; end.time = fn_df$V10 } 
            
            tmp_df <- data.frame(start.time = start.time, end.time = end.time, 
                                 tot.count = nrow(sel_df), 
                                 qa0p5.count = nrow(sel_df %>% filter(qa >= 0.50)), 
                                 qa0p4.count = nrow(sel_df %>% filter(qa >= 0.40)), 
                                 qa0p7.count = nrow(sel_df %>% filter(qa >= 0.70)), 
                                 fn = fn, stringsAsFactors = F)
            tropomi_df <- rbind(tropomi_df, tmp_df)
            return(tropomi_df)
        }  # end if
    }   # end for fn

    if (is.null(tropomi_df)) {
        cat('Cannot find any TROPOMI data that match your time and spatial domain...\n')
        return()
    }
}
# end of subrouti




### ------------------------------- 
# timestr needs to have a format of YYYYMMDD
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
grab.tropomi.no2 <- function(tropomi.path = NULL, timestr = NULL, lon.lat, 
                             tropomi.fn = NULL) {

    library(ncdf4)

    if (is.null(tropomi.fn)) {
        
        if (nchar(timestr) > 8) timestr <- substr(timestr, 1, 8)
        fn <- list.files(tropomi.path, paste0('____', timestr))

        if (length(fn) > 1) { 
            cat('grab.tropomi.no2(): Multiple TROPOMI files;\nX-STILT is looking for the tropomi file that has soundings over your spatial domain\n')
            tropomi.info <- find.tropomi(tropomi.path, timestr, lon.lat)
            fn <- tropomi.info$fn
        }
        if (length(fn) == 0) { cat('NO TROPOMI NO2 files found..\n'); return() }
        fn <- file.path(tropomi.path, fn)

    } else fn <- tropomi.fn 

    # get dimension first, index starts with ZERO
    dat <- nc_open(fn)
    indx_along_track  <- ncvar_get(dat, 'PRODUCT/scanline')
    indx_across_track <- ncvar_get(dat, 'PRODUCT/ground_pixel')
    corner <- ncvar_get(dat, 'PRODUCT/corner')

    ## variables
    lat  <- ncvar_get(dat, 'PRODUCT/latitude')  # centered lat/lon
    lon  <- ncvar_get(dat, 'PRODUCT/longitude')
    lats <- ncvar_get(dat, 'GEOLOCATIONS/latitude_bounds')
    lons <- ncvar_get(dat, 'GEOLOCATIONS/longitude_bounds')

    time <- substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
    timestr <- format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), format = '%Y%m%d%H%M%S')
    time_mtrx <- t(replicate(length(indx_across_track), as.numeric(timestr)))
    qa   <- ncvar_get(dat, 'PRODUCT/qa_value')  # quality assurances
    hsfc <- ncvar_get(dat, 'INPUT_DATA/surface_altitude')
    psfc <- ncvar_get(dat, 'INPUT_DATA/surface_pressure') / 100 # Pa to hPa
    amf_tot <- ncvar_get(dat, 'PRODUCT/air_mass_factor_total')
    amf_tropo <- ncvar_get(dat, 'PRODUCT/air_mass_factor_troposphere')

    # "TM5 layer index of the highest layer in the tropopause"
    tropo_indx <- ncvar_get(dat, 'PRODUCT/tm5_tropopause_layer_index')
    sfc_class <- ncvar_get(dat, 'INPUT_DATA/surface_classification')
    tropo_xno2 <- ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column')
    tropo_xno2_uncert <- ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision')
    tropo_xno2_kernel <- ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision_kernel')

    # assign proper dimensions to variables
    dimnames(lat) <- dimnames(lon) <- dimnames(qa) <- dimnames(time_mtrx) <- 
    dimnames(hsfc) <- dimnames(psfc) <- dimnames(sfc_class) <- 
    dimnames(tropo_indx) <- dimnames(tropo_xno2) <- dimnames(tropo_xno2_uncert) <- 
        list(indx_across_track, indx_along_track)
    
    dimnames(lats) <- dimnames(lons) <- list(corner, indx_across_track, indx_along_track)

    # -------------------------  merge matrix 1
    var_name <- list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                       'psfc', 'sfc_class', 'tropo_k', 'tropo_xno2', 'tropo_xno2_uncert'))
    var_list <- list(lat, lon, time_mtrx, qa, hsfc, psfc, sfc_class, tropo_indx,
                     tropo_xno2, tropo_xno2_uncert)    
    
    var_array <- array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    var_df <- dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] <- c('indx_across_track', 'indx_along_track')
    
    # select data frame based on spatial domain 
    sel_df <- var_df %>% filter(center_lon >= lon.lat$minlat - 0.05, 
                                center_lon <= lon.lat$maxlon + 0.05,
                                center_lat >= lon.lat$minlat - 0.05, 
                                center_lat <= lon.lat$maxlat + 0.05, !is.na(tropo_k)) 
    sel_df <- sel_df %>% mutate(sample = 1 : nrow(sel_df)) 
    
    # calculate tropopause pressure, need pressure levels and tropopause level indx
    # TM5 hybrid A coefficient at upper and lower interface levels in Pa, dim of [34, 2]
    a <- ncvar_get(dat, 'PRODUCT/tm5_constant_a') / 100   # convert to hPa now
    # TM5 hybrid B coefficient at upper and lower interface levels, unitless, dim of [34, 2]
    b <- ncvar_get(dat, 'PRODUCT/tm5_constant_b')   
    # k from surface to top of atmosphere, 0 to 33 
    k <- ncvar_get(dat, 'PRODUCT/layer')  
    # vertices: lower or upper pressure level boundaries, v = 0 for lower level
    v <- ncvar_get(dat, 'PRODUCT/vertices')  

    # grab pressure coefficient for lower and upper level boundaries
    ab_df <- data.frame(a_lower = rep(a[1, ], nrow(sel_df)), 
                        b_lower = rep(b[1, ], nrow(sel_df)), 
                        a_upper = rep(a[2, ], nrow(sel_df)), 
                        b_upper = rep(b[2, ], nrow(sel_df)), 
                        k = rep(k, nrow(sel_df)),   # from sfc to TOA, indx = 0 to 33
                        sample = rep(1 : nrow(sel_df), each = length(a[1, ])))

    add_df <- sel_df %>% left_join(ab_df, by = 'sample') %>% group_by(sample) %>%
                         filter(k == tropo_k) %>% ungroup() %>% 
                         mutate(tropo_lower_pres = a_lower + b_lower * psfc, 
                                tropo_upper_pres = a_upper + b_upper * psfc)
    #print(range(add_df$tropo_lower_pres))

    # ------------------------- merge matrix 2
    loc_name <- list(c('lats', 'lons')); loc_list <- list(lats, lons)
    loc_array <- array(
        data = do.call(cbind, loc_list), 
        dim = c(dim(loc_list[[1]]), length(loc_list)), 
        dimnames = c(dimnames(loc_list[[1]]), loc_name)
    )
    
    loc_df <- dcast(melt(loc_array), Var1 + Var2 + Var3~Var4, value.var = 'value')
    colnames(loc_df)[1:3] <- c('corner', 'indx_across_track', 'indx_along_track')
    loc_df <- loc_df %>% filter(lons >= lon.lat$minlon - 0.05, 
                                lons <= lon.lat$maxlon + 0.05,
                                lats >= lon.lat$minlat - 0.05, 
                                lats <= lon.lat$maxlat + 0.05)

    merge_df <- right_join(add_df, loc_df, by = c('indx_across_track', 'indx_along_track'))
    nc_close(dat)

    return(merge_df)
}
# end of subroutinr



### ------------------------------- 
# timestr should be format of YYYYMMDD, if it has > 8 letters, we will simply truncate it 
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
grab.tropomi.co <- function(tropomi.path = NULL, timestr = NULL, lon.lat, 
                            getakTF = F, tropomi.fn = NULL) {

    library(ncdf4)
    if (is.null(tropomi.fn)) {

        if (nchar(timestr) > 8) timestr <- substr(timestr, 1, 8)
        fn <- list.files(tropomi.path, paste0('____', timestr))

        if (length(fn) > 1) { 
            cat('grab.tropomi.co(): Multiple TROPOMI files;\nX-STILT is looking for the tropomi file that has soundings over your spatial domain\n')
            tropomi.info <- find.tropomi(tropomi.path, timestr, lon.lat)
            fn <- tropomi.info$fn
        }
        if (length(fn) == 0) { cat('NO TROPOMI CO files found..\n'); return() }
        fn <- file.path(tropomi.path, fn)

    } else fn <- tropomi.fn 

    # get dimension first, index starts with ZERO
    dat <- nc_open(fn)
    indx_along_track  <- ncvar_get(dat, 'PRODUCT/scanline')
    indx_across_track <- ncvar_get(dat, 'PRODUCT/ground_pixel')
    corner <- ncvar_get(dat, 'PRODUCT/corner')

    ## variables
    lat <- ncvar_get(dat, 'PRODUCT/latitude')
    lon <- ncvar_get(dat, 'PRODUCT/longitude')
    lats <- ncvar_get(dat, 'GEOLOCATIONS/latitude_bounds')
    lons <- ncvar_get(dat, 'GEOLOCATIONS/longitude_bounds')

    time <- substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
    timestr <- format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), format = '%Y%m%d%H%M%S')
    time_mtrx <- t(replicate(length(indx_across_track), as.numeric(timestr)))
    qa   <- ncvar_get(dat, 'PRODUCT/qa_value')  # quality assurances
    hsfc <- ncvar_get(dat, 'INPUT_DATA/surface_altitude')
    psfc <- ncvar_get(dat, 'INPUT_DATA/surface_pressure') / 100 # Pa to hPa
    sfc_class <- ncvar_get(dat, 'INPUT_DATA/surface_classification')
    xco <- ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column') # mol m-2
    xco_uncert <- ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column_precision')
    xh2o <- ncvar_get(dat, 'DETAILED_RESULTS/water_total_column')

    # assign proper dimensions to variables
    dimnames(lat) <- dimnames(lon) <- dimnames(qa) <- dimnames(time_mtrx) <- 
    dimnames(hsfc) <- dimnames(psfc) <- dimnames(sfc_class) <- dimnames(xco) <- 
    dimnames(xco_uncert) <- dimnames(xh2o) <- list(indx_across_track, indx_along_track)

    dimnames(lats) <- dimnames(lons) <- list(corner, indx_across_track, indx_along_track)

    # -------------------------  merge matrix 1
    var_name <- list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                       'psfc', 'sfc_class', 'xco', 'xco_uncert', 'xh2o'))
    var_list <- list(lat, lon, time_mtrx, qa, hsfc, psfc, sfc_class, xco, xco_uncert, xh2o)    
    var_array <- array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    var_df <- dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] <- c('indx_across_track', 'indx_along_track')
    sel_df <- var_df %>% filter(!is.na(xco), center_lon >= lon.lat$minlon, 
                                center_lon <= lon.lat$maxlon,
                                center_lat >= lon.lat$minlat, 
                                center_lat <= lon.lat$maxlat) 
    
    # ------------------------- merge matrix 2
    loc_name <- list(c('lats', 'lons')); loc_list <- list(lats, lons)
    loc_array <- array(
        data = do.call(cbind, loc_list), 
        dim = c(dim(loc_list[[1]]), length(loc_list)), 
        dimnames = c(dimnames(loc_list[[1]]), loc_name)
    )
    
    loc_df <- dcast(melt(loc_array), Var1 + Var2 + Var3~Var4, value.var = 'value')
    colnames(loc_df)[1:3] <- c('corner', 'indx_across_track', 'indx_along_track')
    loc_df <- loc_df %>% filter(!is.na(lats), lons >= lon.lat$minlon - 0.1, 
                                lons <= lon.lat$maxlon + 0.1,
                                lats >= lon.lat$minlat - 0.1, 
                                lats <= lon.lat$maxlat + 0.1) 
    
    # ----------------------- if grab CO AK and convert from matrix to df
    if (getakTF) {
        cat('reading CO column averaging kernel...\n')
        ak <- ncvar_get(dat, 'DETAILED_RESULTS/column_averaging_kernel')    # ak in meter
        layer <- ncvar_get(dat, 'PRODUCT/layer')    # height in m
        dimnames(ak) <- list(layer, indx_across_track, indx_along_track)
        
        ak.df <- melt(ak) %>% dplyr::rename(hgt = Var1, indx_across_track = Var2, 
                                            indx_along_track = Var3, ak = value) 
        
        # will calculate the surface normalized AK and 
        sel_df <- ak.df %>% na.omit() %>% 
                  left_join(var_df , by = c('indx_across_track', 'indx_along_track')) %>%
                  filter(center_lon >= lon.lat$minlon, center_lon <= lon.lat$maxlon,
                         center_lat >= lon.lat$minlat, center_lat <= lon.lat$maxlat) %>% 
                  ungroup()
    }   # end if


    merge_df <- right_join(sel_df, loc_df, by = c('indx_across_track', 'indx_along_track'))
    nc_close(dat)

    return(merge_df)
}
# end of subroutinr



plot.oco.tropomi <- function(site, timestr, lon.lat, xco2.obs, sif.obs, xco.obs, 
                             xno2.obs, oco.sensor, xco2.qf = T, xco.qa = 0.4, 
                             xno2.qa = 0.7, zoom = 8, plot.dir = NULL) {

    library(ggpubr)
    tropomi.hr <- unique(substr(xco.obs$time_utc, 9, 10))
    oco.hr <- unique(substr(xco2.obs$timestr, 9, 10))

    # plot google map
    m1 <- ggplot.map(map = 'ggmap', zoom = zoom, center.lat = lon.lat$citylat,
                     center.lon = lon.lat$citylon)[[1]] 
    col <- def.col()[-c(1, length(def.col()))]
    #col <- rev(brewer.pal(11, 'RdYlBu'))

    ### --------------------------- plot xCO
    df1 <- xco.obs %>% filter(qa >= xco.qa) 
    zero.indx <- which(df1$corner == 0)
    df1$group <- findInterval(as.numeric(rownames(df1)), zero.indx)

    c1 <- m1 + theme_bw() + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_polygon(data = df1, 
                       aes(lons, lats, fill = xco * 6.02214 * 1E19, group = group), 
                       alpha = 0.7, color = 'white', size = 0.5) + 
          scale_fill_gradientn(name = 'XCO [molec cm-2]', colours = col) +
          labs(title = paste0('TROPOMI XCO [QA >= ', xco.qa, ']\non ', 
                              substr(timestr, 1, 8), ' ', tropomi.hr, ' UTC')) + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                legend.key.width = unit(1.2, 'cm'))


    ### --------------------------- plot XCO2 using vertex lat/lon 
    if (xco2.qf) xco2.obs <- xco2.obs %>% filter(qf == 0)
    c2 <- m1 + theme_bw() + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_polygon(data = xco2.obs, aes(lons, lats, fill = xco2, group = indx), 
                       alpha = 0.7, color = NA, size = 0.5) + 
          scale_fill_gradientn(name = 'XCO2 [ppm]', colours = col) +
          labs(title = paste(oco.sensor, 'XCO2 [QF = 0] for', site, '\non', 
                             substr(timestr, 1, 8), oco.hr, 'UTC')) + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                legend.key.width = unit(1.2, 'cm'))


    ### --------------------------- plot tropo xNO2
    df3 <- xno2.obs %>% filter(qa >= xno2.qa) %>% 
           mutate(tropo_xno2 = ifelse(tropo_xno2 < 0, 0, tropo_xno2)) 
    zero.indx <- which(df3$corner == 0)
    df3$group <- findInterval(as.numeric(rownames(df3)), zero.indx)

    n1 <- m1 + theme_bw() + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_polygon(data = df3, aes(lons, lats, fill = tropo_xno2 * 6.02214 * 1E19, group = group), 
                       alpha = 0.7, color = 'white', size = 0.5) + 
          scale_fill_gradientn(name = 'Tropospheric NO2\n[molec cm-2]', colours = col) +
          labs(title = paste0('Tropospheric NO2 [QA >= ', xno2.qa, ']\non ', 
                              substr(timestr, 1, 8), ' ', tropomi.hr, ' UTC')) + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                legend.key.width = unit(1.2, 'cm'))

    ### plot sfc altitude 
    h1 <- m1 + theme_bw() +
          geom_polygon(data = df3, aes(lons, lats, fill = hsfc, group = group), 
                       alpha = 0.7, color = 'white', size = 0.5) + 
          scale_fill_gradientn(name = 'Zsfc [m]', colours = col) +
          labs(title = 'TROPOMI-retrieved Zsfc', x = 'LONGITUDE', y = 'LATITUDE') + 
          theme(legend.position = 'bottom', legend.key.height = unit(0.5, 'cm'), 
                legend.key.width = unit(1.2, 'cm'))

    ccnh <- ggarrange(c2, c1, n1, h1, ncol = 4)

    ### --------------------------- plot SIF and land cover
    height = 6; all = ccnh
    fn <- paste0(oco.sensor, '_XCO2_XCO_tNO2_', site, '_', substr(timestr, 1, 8), '.png')

    if (!is.null(sif.obs)) {
        if (!nrow(sif.obs) == 0) {
            melt.sif <- sif.obs %>% 
                    dplyr::select('timestr', 'lat', 'lon', 'sif757', 'sif771', 'avg.sif') %>% 
                    melt(id.var = c('timestr', 'lat', 'lon'))

            title <- paste(oco.sensor, 'SIF [W/m2/sr/µm] and IGBP for', site, 'on', timestr)
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
    
    if (is.null(plot.dir)) plot.dir <- './plot'
    ggsave(all, filename = file.path(plot.dir, fn), width = 16, height = height)

}


