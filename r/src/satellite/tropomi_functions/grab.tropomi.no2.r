

### ------------------------------- 
# timestr needs to have a format of YYYYMMDD
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
grab.tropomi.no2 = function(tropomi.path = NULL, timestr = NULL, lon.lat, 
                            tropomi.fn = NULL) {

    library(ncdf4)
    if (is.null(tropomi.fn)) {
        
        if (nchar(timestr) > 8) timestr = substr(timestr, 1, 8)
        fn = list.files(tropomi.path, paste0('____', timestr), 
                        full.names = T, recursive = T)

        if (length(fn) > 1) { 
            cat('grab.tropomi.no2(): Multiple TROPOMI files;\nLooking for the tropomi file that has soundings over your spatial domain\n')
            tropomi.info = find.tropomi(tropomi.path, timestr, lon.lat)
            fn = tropomi.info$fn
        }
        if (length(fn) == 0) { cat('NO TROPOMI NO2 files found..\n'); return() }
    } else fn = tropomi.fn 

    # -------------------------------------------------------------------------
    # get dimension first, index starts with ZERO
    dat = nc_open(fn)
    indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')
    indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')
    corner = ncvar_get(dat, 'PRODUCT/corner')

    ## variables
    lat  = ncvar_get(dat, 'PRODUCT/latitude')  # centered lat/lon
    lon  = ncvar_get(dat, 'PRODUCT/longitude')
    lats = ncvar_get(dat, 'GEOLOCATIONS/latitude_bounds')
    lons = ncvar_get(dat, 'GEOLOCATIONS/longitude_bounds')

    time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
    timestr = format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), format = '%Y%m%d%H%M%S')
    time_mtrx = t(replicate(length(indx_across_track), as.numeric(timestr)))
    qa   = signif(ncvar_get(dat, 'PRODUCT/qa_value'), 1)  # quality assurances
    hsfc = ncvar_get(dat, 'INPUT_DATA/surface_altitude')
    psfc = ncvar_get(dat, 'INPUT_DATA/surface_pressure') / 100 # Pa to hPa
    amf_tot = ncvar_get(dat, 'PRODUCT/air_mass_factor_total')
    amf_tropo = ncvar_get(dat, 'PRODUCT/air_mass_factor_troposphere')

    # "TM5 layer index of the highest layer in the tropopause"
    tropo_indx = ncvar_get(dat, 'PRODUCT/tm5_tropopause_layer_index')
    sfc_class  = ncvar_get(dat, 'INPUT_DATA/surface_classification')
    no2_vcd_tropo = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column')
    no2_vcd_tropo_uncert = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision')
    #no2_vcd_tropo_kernel = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision_kernel')

    # assign proper dimensions to variables
    dimnames(lat) = dimnames(lon) = dimnames(qa) = dimnames(time_mtrx) = 
    dimnames(hsfc) = dimnames(psfc) = dimnames(sfc_class) = 
    dimnames(tropo_indx) = dimnames(no2_vcd_tropo) = 
    dimnames(no2_vcd_tropo_uncert) = list(indx_across_track, indx_along_track)
    
    dimnames(lats) = dimnames(lons) = list(corner, indx_across_track, indx_along_track)


    # -------------------------  merge matrix 1 -------------------------------
    var_name = list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                      'psfc', 'sfc_class', 'tropo_k', 'no2_vcd_tropo', 
                      'no2_vcd_tropo_uncert'))
    var_list = list(lat, lon, time_mtrx, qa, hsfc, psfc, sfc_class, tropo_indx,
                    no2_vcd_tropo, no2_vcd_tropo_uncert)    
    
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    library(reshape2)
    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] = c('indx_across_track', 'indx_along_track')
    

    # --------------- select data frame based on spatial domain 
    sel_df = var_df %>% filter(center_lon >= lon.lat$minlon, 
                               center_lon <= lon.lat$maxlon,
                               center_lat >= lon.lat$minlat, 
                               center_lat <= lon.lat$maxlat, !is.na(no2_vcd_tropo)) 
    sel_df = sel_df %>% mutate(sample = 1 : nrow(sel_df)) 
    

    # calculate tropopause pressure, need pressure levels and tropopause level indx
    # TM5 hybrid A coefficient at upper and lower interface levels in Pa, dim of [34, 2]
    a = ncvar_get(dat, 'PRODUCT/tm5_constant_a') / 100   # convert to hPa now
    # TM5 hybrid B coefficient at upper and lower interface levels, unitless, dim of [34, 2]
    b = ncvar_get(dat, 'PRODUCT/tm5_constant_b')   
    # k from surface to top of atmosphere, 0 to 33 
    k = ncvar_get(dat, 'PRODUCT/layer')  
    # vertices: lower or upper pressure level boundaries, v = 0 for lower level
    v = ncvar_get(dat, 'PRODUCT/vertices')  

    # grab pressure coefficient for lower and upper level boundaries
    ab_df = data.frame(a_lower = rep(a[1, ], nrow(sel_df)), 
                       b_lower = rep(b[1, ], nrow(sel_df)), 
                       a_upper = rep(a[2, ], nrow(sel_df)), 
                       b_upper = rep(b[2, ], nrow(sel_df)), 
                       k = rep(k, nrow(sel_df)),   # from sfc to TOA, indx = 0 to 33
                       sample = rep(1 : nrow(sel_df), each = length(a[1, ])))

    add_df = sel_df %>% left_join(ab_df, by = 'sample') %>% group_by(sample) %>%
                        filter(k == tropo_k) %>% ungroup() %>% 
                        mutate(tropo_lower_pres = a_lower + b_lower * psfc, 
                               tropo_upper_pres = a_upper + b_upper * psfc)
    #print(range(add_df$tropo_lower_pres))
    nc_close(dat)


    # ------------------------- merge matrix 2 -------------------------------
    loc_name = list(c('lats', 'lons')); loc_list = list(lats, lons)
    loc_array = array(
        data = do.call(cbind, loc_list), 
        dim = c(dim(loc_list[[1]]), length(loc_list)), 
        dimnames = c(dimnames(loc_list[[1]]), loc_name)
    )
    
    loc_df = dcast(melt(loc_array), Var1 + Var2 + Var3~Var4, value.var = 'value')
    colnames(loc_df)[1:3] = c('corner', 'indx_across_track', 'indx_along_track')
    loc_df = loc_df %>% filter(lons >= lon.lat$minlon - 0.05, 
                               lons <= lon.lat$maxlon + 0.05,
                               lats >= lon.lat$minlat - 0.05, 
                               lats <= lon.lat$maxlat + 0.05, !is.na(lats))

    merge_df = right_join(add_df, loc_df, by = c('indx_across_track', 
                                                 'indx_along_track')) %>% 
               filter(!is.na(no2_vcd_tropo))

    # for TROPOMI vertical column density (VCD) convert mole m-2 to ppb first
    # use surface pressure and water column density 
    # to calculate the total column density of dry air in mol m-2 
    g = 9.8        	         # m s-2
    Mair = 29 / 1E3        	 # kg mol-1

    # TROPOMI pres levels: p_lower = a_lower + b_lower * psfc -> psfc = 0 + 1 * psfc
    # so for tropospheric atmosphere, use both lower bound of tropopause layer
    # and the lower bound for the surface layer

    # *** DW, 03/09/2021, no VCD of H2O vapor, only slant column density, 
    # ignore H2Ov, which is not a perfect idea, need to fix it in the future ***
    merge_df = merge_df %>% 
               mutate(air_vcd = 1 / g / Mair * psfc * 100, 
                      air_vcd_tropo = 1 / g / Mair * (psfc - tropo_lower_pres) * 100, 

                      # final tropospheric column in ppb
                      tno2 = no2_vcd_tropo / air_vcd_tropo * 1E9, 
                      tno2_uncert = no2_vcd_tropo_uncert / air_vcd_tropo * 1E9)

    return(merge_df)
}
# end of subroutinr







### ------------------------------- 
# timestr needs to have a format of YYYYMMDD
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
grab.tropomi.no2.lite = function(tropomi.path = NULL, timestr = NULL, lon.lat, 
                                  tropomi.fn = NULL) {

    library(ncdf4)
    if (is.null(tropomi.fn)) {
        
        if (nchar(timestr) > 8) timestr = substr(timestr, 1, 8)
        fn = list.files(tropomi.path, paste0('____', timestr))

        if (length(fn) > 1) { 
            cat('grab.tropomi.no2(): Multiple TROPOMI files;\nLooking for the tropomi file that has soundings over your spatial domain\n')
            tropomi.info = find.tropomi(tropomi.path, timestr, lon.lat)
            fn = tropomi.info$fn
        }
        if (length(fn) == 0) { cat('NO TROPOMI NO2 files found..\n'); return() }
        fn = file.path(tropomi.path, fn)

    } else fn = tropomi.fn 

    # get dimension first, index starts with ZERO
    dat = nc_open(fn)
    indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')
    indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')
    corner = ncvar_get(dat, 'PRODUCT/corner')

    ## variables
    lat  = ncvar_get(dat, 'PRODUCT/latitude')  # centered lat/lon
    lon  = ncvar_get(dat, 'PRODUCT/longitude')
    lats = ncvar_get(dat, 'GEOLOCATIONS/latitude_bounds')
    lons = ncvar_get(dat, 'GEOLOCATIONS/longitude_bounds')

    time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
    timestr = format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), format = '%Y%m%d%H%M%S')
    time_mtrx = t(replicate(length(indx_across_track), as.numeric(timestr)))
    qa   = ncvar_get(dat, 'PRODUCT/qa_value')  # quality assurances
    hsfc = ncvar_get(dat, 'INPUT_DATA/surface_altitude')
    psfc = ncvar_get(dat, 'INPUT_DATA/surface_pressure') / 100 # Pa to hPa
    amf_tot = ncvar_get(dat, 'PRODUCT/air_mass_factor_total')
    amf_tropo = ncvar_get(dat, 'PRODUCT/air_mass_factor_troposphere')

    # "TM5 layer index of the highest layer in the tropopause"
    tropo_indx = ncvar_get(dat, 'PRODUCT/tm5_tropopause_layer_index')
    sfc_class = ncvar_get(dat, 'INPUT_DATA/surface_classification')
    no2_vcd_tropo = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column')
    no2_vcd_tropo_uncert = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision')
    #no2_vcd_tropo_kernel = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column_precision_kernel')
    nc_close(dat)

    # assign proper dimensions to variables
    dimnames(lat) = dimnames(lon) = dimnames(qa) = dimnames(time_mtrx) = 
    dimnames(hsfc) = dimnames(psfc) = dimnames(sfc_class) = 
    dimnames(tropo_indx) = dimnames(no2_vcd_tropo) = 
    dimnames(no2_vcd_tropo_uncert) = list(indx_across_track, indx_along_track)
    
    # -------------------------  merge matrix 1
    var_name = list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                      'psfc', 'sfc_class', 'tropo_k', 'no2_vcd_tropo', 
                      'no2_vcd_tropo_uncert'))
    var_list = list(lat, lon, time_mtrx, qa, hsfc, psfc, sfc_class, tropo_indx,
                    no2_vcd_tropo, no2_vcd_tropo_uncert)    
    
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] = c('indx_across_track', 'indx_along_track')
    
    # select data frame based on spatial domain 
    sel_df = var_df %>% filter(center_lon >= lon.lat$minlon, 
                               center_lon <= lon.lat$maxlon,
                               center_lat >= lon.lat$minlat, 
                               center_lat <= lon.lat$maxlat, 
                               !is.na(no2_vcd_tropo)) 

    # ------------------------- merge matrix 2
    dimnames(lats) = dimnames(lons) = list(corner, indx_across_track, indx_along_track)
    loc_name = list(c('lats', 'lons')); loc_list = list(lats, lons)
    loc_array = array(
        data = do.call(cbind, loc_list), 
        dim = c(dim(loc_list[[1]]), length(loc_list)), 
        dimnames = c(dimnames(loc_list[[1]]), loc_name)
    )
    
    loc_df = dcast(melt(loc_array), Var1 + Var2 + Var3~Var4, value.var = 'value')
    colnames(loc_df)[1:3] = c('corner', 'indx_across_track', 'indx_along_track')
    
    # for lat/lon in the corners, allow some extra buffer, e.g., 0.1 deg here
    loc_df = loc_df %>% filter(!is.na(lats), lons >= lon.lat$minlon - 0.05, 
                                              lons <= lon.lat$maxlon + 0.05,
                                              lats >= lon.lat$minlat - 0.05, 
                                              lats <= lon.lat$maxlat + 0.05)

    merge_df = right_join(sel_df, loc_df, by = c('indx_across_track', 
                                                  'indx_along_track')) %>% 
                filter(!is.na(no2_vcd_tropo))
    

    return(merge_df)
}
# end of subroutinr

