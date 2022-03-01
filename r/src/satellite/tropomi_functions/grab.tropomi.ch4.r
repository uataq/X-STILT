
### ------------------------------- 
# timestr should be format of YYYYMMDD, if it has > 8 letters, we will simply truncate it 
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
grab.tropomi.ch4 = function(tropomi.path = NULL, timestr = NULL, lon_lat, 
                            getakTF = F, tropomi.fn = NULL) {

    library(ncdf4)
    if (is.null(tropomi.fn)) {
        if (nchar(timestr) > 8) timestr = substr(timestr, 1, 8)
        fn = list.files(tropomi.path, paste0('____', timestr), full.names = T, recursive = T)

        if (length(fn) > 1) { 
            tropomi.info = find.tropomi(tropomi.path, timestr, lon_lat)
            fn = tropomi.info$fn
        }
        if (length(fn) == 0) { cat('NO TROPOMI CH4 files found..\n'); return() }

        # multiple files happen if the given lon_lat is large
        if (length(fn) > 1) cat('Found multiple tracks that contain soundings in the given spatial domain...loop over each file to extract obs\n')
    } else fn = tropomi.fn 

    # get dimension first, index starts with ZERO
    all_df = NULL 
    for (f in fn) {  # if multiple files over the same domain

    dat = nc_open(f)
    indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')
    indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')
    corner = ncvar_get(dat, 'PRODUCT/corner')

    ## variables, 13 vertical levels for methane
    lat  = ncvar_get(dat, 'PRODUCT/latitude')
    lon  = ncvar_get(dat, 'PRODUCT/longitude')
    lats = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds')
    lons = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds')

    time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
    timestr = format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), 
                     format = '%Y%m%d%H%M%S')
    time_mtrx = t(replicate(length(indx_across_track), as.numeric(timestr)))

    # A continuous quality descriptor, varying between 0 (no data) and 1 (full quality data). 
    # Recommend to ignore data with qa_value < 0.5"
    qa   = signif(ncvar_get(dat, 'PRODUCT/qa_value'), 1)    # quality assurances
    hsfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude')
    psfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure') / 100   # Pa to hPa

    dp = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/pressure_interval') / 100  # Pa to hPa
    sfc_class = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification')

    sif = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/fluorescence')  # mol s-1 m-2 nm-1 sr-1


    # Methane in column averaged dry air mixing ratio of methane, 1E-9, ppb
    ch4 = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio')     
    ch4_bc = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_bias_corrected')
    ch4_uncert = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_precision')

    # Vertically integrated H2O column from weak or strong bands
    h2o_weak = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/water_weak_twoband_total_column')

    h2o_strong = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/water_strong_twoband_total_column')

    # assign proper dimensions to variables
    dimnames(lat) = dimnames(lon) = dimnames(qa) = dimnames(time_mtrx) = 
    dimnames(hsfc) = dimnames(psfc) = dimnames(sfc_class) = dimnames(dp) = 
    dimnames(ch4) = dimnames(ch4_bc) = dimnames(ch4_uncert) = dimnames(sif) = 
    dimnames(h2o_weak) = dimnames(h2o_strong) = 
        list(indx_across_track, indx_along_track)

    dimnames(lats) = dimnames(lons) = 
        list(corner, indx_across_track, indx_along_track)

    # -------------------------  merge matrix 1 ---------------------------
    var_name = list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                      'psfc', 'dp', 'sfc_class', 'ch4', 'ch4_bc', 'ch4_uncert', 
                      'h2o_weak', 'h2o_strong', 'sif'))
    var_list = list(lat, lon, time_mtrx, qa, hsfc, psfc, dp, sfc_class, ch4, 
                    ch4_bc, ch4_uncert, h2o_weak, h2o_strong, sif)    
    
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] = c('indx_across_track', 'indx_along_track')
    sel_df = var_df %>% filter(center_lon >= lon_lat$minlon, 
                               center_lon <= lon_lat$maxlon,
                               center_lat >= lon_lat$minlat, 
                               center_lat <= lon_lat$maxlat, !is.na(ch4_bc)) 
    
    # ------------------------- merge matrix 2 -----------------------------
    loc_name = list(c('lats', 'lons')); loc_list = list(lats, lons)
    loc_array = array(
        data = do.call(cbind, loc_list), 
        dim = c(dim(loc_list[[1]]), length(loc_list)), 
        dimnames = c(dimnames(loc_list[[1]]), loc_name)
    )
    
    loc_df = dcast(melt(loc_array), Var1 + Var2 + Var3~Var4, 
                   value.var = 'value')
    colnames(loc_df)[1:3] = c('corner', 'indx_across_track', 'indx_along_track')
    
    # for lat/lon in the corners, allow some extra buffer, e.g., 0.1 deg here
    loc_df = loc_df %>% filter(!is.na(lats), lons >= lon_lat$minlon - 0.1, 
                                             lons <= lon_lat$maxlon + 0.1,
                                             lats >= lon_lat$minlat - 0.1, 
                                             lats <= lon_lat$maxlat + 0.1) 
    
    # ----------------------- if grab CH4 AK --------------------------------
    if (getakTF) {

        cat('reading CH4 column averaging kernel...\n')

        # Column averaging kernel for the methane retrieval, unitless
        ak = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel') 

        # 12 layers and 13 levels for CH4 retrieval
        layer = ncvar_get(dat, 'PRODUCT/layer')    
        dimnames(ak) = list(layer, indx_across_track, indx_along_track)
        
        ak.df = melt(ak) %>% rename(hgt = Var1, indx_across_track = Var2, 
                                    indx_along_track = Var3, ak = value) 
        
        # will calculate the surface normalized AK and 
        sel_df = ak.df %>% na.omit() %>% 
                 left_join(var_df, 
                            by = c('indx_across_track', 'indx_along_track')) %>%
                 filter(center_lon >= lon_lat$minlon, 
                        center_lon <= lon_lat$maxlon,
                        center_lat >= lon_lat$minlat, 
                        center_lat <= lon_lat$maxlat) %>% 
                 ungroup()
    }   # end if
    nc_close(dat)

    # no need for conversion, since XCH4 is already in column mixing ratio [ppb]
    merge_df = right_join(sel_df, loc_df, by = c('indx_across_track', 
                                                 'indx_along_track')) %>% 
               filter(!is.na(ch4_bc))

    #return(merge_df)
    all_df = rbind(all_df, merge_df)
    }

    return(all_df)
}
# end of subroutinr

