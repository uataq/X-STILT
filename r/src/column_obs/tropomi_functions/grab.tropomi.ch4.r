
### ------------------------------- 
# timestr should be format of YYYYMMDD, if it has > 8 letters, we will simply truncate it 
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
# add modeled wind data, DW, 03/13/2023 

grab.tropomi.ch4 = function(tropomi.path = NULL, timestr = NULL, lon_lat, 
                            getakTF = F, tropomi.fn = NULL) {
    
    library(ncdf4); library(dplyr)
    if (is.null(tropomi.fn)) {
        tropomi.info = find.tropomi(tropomi.path, timestr, lon_lat, FALSE)
        fn = tropomi.info$fn
        if (length(fn) == 0) { cat('NO TROPOMI CH4 files found..\n'); return() }
    } else fn = tropomi.fn 
    
    if ( unique(grepl('s5p_l2_ch4', fn)) ) 
        df = load.tropomi.ch4.day(fn, lon_lat)

    if ( unique(grepl('S5P_', fn)) ) 
        df = load.tropomi.ch4.hh(fn, getakTF, lon_lat)

    return(df)
}


### ------------------------------- 
# official TROPOMI version, files by every 1-2 hours
load.tropomi.ch4.hh = function(tropomi.fn, getakTF = F, lon_lat) {

    # get dimension first, index starts with ZERO
    all_df = NULL 
    for (f in tropomi.fn) {  # if multiple files over the same domain

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
        qa   = signif(ncvar_get(dat, 'PRODUCT/qa_value'), 1)    # QA
        hsfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude')
        psfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure') / 100   # Pa to hPa

        dp = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/pressure_interval') / 100  # Pa to hPa
        sfc_class = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification')

        sif = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/fluorescence')  # mol s-1 m-2 nm-1 sr-1

        # grab modeled wind at 10 m from ECMWF m s-1
        u10 = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/eastward_wind')
        v10 = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/northward_wind')

        # Methane in column averaged dry air mixing ratio of methane, 1E-9, ppb
        xch4 = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio')     
        xch4_bc = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_bias_corrected')
        xch4_uncert = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_precision')

        # Vertically integrated H2O column from weak or strong bands
        xh2o_weak = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/water_weak_twoband_total_column')

        xh2o_strong = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/water_strong_twoband_total_column')

        # assign proper dimensions to variables
        dimnames(lat) = dimnames(lon) = dimnames(qa) = dimnames(time_mtrx) =    dimnames(hsfc) = dimnames(psfc) = dimnames(sfc_class) = dimnames(dp) = dimnames(u10) = dimnames(v10) = dimnames(xch4) = dimnames(xch4_bc) = dimnames(xch4_uncert) = dimnames(sif) = dimnames(xh2o_weak) = dimnames(xh2o_strong) = list(indx_across_track, indx_along_track)

        dimnames(lats) = dimnames(lons) = 
            list(corner, indx_across_track, indx_along_track)

        # -------------------------  merge matrix 1 ---------------------------
        var_name = list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                        'psfc', 'dp', 'sfc_class', 'u10_ecmwf', 'v10_ecmwf','xch4', 'xch4_bc', 'xch4_uncert', 'xh2o_weak', 'xh2o_strong', 'sif'))
        var_list = list(lat, lon, time_mtrx, qa, hsfc, psfc, dp, sfc_class, 
                        u10, v10, xch4, xch4_bc, xch4_uncert, xh2o_weak, xh2o_strong, sif)
        
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
                                   center_lat <= lon_lat$maxlat, 
                                   !is.na(xch4_bc)) 
        
        # ------------------------- merge matrix 2 -----------------------------
        # lon/lat for the polygon corners
        loc_name = list(c('lats', 'lons')); loc_list = list(lats, lons)
        loc_array = array(
            data = do.call(cbind, loc_list), 
            dim = c(dim(loc_list[[1]]), length(loc_list)), 
            dimnames = c(dimnames(loc_list[[1]]), loc_name)
        )
        
        loc_df = dcast(melt(loc_array), Var1 +Var2 +Var3~Var4, value.var = 'value')
        colnames(loc_df)[1:3] = c('corner', 'indx_across_track', 'indx_along_track')
        
        # for lat/lon in the corners, allow some extra buffer, e.g., 0.1 deg
        loc_df = loc_df %>% filter(!is.na(lats), lons >= lon_lat$minlon - 0.1, 
                                                 lons <= lon_lat$maxlon + 0.1,
                                                 lats >= lon_lat$minlat - 0.1, 
                                                 lats <= lon_lat$maxlat + 0.1) 
        
        # --------------------------------------------------------------------
        # always store the surface AK
        cat('reading CH4 column averaging kernel...\n')

        # Column averaging kernel for the methane retrieval, unitless
        ak = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel') 

        # 12 layers and 13 levels for CH4 retrieval, from TOA to surface
        layer = ncvar_get(dat, 'PRODUCT/layer')   # 0 to 11
        hgts = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/altitude_levels')
        dimnames(ak) = list(layer, indx_across_track, indx_along_track)

        # level 1  <--- Layer 1 ---> level 2, ...downwards...
        # hgt will be the lower hgt for a layer
        ak_df = melt(ak) %>% rename(layer = Var1, indx_across_track = Var2, 
                                    indx_along_track = Var3, ak = value) %>% 
                mutate(lower_level = layer + 2) %>% na.omit()

        hgt_df = melt(hgts) %>% na.omit() %>% 
                 rename(lower_level = Var1, indx_across_track = Var2,
                        indx_along_track = Var3, hgt = value)

        ak_df = left_join(ak_df, hgt_df, 
                          by = c('lower_level', 'indx_across_track', 'indx_along_track'))

        ak_sfc = ak_df %>% filter(lower_level == max(lower_level)) %>% 
                 dplyr::select(-c(layer, lower_level, hgt)) %>% 
                 rename(ak_sfc = ak) 
        
        sel_df = sel_df %>% 
                  left_join(ak_sfc, by = c('indx_across_track', 'indx_along_track'))

        # -------------------- if grab CH4 AK --------------------------------
        if (getakTF) {

            # will calculate the surface normalized AK and 
            sel_df = sel_df %>% left_join(var_df, 
                            by = c('indx_across_track', 'indx_along_track')) %>%
                    filter(center_lon >= lon_lat$minlon, 
                            center_lon <= lon_lat$maxlon,
                            center_lat >= lon_lat$minlat, 
                            center_lat <= lon_lat$maxlat) %>% 
                    ungroup()
        }   # end if
        nc_close(dat)

        merge_df = sel_df %>% left_join(loc_df, by = c('indx_across_track', 'indx_along_track'))
        
        all_df = rbind(all_df, merge_df)
    }

    return(all_df)
}
# end of subroutinr



# if filename contains 's5p_l2_ch4..'
# TROPOMI version from SNOR, files by each day
load.tropomi.ch4.day = function(tropomi.fn, lon_lat) {

    all_df = NULL 
    for ( fn in tropomi.fn) {

        dat = nc_open(fn)
        lat = ncvar_get(dat, 'instrument/latitude_center')
        lon = ncvar_get(dat, 'instrument/longitude_center')
        xch4 = ncvar_get(dat, 'target_product/xch4')
        xch4_uncert = ncvar_get(dat, 'target_product/xch4_precision')
        qa = ncvar_get(dat, 'diagnostics/qa_value')
        time_mtrx = ncvar_get(dat, 'instrument/time')
        time_df = as.data.frame(t(time_mtrx)) %>% 
                  mutate(time_utc = paste0(V1, formatC(V2, width = 2, flag = 0),
                                           formatC(V3, width = 2, flag = 0), 
                                           formatC(V4, width = 2, flag = 0),
                                           formatC(V5, width = 2, flag = 0), 
                                           formatC(V6, width = 2, flag = 0)))

        df = data.frame(center_lon = lon, center_lat = lat, 
                        time_utc = time_df$time_utc, xch4 = xch4, 
                        xch4_uncert = xch4_uncert, qa = qa) 
        nc_close(dat)
        all_df = rbind(all_df, df)
    }

    sel_df =  all_df %>% filter(lon <= 1e36, 
                                lon >= lon_lat$minlon - 0.1, 
                                lon <= lon_lat$maxlon + 0.1,
                                lat >= lon_lat$minlat - 0.1, 
                                lat <= lon_lat$maxlat + 0.1) %>% 
              mutate(xch4 = ifelse(xch4 > 1e36, NA, xch4), 
                     xch4_uncert = ifelse(xch4_uncert > 1e36, NA, xch4_uncert), 
                     qa = ifelse(qa > 1e36, NA, qa))
        
	return(sel_df)
}	

