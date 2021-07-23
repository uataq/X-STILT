
### ------------------------------- 
# timestr should be format of YYYYMMDD, if it has > 8 letters, we will simply truncate it 
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
grab.tropomi.co = function(tropomi.path = NULL, timestr = NULL, lon.lat, 
                           getakTF = F, tropomi.fn = NULL) {

    library(ncdf4)
    if (is.null(tropomi.fn)) {
        if (nchar(timestr) > 8) timestr = substr(timestr, 1, 8)
        fn = list.files(tropomi.path, paste0('____', timestr), full.names = T, recursive = T)

        if (length(fn) > 1) { 
            cat('grab.tropomi.co(): Multiple TROPOMI files;\nLooking for the tropomi file that has soundings over your spatial domain\n')
            tropomi.info = find.tropomi(tropomi.path, timestr, lon.lat)
            fn = tropomi.info$fn
        }
        if (length(fn) == 0) { cat('NO TROPOMI CO files found..\n'); return() }
    } else fn = tropomi.fn 

    # get dimension first, index starts with ZERO
    dat = nc_open(fn)
    indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')
    indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')
    corner = ncvar_get(dat, 'PRODUCT/corner')

    ## variables
    lat  = ncvar_get(dat, 'PRODUCT/latitude')
    lon  = ncvar_get(dat, 'PRODUCT/longitude')
    lats = ncvar_get(dat, 'GEOLOCATIONS/latitude_bounds')
    lons = ncvar_get(dat, 'GEOLOCATIONS/longitude_bounds')

    time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
    timestr = format(as.POSIXct(time, format = '%Y-%m-%dT%H:%M:%S'), format = '%Y%m%d%H%M%S')
    time_mtrx = t(replicate(length(indx_across_track), as.numeric(timestr)))

    # A continuous quality descriptor, varying between 0 (no data) and 1 (full quality data). 
    # Recommend to ignore data with qa_value < 0.5 for TROPOMI CO
    qa   = signif(ncvar_get(dat, 'PRODUCT/qa_value'), 1)  # quality assurances
    hsfc = ncvar_get(dat, 'INPUT_DATA/surface_altitude')
    psfc = ncvar_get(dat, 'INPUT_DATA/surface_pressure') / 100 # Pa to hPa
    sfc_class = ncvar_get(dat, 'INPUT_DATA/surface_classification')
    co_vcd = ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column') # mol m-2
    co_vcd_uncert = ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column_precision')
    h2o_vcd = ncvar_get(dat, 'DETAILED_RESULTS/water_total_column')

    # assign proper dimensions to variables
    dimnames(lat) = dimnames(lon) = dimnames(qa) = dimnames(time_mtrx) = 
    dimnames(hsfc) = dimnames(psfc) = dimnames(sfc_class) = dimnames(co_vcd) = 
    dimnames(co_vcd_uncert) = dimnames(h2o_vcd) = list(indx_across_track, indx_along_track)

    dimnames(lats) = dimnames(lons) = list(corner, indx_across_track, indx_along_track)

    # -------------------------  merge matrix 1
    var_name = list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                     'psfc', 'sfc_class', 'co_vcd', 'co_vcd_uncert', 'h2o_vcd'))
    var_list = list(lat, lon, time_mtrx, qa, hsfc, psfc, sfc_class, co_vcd, 
                    co_vcd_uncert, h2o_vcd)    
    
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] = c('indx_across_track', 'indx_along_track')
    sel_df = var_df %>% filter(center_lon >= lon.lat$minlon, 
                               center_lon <= lon.lat$maxlon,
                               center_lat >= lon.lat$minlat, 
                               center_lat <= lon.lat$maxlat, !is.na(co_vcd)) 
    
    # ------------------------- merge matrix 2
    loc_name = list(c('lats', 'lons')); loc_list = list(lats, lons)
    loc_array = array(
        data = do.call(cbind, loc_list), 
        dim = c(dim(loc_list[[1]]), length(loc_list)), 
        dimnames = c(dimnames(loc_list[[1]]), loc_name)
    )
    
    loc_df = dcast(melt(loc_array), Var1 + Var2 + Var3~Var4, value.var = 'value')
    colnames(loc_df)[1:3] = c('corner', 'indx_across_track', 'indx_along_track')
    
    # for lat/lon in the corners, allow some extra buffer, e.g., 0.1 deg here
    loc_df = loc_df %>% filter(!is.na(lats), lons >= lon.lat$minlon - 0.1, 
                                             lons <= lon.lat$maxlon + 0.1,
                                             lats >= lon.lat$minlat - 0.1, 
                                             lats <= lon.lat$maxlat + 0.1) 
    
    # ----------------------- if grab CO AK and convert from matrix to df
    if (getakTF) {
        cat('reading CO column averaging kernel...\n')
        ak = ncvar_get(dat, 'DETAILED_RESULTS/column_averaging_kernel')    # ak in meter
        layer = ncvar_get(dat, 'PRODUCT/layer')    # height in m
        dimnames(ak) = list(layer, indx_across_track, indx_along_track)
        
        ak.df = melt(ak) %>% dplyr::rename(hgt = Var1, indx_across_track = Var2, 
                                           indx_along_track = Var3, ak = value) 
        
        # will calculate the surface normalized AK and 
        sel_df = ak.df %>% na.omit() %>% 
                  left_join(var_df , by = c('indx_across_track', 'indx_along_track')) %>%
                  filter(center_lon >= lon.lat$minlon, center_lon <= lon.lat$maxlon,
                         center_lat >= lon.lat$minlat, center_lat <= lon.lat$maxlat) %>% 
                  ungroup()
    }   # end if
    nc_close(dat)

    merge_df = right_join(sel_df, loc_df, by = c('indx_across_track', 
                                                  'indx_along_track')) %>% 
               filter(!is.na(co_vcd))
    
    # for TROPOMI vertical column density (VCD) convert mole m-2 to ppb first
    # use surface pressure and water column density 
    # to calculate the total column density of dry air in mol m-2 
    g = 9.8        	        # m s-2
    Mair = 29 / 1E3        	 # kg mol-1
    merge_df = merge_df %>% mutate(air_vcd = 1 / g / Mair * psfc * 100 - h2o_vcd, 
                                   xco = co_vcd / air_vcd * 1E9, 
                                   xco_uncert = co_vcd_uncert / air_vcd * 1E9)

    return(merge_df)
}
# end of subroutinr

