

### ------------------------------- 
# timestr needs to have a format of YYYYMMDD or YYYYMMDDHH
# if one provides tropomi.fn, no need to provide tropomi.path and timestr
grab.tropomi.hcho = function(tropomi.path = NULL, timestr = NULL, lon_lat, 
                             tropomi.fn = NULL) {

    library(ncdf4)
    if (is.null(tropomi.fn)) {
        tropomi.info = find.tropomi(tropomi.path, timestr, lon_lat)
        fn = tropomi.info$fn
        if (length(fn) == 0) {cat('NO TROPOMI HCHO files found..\n'); return()}
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
    lats = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds')
    lons = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds')
    time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC

    if (unique(time) == '') {

        # for HCHO  seconds since 2010-01-01 00:00:00
        ref_time = ncvar_get(dat, 'PRODUCT/time') # reference time in second
        
        # offset from reference start time of measurement
        # milliseconds since @param time, [ground_pix, scaline]
        dt   = ncvar_get(dat, 'PRODUCT/delta_time') * 1E-3 # now in second
        time = matrix(ref_time, nrow = nrow(dt), ncol = ncol(dt)) + dt
        dimnames(time) = list(indx_across_track, indx_along_track)
        
        time_mtrx = time
        time_mtrx[] = format(as.POSIXct(time_mtrx[], 'UTC', 
                             origin = '2010-01-01 00:00:00'),
                             format = '%Y%m%d%H%M%S')

    } else {
        timestr = as.POSIXct(time, 'UTC', format = '%Y-%m-%dT%H:%M:%S')  
        timestr = format(timestr, format = '%Y%m%d%H%M%S')
        time_mtrx = t(replicate(length(indx_across_track), as.numeric(timestr)))
    }

    qa   = signif(ncvar_get(dat, 'PRODUCT/qa_value'), 1)  # quality assurances
    hsfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude')
    psfc = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure') / 100 # Pa to hPa

    sfc_class  = ncvar_get(dat, 'PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification')

    hcho_vcd = ncvar_get(dat, 'PRODUCT/formaldehyde_tropospheric_vertical_column')
    
    hcho_vcd_uncert = ncvar_get(dat, 'PRODUCT/formaldehyde_tropospheric_vertical_column_precision')

    # assign proper dimensions to variables
    dimnames(lat) = dimnames(lon) = dimnames(qa) = dimnames(time_mtrx) = 
    dimnames(hsfc) = dimnames(psfc) = dimnames(sfc_class) = 
    dimnames(hcho_vcd) = dimnames(hcho_vcd_uncert) = 
        list(indx_across_track, indx_along_track)
    
    dimnames(lats) = dimnames(lons) = 
        list(corner, indx_across_track, indx_along_track)


    # -------------------------  merge matrix 1 ----------------------------
    var_name = list(c('center_lat', 'center_lon', 'time_utc', 'qa', 'hsfc', 
                      'psfc', 'sfc_class', 'hcho_vcd', 'hcho_vcd_uncert'))
    var_list = list(lat, lon, time_mtrx, qa, hsfc, psfc, sfc_class, 
                    hcho_vcd, hcho_vcd_uncert)    
    
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions
    
    library(reshape2)
    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value') %>% 
             mutate_if(is.character, as.numeric)
    colnames(var_df)[1:2] = c('indx_across_track', 'indx_along_track')
    

    # --------------- select data frame based on spatial domain 
    sel_df = var_df %>% filter(center_lon >= lon_lat$minlon, 
                               center_lon <= lon_lat$maxlon,
                               center_lat >= lon_lat$minlat, 
                               center_lat <= lon_lat$maxlat, !is.na(hcho_vcd)) 
    
    sel_df = sel_df %>% mutate(sample = 1 : nrow(sel_df)) 
    nc_close(dat)

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
    loc_df = loc_df %>% filter(lons >= lon_lat$minlon - 0.05, 
                               lons <= lon_lat$maxlon + 0.05,
                               lats >= lon_lat$minlat - 0.05, 
                               lats <= lon_lat$maxlat + 0.05, !is.na(lats))

    merge_df = right_join(sel_df, loc_df, 
                          by = c('indx_across_track', 'indx_along_track')) %>% 
               filter(!is.na(hcho_vcd))

    # for TROPOMI vertical column density (VCD) convert mole m-2 to ppb first
    # use surface pressure and water column density 
    # to calculate the total column density of dry air in mol m-2 
    g = 9.8        	         # m s-2
    Mair = 29 / 1E3        	 # kg mol-1

    # *** DW, 03/09/2021, no VCD of H2O vapor, only slant column density, 
    # ignore H2Ov, which is not a perfect idea, need to fix it in the future ***
    merge_df = merge_df %>% 
               mutate(air_vcd = 1 / g / Mair * psfc * 100, 

                      # final tropospheric column in ppb
                      xhcho = hcho_vcd / air_vcd * 1E9, 
                      xhcho_uncert = hcho_vcd_uncert / air_vcd * 1E9)

    return(merge_df)
}
# end of subroutinr

