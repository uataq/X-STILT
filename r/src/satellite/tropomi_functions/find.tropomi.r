
### ------------------------------- 
# script to find TROPOMI data given spatial domain and day, DW, 03/04/2021
# timestr needs to have a format of YYYYMMDD, only one day at a time !!!

find.tropomi = function(tropomi_path, timestr, lon_lat) {

    library(ncdf4)    
    if (nchar(timestr) > 8) timestr = substr(timestr, 1, 8)
    fns = list.files(tropomi_path, paste0('____', timestr), full.names = T, recursive = T)
    if (length(fns) == 0) { cat('NO TROPOMI files found...please check..\n'); return() }
    
    cat('find.tropomi(): Searching for the correct TROPOMI file...\n')
    tropomi_df = NULL
    for (f in 1 : length(fns)) {
        
        if (f %% 3 == 0) 
            cat(paste('# ---- ', signif(f / length(fns) * 100, 3), '% ---- #\n'))

        fn = fns[f]
        dat = nc_open(fn)
        indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')
        indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')
        lat = ncvar_get(dat, 'PRODUCT/latitude')    # [ground_pixel, scanline, time]
        lon = ncvar_get(dat, 'PRODUCT/longitude')   # [ground_pixel, scanline, time]
        time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
        timestr = format(as.POSIXct(time, 'UTC', format = '%Y-%m-%dT%H:%M:%S'), 
                         format = '%Y%m%d%H%M%S')   # [scanline, time]
        time_df = data.frame(timestr, indx_along_track, stringsAsFactors = F)

        # for CO and NO2, 
        # for CH4, A continuous quality descriptor, varying between 0 (no data) 
        # and 1 (full quality data). Recommend to ignore data with qa_value < 0.5
        qa = ncvar_get(dat, 'PRODUCT/qa_value')  # quality assurances
        
        if (unique(grepl('_CO_', fn))) 
            val = ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column')
        
        # bias corrected dry_atmosphere_mole_fraction_of_methane
        if (unique(grepl('_CH4_', fn)))     # unit of 1e-9
            val = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_bias_corrected')
            #val = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio')

        if (unique(grepl('_NO2_', fn))) 
            val = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column')

        dimnames(lat) = dimnames(lon) = dimnames(qa) = dimnames(val) = 
            list(indx_across_track, indx_along_track)

        # merge matrix 
        var_name = list(c('lat', 'lon', 'qa', 'val'))
        var_list = list(lat, lon, qa, val)
        
        #it is assumed all matrices in the list have equal dimensions
        var_array = array(
            data = do.call(cbind, var_list), 
            dim = c(dim(var_list[[1]]), length(var_list)), 
            dimnames = c(dimnames(var_list[[1]]), var_name)
        )
        
        var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
        colnames(var_df)[1:2] = c('indx_across_track', 'indx_along_track')

        # qa_value > 0.75. For most users this is the recommended pixel filter. 
        # This removes cloud-covered scenes (cloud radiance fraction > 0.5), 
        # part of the scenes covered by snow/ice, errors and problematic retrievals.
        # qa_value > 0.50. This adds the good quality retrievals over clouds and over scenes covered by snow/ice. 
        # Errors and problematic retrievals are still filtered out. 
        # In particular this choice is useful for assimilation and model comparison
        sel_df = var_df %>% filter(!is.na(val), 
                                   lon >= lon_lat$minlon, lon <= lon_lat$maxlon,
                                   lat >= lon_lat$minlat, lat <= lon_lat$maxlat) %>% 
                left_join(time_df, by = 'indx_along_track')
        nc_close(dat)

        if (nrow(sel_df) > 0) { # if found one file, return value
            cat(paste0('found the data that match the criteria; see nc file: ', 
                       basename(fn), '; stop searching...\n'))
            fn_df = strsplit.to.df(basename(fn)) 
            if (unique(grepl('_CO_', fn))) { start.time = fn_df$V10; end.time = fn_df$V11 }
            if (unique(grepl('_CH4_', fn))) { start.time = fn_df$V9; end.time = fn_df$V10 } 
            if (unique(grepl('_NO2_', fn))) { start.time = fn_df$V9; end.time = fn_df$V10 } 
            
            overpass.timestr = unique(substr(sel_df$timestr, 1, 14))
            tmp_df = data.frame(start.time = start.time, end.time = end.time, 
                                overpass.start.time = min(overpass.timestr), 
                                overpass.end.time = max(overpass.timestr), 
                                tot.count = nrow(sel_df), 
                                qa0p5.count = nrow(sel_df %>% filter(qa >= 0.50)), 
                                qa0p4.count = nrow(sel_df %>% filter(qa >= 0.40)), 
                                qa0p7.count = nrow(sel_df %>% filter(qa >= 0.70)), 
                                fn = fn, stringsAsFactors = F)
            
            # if there is no soundings with QA >= 0.4, look for the next TROPOMI file
            # DW, 06/02/2021
            if (tmp_df$qa0p4.count == 0) next 
            tropomi_df = rbind(tropomi_df, tmp_df)
            return(tropomi_df)
        }  # end if
        
    }   # end for fn

    if (is.null(tropomi_df)) {
        cat('Cannot find any TROPOMI data that match your time and spatial domain...\n')
        return()
    }
}
# end of subrouti
