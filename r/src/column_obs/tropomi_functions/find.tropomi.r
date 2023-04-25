
### ------------------------------- 
# script to find TROPOMI data given spatial domain and day, DW, 03/04/2021
# timestr needs to have a format of YYYYMMDD or YYYYMMDDHHmm, only one day at a time !!!

# this function finds all tracks that have soundings in lon_lat, DW, 10/27/2021
# one need to select the track if you only want one with most sounding numbers, DW, 11/05/2021
find.tropomi = function(tropomi_path, timestr, lon_lat, bufferTF = F) {

    library(ncdf4)    
    if (nchar(timestr) > 8) {
        YYYYMMDD = substr(timestr, 1, 8)
        HH = substr(timestr, 9, nchar(timestr))
    } else { YYYYMMDD = timestr; HH = NA } 
    
    fns = list.files(tropomi_path, paste0('___', YYYYMMDD), 
                     full.names = T, recursive = T)
    if (length(fns) == 0) { 
        cat('NO TROPOMI files found...please check..\n'); return() }
    
    # if hour and beyond is given, further select some of the TROPOMI files
    # DW, 08/19/2021
    if ( !is.na(HH) ) {

        cat('\nfind.tropomi(): HH is given, narrowing down TROPOMI files...\n')
        datestr = as.POSIXct(paste0(YYYYMMDD, HH), 'UTC', format = '%Y%m%d%H')
        HH_fns = strsplit.to.df(basename(fns))
        if (unique(HH_fns$V5) == 'CO') {
            HH_fns = HH_fns %>% dplyr::select(mn = V10, mx = V11)
        } else HH_fns = HH_fns %>% dplyr::select(mn = V9, mx = V10)

        HH_fns = HH_fns %>% 
                 mutate(date_mn = as.POSIXct(mn, 'UTC',format ='%Y%m%dT%H%M%S'),
                        date_mx = as.POSIXct(mx, 'UTC',format ='%Y%m%dT%H%M%S'))
        
        # bug fixed for selecting files based on `timestr`
        exact_indx = unique(
                            c(which(datestr >= HH_fns$date_mn &
                                    datestr <= HH_fns$date_mx), 
                              which(datestr + 59 * 60 >= HH_fns$date_mn &
                                    datestr + 59 * 60 <= HH_fns$date_mx))
                            )

        # add buffer by 1 nc file to make sure, +/- file before and after 
        if (bufferTF) {
            indices = c(min(exact_indx) - 1, exact_indx, max(exact_indx) + 1)
            indices = indices[indices != 0 | indices > length(fns)]
        } else indices = exact_indx
        
        fns = fns[indices]
        if (length(fns) == 0) {
            cat('Cannot find any TROPOMI files that match your time...\n')
            return()
        }
        
    }   # end if HH selection
    

    tropomi_df = NULL
    for (f in 1 : length(fns)) {

        if (f == 2) 
            cat('find.tropomi(): Searching for the correct TROPOMI file...\n')
        if (f %% 3 == 0) 
            cat(paste('# ----', signif(f / length(fns) * 100, 3), '% --- #\n'))
        
        fn = fns[f]
        dat = nc_open(fn)
        indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')
        indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')

        # [ground_pixel, scanline, time]
        lat = ncvar_get(dat, 'PRODUCT/latitude')    
        lon = ncvar_get(dat, 'PRODUCT/longitude') 
        time = substr(ncvar_get(dat, 'PRODUCT/time_utc'), 1, 19)  # UTC
        
        if ('' %in% time) {
            # for HCHO  seconds since 2010-01-01 00:00:00
            ref_time = ncvar_get(dat, 'PRODUCT/time') # reference time in second
            
            # offset from reference start time of measurement
            # milliseconds since @param time, [ground_pix, scaline]
            dt   = ncvar_get(dat, 'PRODUCT/delta_time') * 1E-3 # now in second
            time = matrix(ref_time, nrow = nrow(dt), ncol = ncol(dt)) + dt
            dimnames(time) = list(indx_across_track, indx_along_track)
            
            time_df = melt(time) %>% 
                      rename(indx_across_track = Var1, 
                             indx_along_track = Var2, timestr = value) %>% 
                      mutate(timestr = as.POSIXct(timestr, 'UTC', 
                                        origin = '2010-01-01 00:00:00'), 
                             timestr = format(timestr, format = '%Y%m%d%H%M%S'))

        } else {
            timestr = as.POSIXct(time, 'UTC', format = '%Y-%m-%dT%H:%M:%S')  
            timestr = format(timestr, format = '%Y%m%d%H%M%S')
            time_df = data.frame(timestr, indx_along_track, 
                                 stringsAsFactors = F)
        }

        # for CO and NO2, 
        # for CH4, A continuous quality descriptor, varying between 0 (no data) 
        # and 1 (full quality data).Recommend to ignore data with qa_value < 0.5
        qa = ncvar_get(dat, 'PRODUCT/qa_value')  # quality assurances
        
        if (unique(grepl('_CO_', fn))) 
            val = ncvar_get(dat, 'PRODUCT/carbonmonoxide_total_column')
        
        # bias corrected dry_atmosphere_mole_fraction_of_methane
        if (unique(grepl('_CH4_', fn)))     # unit of 1e-9
            val = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio_bias_corrected')
            #val = ncvar_get(dat, 'PRODUCT/methane_mixing_ratio')

        if (unique(grepl('_NO2_', fn))) 
            val = ncvar_get(dat, 'PRODUCT/nitrogendioxide_tropospheric_column')

        if (unique(grepl('_HCHO_', fn))) 
            val = ncvar_get(dat, 'PRODUCT/formaldehyde_tropospheric_vertical_column')

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
                                   lon >= lon_lat$minlon, 
                                   lon <= lon_lat$maxlon,
                                   lat >= lon_lat$minlat, 
                                   lat <= lon_lat$maxlat) 
                                   
        if (ncol(time_df) == 2) 
            sel_df = sel_df %>% left_join(time_df, by = 'indx_along_track')
        
        if (ncol(time_df) > 2) 
            sel_df = sel_df %>% 
                     left_join(time_df, by = c('indx_along_track',          
                                               'indx_across_track'))
        
        nc_close(dat)

        if (nrow(sel_df) > 0) { # if found one file, return value
            cat(paste0('found the data that match the criteria; see nc file: ', 
                       basename(fn), '\n'))
            fn_df = strsplit.to.df(basename(fn)) 
            
            if (unique(grepl('_CO_', fn))) { 
                start.time = fn_df$V10; end.time = fn_df$V11 
            } else if (unique(grepl('_CH4_', fn)) ) { 
                start.time = fn_df$V9; end.time = fn_df$V10 
            } else if (unique(grepl('_NO2_', fn))) { 
                start.time = fn_df$V9; end.time = fn_df$V10 
            } else if (unique(grepl('_HCHO_', fn))) { 
                start.time = fn_df$V8; end.time = fn_df$V9 
            } 
            
            overpass.timestr = unique(substr(sel_df$timestr, 1, 14))
            tmp_df = data.frame(start.time = start.time, end.time = end.time, 
                                overpass.start.time = min(overpass.timestr), 
                                overpass.end.time = max(overpass.timestr), 
                                tot.count = nrow(sel_df), 
                                qa0p5.count = nrow(sel_df %>% filter(qa>= 0.5)),
                                qa0p4.count = nrow(sel_df %>% filter(qa>= 0.4)),
                                qa0p7.count = nrow(sel_df %>% filter(qa>= 0.7)),
                                fn = fn, stringsAsFactors = F)
            
            # if no soundings with QA >= 0.4, look for the next TROPOMI file
            # DW, 06/02/2021
            if (tmp_df$qa0p4.count == 0) next 
            tropomi_df = rbind(tropomi_df, tmp_df)
        }  # end if
    }   # end for fn

    return(tropomi_df)

    if (is.null(tropomi_df)) {
        cat('Cannot find any TROPOMI data that match your time and spatial domain...\n')
        return()
    }
}
# end of subrouti
