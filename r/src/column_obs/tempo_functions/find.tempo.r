### -------------------------------
# script to find TEMPO data given spatial domain and day, mimic find.tropomi(),
# DW, 10/20/2026
# timestr needs to have a format of YYYYMMDD or YYYYMMDDHH, only one day at a time !!!
# works for both L2 (granule) and L3 (scan) files of NO2 or HCHO, e.g.,
#   TEMPO_NO2_L2_V04_20240725T103050Z_S001G01.nc
#   TEMPO_NO2_L3_V04_20240728T103139Z_S001.nc

# this function finds all granules/scans that have soundings in lon_lat
# one need to select the scan if you only want one with most sounding numbers
# TEMPO QA: main_data_quality_flag, 0 normal, 1 suspicious, 2 bad

if (F) {
    tempo_path = obs_path
}

find.tempo = function(tempo_path, timestr, lon_lat, bufferTF = F){

    library(ncdf4)
    if (nchar(timestr) > 8) {
        YYYYMMDD = substr(timestr, 1, 8)
        HH = substr(timestr, 9, 10)
    } else { YYYYMMDD = timestr; HH = NA }

    fns = list.files(tempo_path, 'TEMPO_.*\\.nc$', full.names = T)
    fns = fns[grepl(paste0('_', YYYYMMDD, 'T'), basename(fns))]

    if (length(fns) == 0) {
        cat('find.tempo(): NO TEMPO files found for the desired days...please check..\n')
        return()
    }

    # start time from file name, end time from global attributes
    # e.g., V2 = species, V3 = level, V5 = start time, V6 = scan + granule
    fn_df = strsplit.to.df(gsub('.nc', '', basename(fns)))
    fn_df = fn_df %>%
            dplyr::select(species = V2, level = V3, mn = V5, scan = V6) %>%
            mutate(mn = gsub('Z', '', mn),
                   mx = sapply(fns, USE.NAMES = F, function(x) {
                        dat = nc_open(x)
                        mx = ncatt_get(dat, 0, 'time_coverage_end')$value
                        nc_close(dat)
                        format(as.POSIXct(mx, 'UTC',
                                          format = '%Y-%m-%dT%H:%M:%SZ'),
                               format = '%Y%m%dT%H%M%S') }),
                   date_mn = as.POSIXct(mn, 'UTC', format = '%Y%m%dT%H%M%S'),
                   date_mx = as.POSIXct(mx, 'UTC', format = '%Y%m%dT%H%M%S'))

    # if hour is given, further select some of the TEMPO files that
    # overlap with [HH:00, HH:59]
    if ( !is.na(HH) ) {

        cat('\nfind.tempo(): HH is given, narrowing down TEMPO files...\n')
        datestr = as.POSIXct(paste0(YYYYMMDD, HH), 'UTC', format = '%Y%m%d%H')
        exact_indx = which(fn_df$date_mn <= datestr + 59 * 60 &
                           fn_df$date_mx >= datestr)

        if (length(exact_indx) == 0) {
            cat('Cannot find any TEMPO files that match your time...\n')
            return()
        }

        # add buffer by 1 nc file to make sure, +/- file before and after
        if (bufferTF) {
            indices = c(min(exact_indx) - 1, exact_indx, max(exact_indx) + 1)
            indices = indices[indices > 0 & indices <= length(fns)]
        } else indices = exact_indx

        fns = fns[indices]
        fn_df = fn_df[indices, ]
    }   # end if HH selection


    tempo_df = NULL
    for (f in 1 : length(fns)) {

        if (f == 2)
            cat('find.tempo(): Searching for the correct TEMPO file...\n')
        if (f %% 3 == 0)
            cat(paste('# ----', signif(f / length(fns) * 100, 3), '% --- #\n'))

        fn = fns[f]
        dat = nc_open(fn)

        # vertical column name differs among species
        if (fn_df$species[f] == 'NO2') {
            val_name = 'product/vertical_column_troposphere'
        } else if (fn_df$species[f] == 'HCHO') {
            val_name = 'product/vertical_column'
        } else stop('find.tempo(): only NO2 and HCHO are implemented...\n')

        if (fn_df$level[f] == 'L2') {

            # [xtrack, mirror_step]
            lat = ncvar_get(dat, 'geolocation/latitude')
            lon = ncvar_get(dat, 'geolocation/longitude')
            qa  = ncvar_get(dat, 'product/main_data_quality_flag')
            val = ncvar_get(dat, val_name)

            # [mirror_step], sec since 1980-01-06T00:00:00Z
            sec = ncvar_get(dat, 'geolocation/time')
            timestr = format(as.POSIXct(sec, 'UTC', origin = '1980-01-06'),
                             format = '%Y%m%d%H%M%S')
            timestr = matrix(timestr, nrow = nrow(lat), ncol = ncol(lat),
                             byrow = T)

            sel_indx = which(!is.na(val) &
                             lon >= lon_lat$minlon & lon <= lon_lat$maxlon &
                             lat >= lon_lat$minlat & lat <= lon_lat$maxlat)
            sel_df = data.frame(qa = qa[sel_indx], timestr = timestr[sel_indx],
                                stringsAsFactors = F)
            
        } else if (fn_df$level[f] == 'L3') {

            # regular grid, only read the sub-domain, [lon, lat, time]
            lat = ncvar_get(dat, 'latitude')
            lon = ncvar_get(dat, 'longitude')
            lon_indx = which(lon >= lon_lat$minlon & lon <= lon_lat$maxlon)
            lat_indx = which(lat >= lon_lat$minlat & lat <= lon_lat$maxlat)

            if (length(lon_indx) == 0 | length(lat_indx) == 0) {
                sel_df = data.frame()
            } else {
                st = c(min(lon_indx), min(lat_indx), 1)
                ct = c(length(lon_indx), length(lat_indx), 1)
                qa  = ncvar_get(dat, 'product/main_data_quality_flag', st, ct)
                val = ncvar_get(dat, val_name, st, ct)

                # no pixel-level time for L3, use scan start time
                qa = qa[!is.na(val)]
                sel_df = data.frame(qa = qa,
                                    timestr = rep(format(fn_df$date_mn[f],
                                              format = '%Y%m%d%H%M%S'),
                                              length(qa)),
                                    stringsAsFactors = F)
            }

        } else stop('find.tempo(): only L2 and L3 files are implemented...\n')
        nc_close(dat)

        if (nrow(sel_df) > 0) { # if found one file, return value
            cat(paste0('found the data that match the criteria; see nc file: ', basename(fn), '\n'))

            overpass.timestr = unique(substr(sel_df$timestr, 1, 14))
            tmp_df = data.frame(start.time = fn_df$mn[f],
                                end.time = fn_df$mx[f],
                                overpass.start.time = min(overpass.timestr),
                                overpass.end.time = max(overpass.timestr),
                                tot.count = nrow(sel_df),
                                qa0.count = nrow(sel_df %>% filter(qa == 0)),
                                qa1.count = nrow(sel_df %>% filter(qa <= 1)),
                                fn = fn, stringsAsFactors = F)

            # if no soundings with QA <= 1, look for the next TEMPO file
            if (tmp_df$qa1.count == 0) {
                cat('find.tempo(): NO qualified data with QA <= 1 found\n')
                next
            }
            tempo_df = rbind(tempo_df, tmp_df)
        }   # end if
    }   # end for fn

    if (is.null(tempo_df)) {
        cat('Cannot find any TEMPO data that match your time and spatial domain...\n')
        return()
    }

    return(tempo_df)
}
# end of subroutine
