
# Read S5P-PAL L3 daily gridded TCWV (HARP format, 0.022 deg global grid), 
# DW, 11/09/2026. Files are produced with harpmerge from all L2 orbits of one 
# day, keeping only pixels with qa == 1 (validity > 50), and may have no .nc 
# extension (e.g. a UUID name). Variables: water_vapor_column_density [kg/m2],
# water_vapor_column_density_uncertainty, weight, datetime_start/stop (days 
# since 2000-01-01, whole-day span), latitude, longitude (cell centers).
#
# !!! no per-pixel time in L3. time_utc is ESTIMATED from the S5P ~13:30 local 
#     solar time overpass: UTC hour = 13.5 - lon / 15 (good to about +/- 1 h 
#     across the swath). Use L2 if exact overpass time matters.
# !!! no surface pressure in L3, so no XH2O; only tcwv (kg m-2) and h2o_vcd 
#     (mol m-2) are returned. qa is set to 1 since the L3 only kept qa == 1.
# Output mimics the L2 reader: 4 rows (corner 0..3) per cell built from the 
# grid spacing, center_lon/center_lat, numeric time_utc, qa, tcwv, ...

grab.tropomi.tcwv.l3 = function(tropomi.path = NULL, timestr = NULL, lon_lat, 
                                tropomi.fn = NULL) {
    
    library(ncdf4)
    if (is.null(tropomi.fn)) {      # pick the file(s) whose day matches timestr
        if (is.null(timestr)) stop('grab.tropomi.tcwv.l3(): need timestr or tropomi.fn\n')
        fns = list.files(tropomi.path, full.names = T)
        fn = NULL
        for (f in fns) {
            dat = try(nc_open(f), silent = T)
            if (inherits(dat, 'try-error')) next
            if ('datetime_start' %in% names(dat$var)) {
                day = format(as.POSIXct(ncvar_get(dat, 'datetime_start') * 86400, 
                                        'UTC', origin = '2000-01-01'), '%Y%m%d')
                if (day == substr(timestr, 1, 8)) fn = c(fn, f)
            }
            nc_close(dat)
        }
        if (length(fn) == 0) {cat('NO TROPOMI L3 TCWV files found..\n'); return()}
    } else fn = tropomi.fn 

    Mh2o = 18 / 1E3         # kg mol-1
    all_df = NULL
    for (f in fn) {

        dat = nc_open(f)
        lat = ncvar_get(dat, 'latitude')        # cell centers
        lon = ncvar_get(dat, 'longitude')
        dlat = abs(diff(lat[1:2])); dlon = abs(diff(lon[1:2]))
        ix = which(lon >= lon_lat$minlon & lon <= lon_lat$maxlon)
        iy = which(lat >= lon_lat$minlat & lat <= lon_lat$maxlat)
        if (length(ix) == 0 | length(iy) == 0) { nc_close(dat); next }

        # fields are [longitude, latitude, time] in R
        st = c(min(ix), min(iy), 1); ct = c(length(ix), length(iy), 1)
        get = function(v) ncvar_get(dat, v, start = st, count = ct)
        tcwv = get('water_vapor_column_density')
        tcwv_uncert = get('water_vapor_column_density_uncertainty')
        wgt  = get('weight')
        day0 = as.POSIXct(ncvar_get(dat, 'datetime_start') * 86400, 'UTC', 
                          origin = '2000-01-01')
        day0 = as.POSIXct(format(day0, '%Y-%m-%d'), 'UTC')
        nc_close(dat)

        # matrices are [lon, lat]; as.vector() runs lon fastest
        grd = expand.grid(center_lon = lon[ix], center_lat = lat[iy])
        var_df = data.frame(grd, tcwv = as.vector(tcwv), 
                            tcwv_uncert = as.vector(tcwv_uncert), 
                            weight = as.vector(wgt)) %>% 
                 filter(!is.na(tcwv)) %>% 
                 mutate(qa = 1, h2o_vcd = tcwv / Mh2o, 
                        h2o_vcd_uncert = tcwv_uncert / Mh2o, 
                        # estimated overpass time, ~13:30 local solar time
                        date = day0 + (13.5 - center_lon / 15) * 3600, 
                        time_utc = as.numeric(format(date, '%Y%m%d%H%M%S'))) %>% 
                 dplyr::select(-date)

        # square corners (SW, SE, NE, NW) from grid spacing, corner 0..3
        loc_df = var_df %>% dplyr::select(center_lon, center_lat) %>% 
                 tidyr::crossing(corner = 0:3) %>% 
                 mutate(lons = center_lon + c(-1, 1, 1, -1)[corner + 1] * dlon / 2, 
                        lats = center_lat + c(-1, -1, 1, 1)[corner + 1] * dlat / 2)

        merge_df = loc_df %>% left_join(var_df, by = c('center_lon', 'center_lat'))
        all_df = rbind(all_df, merge_df)
    }   # end for f

    # receptor time as POSIXlt in UTC (moved from get.sensor.obs), DW, 11/09/2026
    if (!is.null(all_df)) 
        all_df$datestr = as.POSIXlt(as.character(all_df$time_utc), 'UTC', 
                                   format = '%Y%m%d%H%M%S')
    closeAllConnections()
    return(all_df)
}
# end of subroutine
