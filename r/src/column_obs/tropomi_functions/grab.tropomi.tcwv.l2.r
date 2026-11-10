
# Read TROPOMI (S5P-PAL) L2 total column water vapor (TCWV) files, DW, 11/09/2026
#   file name, e.g., S5P_PAL__L2__TCWV___20260601T194210_20260601T204034_44739_03_010601_20260603T164247.nc
#
# Differences from the other TROPOMI L2 products:
#   - no PRODUCT/time_utc; time = PRODUCT/time (sec since 2010-01-01) + 
#     PRODUCT/delta_time (msec since the day's midnight), one value per scanline
#   - no PRODUCT/corner variable, no surface_altitude / surface_classification
#   - qa_value is stored as ubyte with scale 0.01 (ncdf4 returns 0..1); the 
#     values are discrete {0, 0.01, 0.12, 0.25, 0.5, 1}; qa == 1 recommended
#   - TCWV in kg m-2; also returned as H2O column in mol m-2 and column-mean 
#     mixing ratio XH2O in ppm (dry air column from psfc and H2O column)
#
# Output layout follows grab.tropomi.co(): 4 rows (corner 0..3) per pixel with
# center_lon/center_lat, lons/lats, numeric time_utc (YYYYMMDDHHMMSS), qa, ...
# Only the sub-domain is read from disk (start/count on [ground_pixel, scanline])
# if tropomi.fn is given, tropomi.path and timestr are not needed

if (F) {
    tropomi.path = sensor_path = obs_path
    tropomi.fn = sensor_fn = obs_fn
}

grab.tropomi.tcwv.l2 = function(tropomi.path = NULL, timestr = NULL, lon_lat, 
                                tropomi.fn = NULL) {
    
    library(ncdf4)
    if (is.null(tropomi.fn)) {
        tropomi.info = find.tropomi(tropomi.path, timestr, lon_lat)
        fn = tropomi.info$fn
        if (length(fn) == 0) {cat('NO TROPOMI TCWV files found..\n'); return()}
    } else fn = tropomi.fn 

    g = 9.8                 # m s-2
    Mair = 29 / 1E3         # kg mol-1
    Mh2o = 18 / 1E3         # kg mol-1

    all_df = NULL
    for (f in fn) {
        
        dat = nc_open(f)
        indx_along_track  = ncvar_get(dat, 'PRODUCT/scanline')      # start at 0
        indx_across_track = ncvar_get(dat, 'PRODUCT/ground_pixel')
        corner = 0:3                            # no PRODUCT/corner var in TCWV
        lat = ncvar_get(dat, 'PRODUCT/latitude')        # [ground_pixel, scanline]
        lon = ncvar_get(dat, 'PRODUCT/longitude')

        # index block whose pixel centers fall inside lon_lat; skip if none
        inTF = !is.na(lon) & !is.na(lat) & 
               lon >= lon_lat$minlon & lon <= lon_lat$maxlon & 
               lat >= lon_lat$minlat & lat <= lon_lat$maxlat
        if (!any(inTF)) { nc_close(dat); next }
        ix = which(apply(inTF, 1, any)); ix = min(ix):max(ix)
        iy = which(apply(inTF, 2, any)); iy = min(iy):max(iy)
        nx = length(ix); ny = length(iy)

        # readers: 2D fields are stored as [ground_pixel, scanline, time], 
        # corner fields as [corner, ground_pixel, scanline, time]
        get  = function(v) ncvar_get(dat, v, start = c(min(ix), min(iy), 1), 
                                             count = c(nx, ny, 1))
        get3 = function(v) ncvar_get(dat, v, start = c(1, min(ix), min(iy), 1), 
                                             count = c(4, nx, ny, 1))

        # ----------------------- time, one value per scanline -----------------
        # delta_time units look like 'milliseconds since 2026-06-02 00:00:00'
        dt_unit = ncatt_get(dat, 'PRODUCT/delta_time', 'units')$value
        dt_orig = as.POSIXct(gsub('milliseconds since ', '', dt_unit), 'UTC')
        dt = ncvar_get(dat, 'PRODUCT/delta_time', start = c(min(iy), 1), 
                                                  count = c(ny, 1))
        scan_date = dt_orig + as.numeric(dt) * 1E-3
        scan_time_utc = as.numeric(format(scan_date, '%Y%m%d%H%M%S'))

        # ----------------------- variables ------------------------------------
        qa   = round(get('PRODUCT/qa_value'), 2)   # ubyte x 0.01f; round so 1 is 1
        tcwv = get('PRODUCT/total_column_water_vapor')              # kg m-2
        tcwv_uncert = get('PRODUCT/total_column_water_vapor_precision')

        psfc = get('PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure') / 100  # hPa
        cld_frac = get('PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction')
        cld_pres = get('PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_pressure') / 100 # hPa
        sfc_alb  = get('PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_albedo')
        amf_tot  = get('PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_total')
        sza = get('PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle')
        saa = get('PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_azimuth_angle')
        vza = get('PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_zenith_angle')
        lons = get3('PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds')
        lats = get3('PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds')
        nc_close(dat)

        # ----------------------- pixel-level data frame -----------------------
        # matrices are [ground_pixel, scanline]; as.vector() runs ground_pixel fastest
        var_df = data.frame(
            indx_across_track = rep(indx_across_track[ix], times = ny), 
            indx_along_track  = rep(indx_along_track[iy], each = nx), 
            center_lon = as.vector(lon[ix, iy]), 
            center_lat = as.vector(lat[ix, iy]), 
            time_utc = rep(scan_time_utc, each = nx), 
            qa = as.vector(qa), 
            psfc = as.vector(psfc), 
            sza = as.vector(sza), 
            saa = as.vector(saa), 
            vza = as.vector(vza), 
            cld_frac = as.vector(cld_frac), 
            cld_pres = as.vector(cld_pres), 
            sfc_alb = as.vector(sfc_alb), 
            amf_tot = as.vector(amf_tot), 
            tcwv = as.vector(tcwv), 
            tcwv_uncert = as.vector(tcwv_uncert)
        ) %>% filter(as.vector(inTF[ix, iy]), !is.na(tcwv))

        # H2O column in mol m-2, dry-air column in mol m-2, XH2O in ppm
        var_df = var_df %>% 
                 mutate(h2o_vcd = tcwv / Mh2o, 
                        h2o_vcd_uncert = tcwv_uncert / Mh2o, 
                        air_vcd = psfc * 100 / g / Mair - h2o_vcd, 
                        xh2o = h2o_vcd / air_vcd * 1E6, 
                        xh2o_uncert = h2o_vcd_uncert / air_vcd * 1E6)

        # ----------------------- corner data frame ----------------------------
        # [corner, ground_pixel, scanline], corner runs fastest
        loc_df = data.frame(
            corner = rep(corner, times = nx * ny), 
            indx_across_track = rep(rep(indx_across_track[ix], each = 4), times = ny), 
            indx_along_track  = rep(indx_along_track[iy], each = 4 * nx), 
            lons = as.vector(lons), 
            lats = as.vector(lats)
        )

        merge_df = loc_df %>% 
                   inner_join(var_df, by = c('indx_across_track', 'indx_along_track'))
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
