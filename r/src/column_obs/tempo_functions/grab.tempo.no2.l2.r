
# Read TEMPO NO2 L2 granules, mirroring grab.tempo.no2.l3():
#   - only the sub-domain is read from disk (start/count on [xtrack, mirror_step])
#   - same output columns as L3 (qa, hsfc, pblh, psfc, ptrp, sza, saa, *_vcd,
#     amf_*, ecf, cld_pres, sec, datestr, time_utc) 
#   - plus L2-only fields: pixel corners (corner, lons, lats), xtrack, 
#     mirror_step, relative azimuth angle (raa) and wind speed (wspd)
# NOTE: averaging kernels are NOT returned here; X-STILT applies the AK per 
#       particle at footprint-weighting time (see get.wgt.tropomi.func for the
#       TROPOMI equivalent), using the scattering_weights / amf_troposphere /
#       Eta pressure fields read directly from the L2 file at the receptor pixel
# DW, 11/09/2026

if (F) {
    tempo.path = sensor_path = obs_path
    tempo.fn = sensor_fn = obs_fn
}

grab.tempo.no2.l2 = function(tempo.path = NULL, timestr = NULL, 
                             lon_lat = NULL, tempo.fn = NULL) {
    
    if (is.null(tempo.fn)) 
        stop('Missing TEMPO L2 NO2 file for your region...\n')

    out_df = NULL
    for (fn in tempo.fn) {

        dat = nc_open(fn)
        xtr = ncvar_get(dat, 'xtrack')
        mrr = ncvar_get(dat, 'mirror_step')
        cor = dat$dim$corner$vals
        lat = ncvar_get(dat, 'geolocation/latitude')     # [xtrack, mirror_step]
        lon = ncvar_get(dat, 'geolocation/longitude')

        # find the xtrack / mirror_step index block whose pixel centers fall 
        # inside lon_lat; skip the granule if nothing overlaps
        inTF = !is.na(lon) & !is.na(lat) & 
               lon >= lon_lat$minlon & lon <= lon_lat$maxlon & 
               lat >= lon_lat$minlat & lat <= lon_lat$maxlat
        if (!any(inTF)) { nc_close(dat); next }
        ix = which(apply(inTF, 1, any)); ix = min(ix):max(ix)
        iy = which(apply(inTF, 2, any)); iy = min(iy):max(iy)
        nx = length(ix)
        ny = length(iy)

        # readers for 2D [xtrack, mirror_step] and 3D [level, xtrack, mirror_step]
        get  = function(v) ncvar_get(dat, v, start = c(min(ix), min(iy)), 
                                             count = c(nx, ny))
        get3 = function(v) ncvar_get(dat, v, start = c(1, min(ix), min(iy)), 
                                             count = c(-1, nx, ny))

        # troposphere nitrogen dioxide vertical column, molecules/cm^2
        tno2_vcd     = get('product/vertical_column_troposphere')
        tno2_vcd_err = get('product/vertical_column_troposphere_uncertainty')

        # main data quality flag, 0 normal, 1 suspicious, 2 bad
        qa = get('product/main_data_quality_flag')

        # effective cloud fraction from cloud retrieval
        ecf  = get('support_data/eff_cloud_fraction')
        psfc = get('support_data/surface_pressure')       # hPa
        ptrp = get('support_data/tropopause_pressure')    # hPa
        hsfc = get('support_data/terrain_height')         # in m
        pblh = get('support_data/pbl_height')             # in m
        sza  = get('geolocation/solar_zenith_angle')
        saa  = get('geolocation/solar_azimuth_angle') # needed for slant columns
        raa  = get('geolocation/relative_azimuth_angle')

        # total and stratospheric NO2 columns, molecules/cm^2
        xno2_vcd = get('support_data/vertical_column_total')
        sno2_vcd = get('product/vertical_column_stratosphere')

        # air mass factors
        amf_tot  = get('support_data/amf_total')
        amf_trp  = get('support_data/amf_troposphere')
        amf_strt = get('support_data/amf_stratosphere')

        # cloud pressure for AMF computation, hPa
        cld_pres = get('support_data/amf_cloud_pressure')

        # wind speed used in the retrieval, m/s
        wspd = get('support_data/wind_speed')

        # pixel corners (SW, SE, NE, NW), [corner, xtrack, mirror_step]
        lons = get3('geolocation/longitude_bounds')
        lats = get3('geolocation/latitude_bounds')

        # per mirror_step time, sec since 1980-01-06T00:00:00Z
        sec = ncvar_get(dat, 'geolocation/time', start = min(iy), count = ny)

        nc_close(dat)

        # ------------------------- pixel-level data frame ---------------------
        # matrices are [xtrack, mirror_step], so as.vector() runs xtrack fastest
        var_df = data.frame(
            xtrack = rep(xtr[ix], times = ny), 
            mirror_step = rep(mrr[iy], each = nx), 
            lon = as.vector(lon[ix, iy]), 
            lat = as.vector(lat[ix, iy]), 
            qa = as.vector(qa),
            hsfc = as.vector(hsfc),
            pblh = as.vector(pblh),
            psfc = as.vector(psfc),
            ptrp = as.vector(ptrp),
            sza = as.vector(sza),
            saa = as.vector(saa),
            raa = as.vector(raa),
            tno2_vcd = as.vector(tno2_vcd),
            tno2_vcd_err = as.vector(tno2_vcd_err),
            sno2_vcd = as.vector(sno2_vcd),
            xno2_vcd = as.vector(xno2_vcd),
            amf_tot = as.vector(amf_tot),
            amf_trp = as.vector(amf_trp),
            amf_strt = as.vector(amf_strt),
            ecf = as.vector(ecf),
            cld_pres = as.vector(cld_pres), 
            wspd = as.vector(wspd), 
            sec = rep(sec, each = nx)
        ) %>% filter(as.vector(inTF[ix, iy]))    # keep centers inside lon_lat

        pix_date = as.POSIXct(var_df$sec, 'GMT', 
                              origin = '1980-01-06T00:00:00Z')
        
        var_df = var_df %>% 
                 mutate(datestr = as.POSIXlt(pix_date, 'UTC'), 
                        time_utc = as.numeric(format(pix_date, '%Y%m%d%H%M%S')))

        # ------------------------- corner data frame --------------------------
        # [corner, xtrack, mirror_step], corner runs fastest
        loc_df = data.frame(
            corner = rep(cor, times = nx * ny), 
            xtrack = rep(rep(xtr[ix], each = length(cor)), times = ny), 
            mirror_step = rep(mrr[iy], each = length(cor) * nx), 
            lons = as.vector(lons), 
            lats = as.vector(lats)
        )

        # ---------------------- merge ---------------------------------------
        # 4 rows (corners) per pixel, like the TROPOMI readers
        gran_df = loc_df %>% inner_join(var_df, by = c('xtrack', 'mirror_step'))
        out_df = rbind(out_df, gran_df)
    }   # end for fn

    closeAllConnections()
    return(out_df)
}
