
# if ak_pres is NA, default is to pick the very bottom pressure layer for wgts
if (F) {
    tempo.path = sensor_path
    tempo.fn = sensor_fn 
}

grab.tempo.no2.l3 = function(tempo.path = NULL, timestr = NULL, 
                             lon_lat = NULL, tempo.fn = NULL){
    
    if (is.null(tempo.fn)) stop('Missing TEMPO L3 NO2 file for your region...\n')

    out_df = NULL
    for (fn in tempo.fn) {
        
        dat = nc_open(fn)
        lat = ncvar_get(dat, 'latitude')            # pixel center
        lon = ncvar_get(dat, 'longitude')

        # find lon/lat indices, xtrack and mirror_step that fit lon/lat
        ix = which(lon >= lon_lat$minlon & lon <= lon_lat$maxlon)
        iy = which(lat >= lon_lat$minlat & lat <= lon_lat$maxlat)
        st = c(min(ix), min(iy), 1)
        ct = c(length(ix), length(iy), 1)
        get = function(v) ncvar_get(dat, v, start = st, count = ct)

        # troposphere nitrogen dioxide vertical column, molecules/cm^2
        tno2_vcd     = get('product/vertical_column_troposphere')
        tno2_vcd_err = get('product/vertical_column_troposphere_uncertainty')

        # main data quality flag, 0 normal, 1 suspicious, 2 bad
        qa = get('product/main_data_quality_flag')

        # effective cloud fraction from cloud retrieval
        ecf  = get('support_data/eff_cloud_fraction')
        psfc = get('support_data/surface_pressure')
        ptrp = get('support_data/tropopause_pressure')
        hsfc = get('support_data/terrain_height')      # in m
        pblh = get('support_data/pbl_height') # in m
        sza = get('geolocation/solar_zenith_angle')
        saa = get('geolocation/relative_azimuth_angle')

        # total and stratospheric NO2 columns, molecules/cm^2
        xno2_vcd = get('support_data/vertical_column_total')
        sno2_vcd = get('product/vertical_column_stratosphere')

        # air mass factors
        amf_tot  = get('support_data/amf_total')
        amf_trp  = get('support_data/amf_troposphere')
        amf_strt = get('support_data/amf_stratosphere')

        # cloud pressure for AMF computation, hPa
        cld_pres = get('support_data/amf_cloud_pressure')
        sel_lon = lon[ix]
        sel_lat = lat[iy]
        sec = ncvar_get(dat, 'time')

        # scan start time and end time (departure by 1 hour for entire scan)
        sec0 = ncatt_get(dat, 0, 'time_coverage_start_since_epoch')$value
        sec1 = ncatt_get(dat, 0, 'time_coverage_end_since_epoch')$value

        # east/west edge of observed columns; weight is [lon, lat] in R
        w = ncvar_get(dat, 'weight')
        obs_lon = lon[rowSums(w > 0, na.rm = TRUE) > 0]
        lon_e = max(obs_lon)
        lon_w = min(obs_lon)

        var_df = data.frame(
            lon = rep(sel_lon, each = length(sel_lat)),
            lat = rep(sel_lat, times = length(sel_lon)),
            qa = as.vector(qa),
            hsfc = as.vector(hsfc),
            pblh = as.vector(pblh),
            psfc = as.vector(psfc),
            ptrp = as.vector(ptrp),
            sza = as.vector(sza),
            saa = as.vector(saa),
            tno2_vcd = as.vector(tno2_vcd),
            tno2_vcd_err = as.vector(tno2_vcd_err),
            sno2_vcd = as.vector(sno2_vcd),
            xno2_vcd = as.vector(xno2_vcd),
            amf_tot = as.vector(amf_tot),
            amf_trp = as.vector(amf_trp),
            amf_strt = as.vector(amf_strt),
            ecf = as.vector(ecf),
            cld_pres = as.vector(cld_pres)
        )

        # per-pixel seconds since 1980-01-06, east pixels earliest
        # these pixel seconds are estimated for specific lon_lat region
        pix_sec = sec0 + (lon_e - var_df$lon) / (lon_e - lon_w) * (sec1 - sec0)
        pix_date = as.POSIXct(pix_sec, 'GMT', origin = '1980-01-06T00:00:00Z')
        
        var_df = var_df %>% mutate(sec = pix_sec, 
                                   datestr = as.POSIXlt(pix_date, 'UTC'),
                                   time_utc = as.numeric(format(pix_date, '%Y%m%d%H%M%S')))
        
        out_df = rbind(out_df, var_df)
        nc_close(dat)
    }

    closeAllConnections()
    return(out_df)
}