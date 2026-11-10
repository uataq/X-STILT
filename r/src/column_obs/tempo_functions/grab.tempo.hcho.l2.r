
grab.tempo.hcho.l2 = function(tempo.path = NULL, timestr = NULL, lon_lat = NULL,
                              tempo.fn = NULL, ak_pres = c(NA, 900)[1]) {
    
    if (is.null(tempo.fn)) {
        
    } 

    dat = nc_open(tempo.fn)
    mrr = ncvar_get(dat, 'mirror_step')
    xtr = ncvar_get(dat, 'xtrack')
    cor = dat$dim$corner$vals
    swt = dat$dim$swt_level$vals
    lat = ncvar_get(dat, 'geolocation/latitude')            # pixel center
    lats = ncvar_get(dat, 'geolocation/latitude_bounds')    # (SW,SE,NE,NW)
    lon = ncvar_get(dat, 'geolocation/longitude')
    lons = ncvar_get(dat, 'geolocation/longitude_bounds')
    
    # ------------------------- corner matrix -----------------------------
    # lon/lat for the polygon corners
    dimnames(lats) = dimnames(lons) = list(cor, xtr, mrr)
    loc_name = list(c('lats', 'lons')); loc_list = list(lats, lons)
    loc_array = array(
        data = do.call(cbind, loc_list), 
        dim = c(dim(loc_list[[1]]), length(loc_list)), 
        dimnames = c(dimnames(loc_list[[1]]), loc_name)
    )
    loc_df = dcast(melt(loc_array), Var1 +Var2 +Var3~Var4, value.var = 'value')
    colnames(loc_df)[1:3] = c('corner', 'xtrack', 'mirror_step')

    # figure out the xtrack and mirror_step that fit lon/lat
    if (!is.null(lon_lat)) {    
        sel_loc_df = loc_df %>% 
                     filter(lons >= lon_lat$minlon - 0.01, 
                            lons <= lon_lat$maxlon + 0.01,
                            lats >= lon_lat$minlat - 0.01, 
                            lats <= lon_lat$maxlat + 0.01) 
        if (nrow(sel_loc_df) == 0) return()

        sel_xtr_indx = which(xtr %in% unique(sel_loc_df$xtrack))
        sel_mrr_indx = which(mrr %in% unique(sel_loc_df$mirror_step))
        sel_loc_df = loc_df %>% 
                     filter(xtrack %in% unique(sel_loc_df$xtrack), 
                            mirror_step %in% unique(sel_loc_df$mirror_step))
    } else { 
        sel_xtr_indx = 1:length(xtr)
        sel_mrr_indx = 1:length(mrr) 
        sel_loc_df = loc_df
    }
    
    sel_lon  = lon[sel_xtr_indx, sel_mrr_indx]
    sel_lat  = lat[sel_xtr_indx, sel_mrr_indx]
    #sel_lons = lons[, sel_xtr_indx, sel_mrr_indx]
    #sel_lats = lats[, sel_xtr_indx, sel_mrr_indx]

    # ------------------------- loading data 
    # 0 normal, 1 suspicious, 2 bad
    qa = ncvar_get(dat, 'product/main_data_quality_flag')[sel_xtr_indx, sel_mrr_indx]

    sza = ncvar_get(dat, 'geolocation/solar_zenith_angle')[sel_xtr_indx, sel_mrr_indx]
    
    saa = ncvar_get(dat, 'geolocation/solar_azimuth_angle')[sel_xtr_indx, sel_mrr_indx]
    
    hsfc = ncvar_get(dat, 'support_data/terrain_height')[sel_xtr_indx, sel_mrr_indx]  # in m
    #qf = ncvar_get(dat, 'support_data/ground_pixel_quality_flag')

    # surface and tropopause pressure, hPa, [xtrack,mirror_step]
    psfc = ncvar_get(dat, 'support_data/surface_pressure')[sel_xtr_indx, sel_mrr_indx]

    albedo = ncvar_get(dat, 'support_data/albedo')[sel_xtr_indx, sel_mrr_indx]

    # [xtrack,mirror_step], molecules/cm^2
    hcho_vcd = ncvar_get(dat, 'product/vertical_column')[sel_xtr_indx, sel_mrr_indx]

    hcho_vcd_err = ncvar_get(dat, 'product/vertical_column_uncertainty')[sel_xtr_indx, sel_mrr_indx]

    # xtrack dependent HCHO slant column background correction, molecules/cm^2
    hcho_bg_corr = ncvar_get(dat, 'support_data/background_correction')[sel_xtr_indx, sel_mrr_indx]

    # total formaldehyde air mass factor (AMF) calculated from surface to top of atmosphere
    amf_tot = ncvar_get(dat, 'support_data/amf')[sel_xtr_indx, sel_mrr_indx]

    # -------- profiles =
    # vertical profile of formaldehyde partial column, molecules/cm^2
    hcho_prof = ncvar_get(dat, 'support_data/gas_profile')[, sel_xtr_indx, sel_mrr_indx]

    # scattering weights [scat levels, xtrack, mirror_step]
    sct_wgt = ncvar_get(dat, 'support_data/scattering_weights')[, sel_xtr_indx, sel_mrr_indx]

    temp_pri = ncvar_get(dat, 'support_data/temperature_profile')[, sel_xtr_indx, sel_mrr_indx]

    # [mirror_step], sec since1980-01-06T00:00:00Z
    sec = ncvar_get(dat, 'geolocation/time')[sel_mrr_indx] 

    # assign proper dimensions to variables
    time_df = data.frame(sec, mirror_step = mrr[sel_mrr_indx]) %>% 
              mutate(date = as.POSIXct(sec, 'GMT', origin = '1980-01-06T00:00:00Z'))



    # assigning dimensions --------
    dimnames(sct_wgt) = dimnames(hcho_prof) = dimnames(temp_pri) = 
        list(swt, xtr[sel_xtr_indx], mrr[sel_mrr_indx])

    dimnames(sel_lat) = dimnames(sel_lon) = dimnames(qa) = dimnames(hsfc) = 
    dimnames(psfc) = dimnames(albedo) = dimnames(hcho_vcd) = 
    dimnames(hcho_vcd_err) = dimnames(hcho_bg_corr) = dimnames(amf_tot) = 
        list(xtr[sel_xtr_indx], mrr[sel_mrr_indx])

    
    # also retrieve Eta values for pressure levels [73 levels] --------
    eta_vals = ncatt_get(dat, 'support_data/surface_pressure')
    plev = array(NA, dim = c(length(eta_vals$Eta_A), dim(psfc)), 
                     dimnames = list(1:length(eta_vals$Eta_A), 
                                     rownames(psfc), colnames(psfc)))

    for ( ll in 1:length(eta_vals$Eta_A) )  # from surface to up 
        plev[ll, , ] = eta_vals$Eta_A[ll] + eta_vals$Eta_B[ll] * psfc
    # plev_df = melt(plev)
    # colnames(plev_df) = c('level', 'xtrack', 'mirror_step', 'pres')


    # -------------------------  merge matrix 1 ---------------------------
    var_name = list(c('lat', 'lon', 'qa', 'hsfc', 'psfc', 'albedo', 'hcho_vcd', 
                      'hcho_vcd_err', 'hcho_bg_corr', 'amf_tot'))
    var_list = list(sel_lat, sel_lon, qa, hsfc, psfc, albedo, hcho_vcd, 
                    hcho_vcd_err, hcho_bg_corr, amf_tot)
    var_array = array(  
        data = do.call(cbind, var_list), 
        dim = c(dim(var_list[[1]]), length(var_list)), 
        dimnames = c(dimnames(var_list[[1]]), var_name)
    )   # assuming all matrices in the list have equal dimensions

    var_df = dcast(melt(var_array), Var1 + Var2~Var3, value.var = 'value')
    colnames(var_df)[1:2] = c('xtrack', 'mirror_step')


    # ------------------------- merge matrix 3 -----------------------------
    # prior and scattering profiles
    prf_name = list(c('sct_wgt', 'hcho_prof', 'temp_pri', 'pres_lower'))
    prf_list = list(sct_wgt, hcho_prof, temp_pri, plev[-73, , ])
    prf_array = array(
        data = do.call(cbind, prf_list), 
        dim = c(dim(prf_list[[1]]), length(prf_list)), 
        dimnames = c(dimnames(prf_list[[1]]), prf_name)
    )
    prf_df = dcast(melt(prf_array), Var1 +Var2 +Var3~Var4, value.var = 'value')
    colnames(prf_df)[1:3] = c('layer', 'xtrack', 'mirror_step')


    # cropping ------
    sel_df = var_df %>% left_join(time_df, by = 'mirror_step') %>%
             right_join(sel_loc_df, by = c('xtrack', 'mirror_step')) 

    # melt scattering weights and calculate AK using AMF ----------------
    amf_df = unique(sel_df %>% dplyr::select('qa', contains('amf'), 
                                             'xtrack', 'mirror_step'))
    
    # The tropospheric averaging kernel can be obtained by scaling the kernel by
    # M/Mtrop (see [RD20]) and setting all elements of the kernel to zero above 
    # the tropopause layer
    wgt_df = prf_df %>% right_join(amf_df, by = c('xtrack', 'mirror_step')) %>%
             na.omit() %>% mutate(ak_tot = sct_wgt / amf_tot)
    

    # only pick the bottom few levels for surface AK 
    if (is.na(ak_pres)) {
        wgt_sfc_df = wgt_df %>% group_by(xtrack, mirror_step) %>% 
                     filter(pres_lower == max(pres_lower))
    } else {
        wgt_sfc_df = wgt_df %>% group_by(xtrack, mirror_step) %>% 
                     filter(pres_lower > ak_pres)
    }

    # take the vertical average of trop AK, scattering weights, and prior NO2
    wgt_avg_df = wgt_sfc_df %>% 
                 summarise(ak_tot = mean(ak_tot), sct_wgt = mean(sct_wgt), 
                           hcho_vcd_prior = mean(hcho_prof), 
                           pres = mean(pres_lower), .groups = 'drop') %>% ungroup() 

    out_df = sel_df %>% left_join(wgt_avg_df, by = c('xtrack', 'mirror_step'))


    if (F) {

        avg_p_df = data.frame(layer = 1:72, 
                        p = apply(plev, 1, function(x) mean(x, na.rm = T))[-73])
        avg_wgt_df = wgt_df %>% na.omit() %>% group_by(layer, qa) %>%
                     summarise_all(mean) %>% ungroup() %>%
                     left_join(avg_p_df, by = 'layer') %>% 
                     mutate(qf = ifelse(qa == 0, '0 - normal', 
                                 ifelse(qa == 1, '1 - suspicious', '2 - bad')))

        p0 = ggplot(data = avg_wgt_df, aes(color = as.factor(qf))) + 
             theme_bw() + ylim(1000, 0) + labs(y = NULL)
        t1 = p0 + geom_point(aes(temp_pri - 273.15, p))
        n1 = p0 + geom_point(aes(hcho_prof, p)) + 
                  xlab('HCHO profile\n[molecules/cm^2]')
        w1 = p0 + geom_point(aes(sct_wgt, p)) 
        m1 = p0 + geom_point(aes(amf_tot, p))
        a1 = p0 + geom_vline(xintercept = 1, linetype = 2) + 
                  geom_hline(yintercept = 850, linetype = 2) +
                  geom_point(aes(ak_tot, p))
        pp = ggarrange(t1, n1, w1, m1, a1, nrow = 1, ncol = 5,common.legend = T)
        ggsave(pp, filename = 'wgt_hcho_profile_eg.png', width = 14, height = 6,
                   bg = 'white')

    }

    nc_close(dat)
    closeAllConnections()
    return(out_df)
}