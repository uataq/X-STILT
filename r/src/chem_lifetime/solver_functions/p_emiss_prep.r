
# e.g., epa_tz = 'MST', insert the local standard time for your power plants
# add a module for emission uncertainty, DW, 04/15/2022
# also perturb CO2 and CO emissions if perturbTF = T, DW, 08/11/2022 

p_emiss_prep = function(site, output, p, timestr, aq_invent, ghg_invent, eno_fn,
                        eco_fn, ech4_fn, eco2_fn, epa_lst = NA, eno_sf = 1, 
                        eco_sf = 1, ech4_sf = 1, eco2_sf = 1, tno2_fn, 
                        xco_fn = NA, xch4_fn = NA, xco2_fn = NA, perturbTF = F,
                        perturb_indx = NA, perturb_fn = NA) {
    
    library(dplyr)
    receptor = output$receptor
    crs0 = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'

    # if emission domain is not 
    # if ( is.na(xmin) | is.na(xmax) | is.na(ymin) | is.na(ymax) ) {
    #     xmin = min(p$long); xmax = max(p$long)
    #     ymin = min(p$lati); ymax = max(p$lati)
        
    #     if ( ghg_invent == 'vulcan' | aq_invent == 'vulcan' ) {
    #         xmin = min(p$long) - 1; xmax = max(p$long) + 1
    #         ymin = min(p$lati) - 1; ymax = max(p$lati) + 1
    #     }  
    # } 

    # ----------------------------------------------------------------------
    # loading pre-processed emissions, see function preproc_emiss()
    eno_rt = raster(eno_fn)
    eco_rt = raster(eco_fn)
    ech4_rt = raster(ech4_fn)
    eco2_rt = raster(eco2_fn)
    
    # ----------------------------------------------------------------------
    # for lifetime look up relationship, need to know SZA -----------------
    # grab initial footprint, not the weighted footprint - need to re-weight
    cat('p_emiss_prep(): subset particles and add SZA...\n') 
    if ('foot_before_weight' %in% colnames(p)) 
        p = p %>% dplyr::select(-foot) %>% rename(foot = foot_before_weight)

    # 0.5. select trajectories for speeding up, only 12 hours back is sufficient
    library(GeoLight)
    pex = p %>% filter(!is.na(recp_lat)) %>% 
          dplyr::select(indx, time, date, lati, long, zagl, zsfc, foot, mlht, 
                        pres, temz, dens, xhgt, any_of(c('sphu', 'rhfr'))) %>% 
          mutate(psza = zenith(solar(date), long, lati)) %>% arrange(time)


    # assign emission to paticle locations using raster::extract --------------
    cat('p_emiss_prep(): assigning emissions to particles...\n')
    p_sp = pex[, c('long', 'lati')]
    coordinates(p_sp) = ~long + lati; crs(p_sp) = crs0
    suppressWarnings({      
        # suppress the following raster warning, when emission projection does not match the particle projection, e.g., for dealing with Vulcan
        pex$eno = raster::extract(eno_rt, p_sp, method = 'simple')
        pex$eco = raster::extract(eco_rt, p_sp, method = 'simple') 
        pex$ech4 = raster::extract(ech4_rt, p_sp, method = 'simple') 
        pex$eco2 = raster::extract(eco2_rt, p_sp, method = 'simple') 
    })

    pex = pex %>% mutate(eno = ifelse(is.na(eno), 0, eno), 
                         eco = ifelse(is.na(eco), 0, eco), 
                         ech4 = ifelse(is.na(ech4), 0, ech4), 
                         eco2 = ifelse(is.na(eco2), 0, eco2))
    cat(paste('p_emiss_prep(): max ENO of', signif(max(pex$eno), 3), 
              'umol-ENO m-2 s-1...\n'))


    # ----------------------------------------------------------------------
    # only scale the max value for PP based on EPA, DW, 09/09/2021
    # correct for bottom-up prior emissions if needed (for optimization)
    if ( aq_invent == 'epa' ) {    # only scale emiss for PP based on EPA

        epa_co2_df = epa_lst$epa_co2_df
        epa_no_df = epa_lst$epa_no_df 
        if (is.null(epa_co2_df) | is.null(epa_no_df)) 
            stop('p_emiss_prep(): PP emiss or EPA file is missing when trying to scale emissions based on EPA...\n')

        epa_df = full_join(epa_co2_df %>% rename(epa_co2_pp = epa_pp), 
                           epa_no_df %>% rename(epa_no_pp = epa_pp), 
                           by = c('timestr_lst', 'timestr_utc'))
        
        # re-locate the emission grid where PP falls into 
        # x, y coordinate in lower left corner now, DW, 03/22/2022
        names(eno_rt) = 'enox'
        eno_res = res(eno_rt)[1]
        eno_xy = as.data.frame(eno_rt, xy = T) %>% 
                 filter(enox == epa_lst$epa_pp) %>%
                 mutate(xmn = x - eno_res / 2, ymn = y - eno_res / 2, 
                        xmx = x + eno_res / 2, ymx = y + eno_res / 2)
        
        # correct emission if particles may pass over PP
        if (nrow(eno_xy) == 1) {
            pex = pex %>% 
                  mutate(timestr_p = format(date, format = '%Y%m%d%H')) %>% 
                  left_join(epa_df, by = c('timestr_p' = 'timestr_utc')) %>% 

                  # if particles fall within the PP grids
                  # if so, scale emissions by sf
                  mutate(ppTF = long >= eno_xy$xmn & long < eno_xy$xmx &
                                lati >= eno_xy$ymn & lati < eno_xy$ymx, 
                         sf_pp_no  = ifelse(ppTF, sf_pp_no, 1), 
                         sf_pp_co2 = ifelse(ppTF, sf_pp_co2, 1), 
                         eno = eno * sf_pp_no, eco2 = eco2 * sf_pp_co2) %>% 
                  dplyr::select(-c(timestr_p,epa_co2_pp,epa_no_pp, timestr_lst))

            cat(paste('p_emiss_prep(): scaling ENOx by', 
                      signif(min(pex$sf_pp_no, na.rm = T), 3), '-', 
                      signif(max(pex$sf_pp_no, na.rm = T), 3), '...\n\n'))

        } else cat('p_emiss_prep(): No trajec interact with the PP emission...no need to scale ENOx...\n\n') 
        # if e*_xy not found, it means that particles will not interact with fluxes from PP, thus, no need to adjust emissions

    }  # end if scaling PP emissions based on EPA PP
    

    # ----------------------------------------------------------------------
    # perturb the emission if for local ensemble KF, DW, 04/15/2022
    # find particle that fall within the emission perturbation grid
    if ( perturbTF ) {
        cat('p_emiss_prep(): obtaining perturbation scaling factor...\n')
        
        # load perturbation scaling factor, DW, 04/15/2022
        sf_df = readRDS(perturb_fn) %>% filter(ens == perturb_indx)
        sf_rt = rasterFromXYZ(sf_df[, c('x', 'y', 'sf')])
        pex$sf_perturb = raster::extract(sf_rt, p_sp, method = 'simple') 
        
        # if NA meaning no need to perturb the NOx emissions
        # also perturb CO2 and CO emissions if perturbTF = T, DW, 08/11/2022 
        pex = pex %>% mutate(sf_perturb = ifelse(is.na(sf_perturb), 
                                                 1, sf_perturb), 
                             #eno_init = eno,    # preserve original emission
                             eno = eno * sf_perturb, 
                             eco = eco * sf_perturb, 
                             ech4 = ech4 * sf_perturb,
                             eco2 = eco2 * sf_perturb) 
    }   # end if


    ### -----------------------------------------------------------------------
    # use initial footprint, AK, and a normalized PWF, DW, 08/25/2021
    # default PWF is weighted based on the entire column, modify it
    g = 9.8        	         # m s-2
    Mair = 29 / 1E3        	 # kg mol-1 

    # tropospheric column according to TROPOMI tropopause
    if (!file.exists(tno2_fn)) 
        stop('p_emiss_prep(): TROPOMI NO2 file not found at the moment...\nplease double check')
    no2_prof = get.wgt.tropomi.func(output, tno2_fn, 'NO2')
    if ( is.null(no2_prof) ) return()
    
    no2_df = no2_prof$combine.prof %>% filter(stiltTF)
    no2_info = no2_prof$trp.info

    # norm PWFs have an average val of 1 from sfc to the top stilt lev
    # meaning that NOx from lower levels will be weighted more
    # update AK.PWF for NOx calculation
    wgt_df = no2_df %>% 
             mutate(pwf_norm = dp / mean(dp)) %>% rename(ak_tno2 = ak.norm)
    
    # calculate dry-air column density for different columns 
    # STILT column, TROPOMI tropospheric column, TROPOMI total column 
    no2_info$stilt_psfc = receptor$psfc
    no2_info$air_vcd_stilt = sum((1 - wgt_df$q) / g / Mair * wgt_df$dp * 100)
    no2_info$h2ov_vcd_stilt = sum(wgt_df$q / g / Mair * wgt_df$dp * 100)
    
    no2_info$air_vcd_tropo = unique(no2_df$xdry.tot)
    no2_info$tno2_ppb = no2_info$no2_vcd_tropo / no2_info$air_vcd_tropo * 1e9


    ### -----------------------------------------------------------------------
    if (!is.na(xco_fn)) {      # also get total column dry air density 
        if (!file.exists(xco_fn)) 
            stop('p_emiss_prep(): TROPOMI CO file not found at the moment...\nplease double check or re-run the code')
        
        co_prof = get.wgt.tropomi.func(output, xco_fn, 'CO', tropo.no2TF = F)

        if ( !is.null(co_prof)) {
            co_df = co_prof$combine.prof 
            co_info = co_prof$trp.info 
            co_info$air_vcd_total = unique(co_df$xdry.tot)
            co_info$xco_ppb = co_info$co_vcd / co_info$air_vcd_total * 1e9

            # interpolate AK of XCO to TROPOMI NO2 levels
            wgt_df = wgt_df %>% 
                     mutate(ak_xco = approx(co_df$lower_pres, co_df$ak.norm, 
                                            xout = wgt_df$lower_pres)$y, 
                            ak_xco = ifelse(lower_pres > max(co_df$lower_pres), 
                                            co_df$ak.norm[1], ak_xco))
                        
        } else co_info = NULL
    } else co_info = NULL # end if


    ### -----------------------------------------------------------------------
    if (!is.na(xch4_fn)) {      # also get total column dry air density 
        if (!file.exists(xch4_fn)) 
            stop('p_emiss_prep(): TROPOMI CO file not found at the moment...\nplease double check or re-run the code')
        
        ch4_prof = get.wgt.tropomi.func(output, xch4_fn, 'CH4', tropo.no2TF = F)
        
        if ( !is.null(ch4_prof) ) {
            ch4_df = ch4_prof$combine.prof 
            ch4_info = ch4_prof$trp.info 
            ch4_info$air_vcd_total = unique(ch4_df$xdry.tot)
            ch4_info$xch4_ppb = ch4_info$xch4_bc    # already in PPB

            # interpolate AK of XCO to TROPOMI NO2 levels
            wgt_df = wgt_df %>% 
                     mutate(ak_xch4 = approx(ch4_df$lower_pres, ch4_df$ak.norm, 
                                             xout = wgt_df$lower_pres)$y, 
                            ak_xch4 = ifelse(lower_pres >max(ch4_df$lower_pres),
                                             ch4_df$ak.norm[1], ak_xch4))
                        
        } else ch4_info = NULL
    } else ch4_info = NULL # end if


    ### -----------------------------------------------------------------------
    if (!is.na(xco2_fn)) {   # get xCO2 ak and aprior profiles from OCO

        if (!file.exists(xco2_fn)) 
            stop('p_emiss_prep(): OCO XCO2 file not found at the moment...\nplease double check or re-run the code')
        
        co2_info = get.oco.info(receptor = receptor, oco.fn = xco2_fn, 
                                diff.td = 1e3) 
        
        if ( !is.null(co2_info)) {

            # calculate dry air column density of OCO columns
            co2_info$air_vcd_total = co2_info$oco.psfc * 100 / g / Mair - co2_info$oco.xh2o
            
            wgt_df = wgt_df %>% 
                     mutate(ak_xco2 = approx(co2_info$pres, co2_info$ak.norm, 
                                            xout = wgt_df$lower_pres)$y,
                            ak_xco2 = ifelse(lower_pres > max(co2_info$pres), 
                                            co2_info$ak[1], ak_xco2), 
                            ap_xco2 = approx(co2_info$pres, co2_info$ap, 
                                            xout = wgt_df$lower_pres)$y)
        } else co2_info = NULL
    } else co2_info = NULL


    ### -----------------------------------------------------------------------
    # merge weighting factor into particle and recalculate the trajec-level foot
    pex = pex %>% 
          left_join(wgt_df %>% dplyr::select(indx, pwf_norm, contains('ak_')), 
                    by = 'indx') %>% 
          
          # pressure weighting is applied to STILT levels, not the entire atmos
          # here norm.pwf = dp_model / sum(dp_model), not the PWF from satellite
          mutate(foot_wgt = foot * pwf_norm, 
                 demiss_no = eno * foot_wgt * eno_sf, 
                 demiss_co = eco * foot_wgt * eco_sf, 
                 demiss_ch4 = ech4 * foot_wgt * ech4_sf, 
                 demiss_co2 = eco2 * foot_wgt * eco2_sf) %>% 
          dplyr::select(-c('pwf_norm'))

    # check delt - default is 1 min from trajec file
    times = unique(pex$time)
    delt = min(abs(unique(diff(times))))  # in min
    # if (max(pex$time) != -delt) 
    #     stop(paste0('p_emiss_prep(): missing near-field trajec e.g., with time of -', delt, ' min.\nPLEASE check your xmin, xmax, ymin, ymax...DO NOT subset near-field trajec!!!\n'))

    return(list(pex = pex, no2_info = no2_info, co_info = co_info, 
                           ch4_info = ch4_info, co2_info = co2_info))
}

