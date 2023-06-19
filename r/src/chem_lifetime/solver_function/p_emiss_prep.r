
# inner functions
if (F) { xmin = -113; xmax = -110; ymin = 32; ymax = 34 }

# e.g., epa_tz = 'MST', insert the local standard time for your power plants
# add a module for emission uncertainty, DW, 04/15/2022
# also perturb CO2 and CO emissions if perturbTF = T, DW, 08/11/2022 

p_emiss_prep = function(site, output, p, timestr, eno_fn, eco_fn, eco2_edgar_fn,
                        eco2_odiac_fn = NA, epa_fn = NA, epa_name = NA,
                        epa_tz = NA, no2_fn, co_fn = NA, co2_fn = NA, eno_sf, 
                        eco_sf, eco2_sf, perturbTF = F, perturb_indx = 1, 
                        perturb_fn = NA, aq_invent, ghg_invent, xmin = NA, 
                        xmax = NA, ymin = NA, ymax = NA, nmins = -720) {
    
    library(dplyr)
    receptor = output$receptor
    
    # if emission domain is not 
    if (is.na(xmin) | is.na(xmax) | is.na(ymin) | is.na(ymax)) {
        xmin = min(p$long); xmax = max(p$long)
        ymin = min(p$lati); ymax = max(p$lati)
    }

    # load EDGAR CO2 emissions first --------------------------------------
    cat(paste('p_emiss_prep(): loading ECO2 from', toupper(ghg_invent),'...\n'))
    eco2_lst = load_eco2(timestr, invent = ghg_invent, eco2_edgar_fn, 
                         eco2_odiac_fn, epa_fn, epa_name, epa_tz, xmin, xmax, 
                         ymin, ymax)
    eco2_rt = eco2_lst$edgar_rt     # EDGAR eCO2


    # load EDGAR NO emissions ----------------------------------------------
    cat(paste('p_emiss_prep(): loading ENOx from', toupper(aq_invent), '..\n'))
    # dmn = r_run_time + nmins * 60 
    # dmx = r_run_time
    dmn = min(p$date); dmx = max(p$date)
    eno_lst = load_eno(timestr, invent = aq_invent, emiss_fn = eno_fn, epa_fn, 
                       epa_name, epa_tz, xmin, xmax, ymin, ymax, dmn, dmx)
    eno_rt = eno_lst$emiss_rt    # EDGAR eNO
    eno_pp = eno_lst$emiss_pp    # EDGAR eNO for PP if needs to be scaled 


    # load EDGAR CO emissions ----------------------------------------------
    eco_rt = load_eco(timestr, invent = aq_invent, emiss_fn = eco_fn, 
                      xmin, xmax, ymin, ymax)


    # ----------------------------------------------------------------------
    # if ODIAC, calc ENO based on ECO2_ODIAC and ERNO_EDGAR, DW, 03/03/2022
    if (aq_invent == 'odiac') {

        cat('***ENOx from ODIAC is calc from ODIAC ECO2 and ER from EDGAR...\n')
        eco2_odiac = eco2_lst$odiac_rt
        eco_edgar  = eco_rt 
        eco2_edgar = eco2_rt 
        eno_edgar  = eno_rt 
        
        # EDGAR based emission ratio
        erno_edgar = eno_edgar / eco2_edgar         # EDGAR based ER
        erco_edgar = eco_edgar / eco2_edgar 
        erno_prj = projectRaster(erno_edgar, eco2_odiac, method = 'ngb')
        erco_prj = projectRaster(erco_edgar, eco2_odiac, method = 'ngb')
        eno_odiac = eco2_odiac * erno_prj     # ODIAC-based ENO
        eco_odiac = eco2_odiac * erco_prj     # ODIAC-based ECO
        eno_odiac[eno_odiac < 0] = 0
        eco_odiac[eco_odiac < 0] = 0
        
        # use ODIAC and ODIAC based emission
        eno_rt = eno_odiac
        eco_rt = eco_odiac
        eco2_rt = eco2_odiac
    }


    # ----------------------------------------------------------------------
    # for lifetime look up relationship, need to know SZA -----------------
    # grab initial footprint, not the weighted footprint - need to re-weight
    cat('p_emiss_prep(): subset particles and add SZA...\n') 
    if ('foot_before_weight' %in% colnames(p)) 
        p = p %>% dplyr::select(-foot) %>% rename(foot = foot_before_weight)

    # 0.5. select trajectories for speeding up, only 12 hours back is sufficient
    library(GeoLight)
    pex = p %>% filter(time >= nmins, long >= xmin, long <= xmax, lati >= ymin, 
                       lati <= ymax, !is.na(recp_lat)) %>% 
          dplyr::select(indx, time, date, lati, long, zagl, zsfc, foot, mlht, 
                        pres, temz, dens, xhgt, any_of(c('sphu', 'rhfr'))) %>% 
          mutate(psza = zenith(solar(date), long, lati)) %>% arrange(time)


    # assign emission to paticle locations using raster::extract --------------
    cat('p_emiss_prep(): assigning emissions to particles...\n')
    p_sp = pex[, c('long', 'lati')]
    coordinates(p_sp) = ~long + lati; crs(p_sp) = crs(eno_rt)
    pex$eno = raster::extract(eno_rt, p_sp, method = 'simple') 
    pex$eco = raster::extract(eco_rt, p_sp, method = 'simple') 
    pex$eco2 = raster::extract(eco2_rt, p_sp, method = 'simple') 
    pex = pex %>% mutate(eno = ifelse(is.na(eno), 0, eno), 
                         eco = ifelse(is.na(eco), 0, eco), 
                         eco2 = ifelse(is.na(eco2), 0, eco2))
    

    # ----------------------------------------------------------------------
    # only scale the max value for PP based on EPA, DW, 09/09/2021
    # correct for bottom-up prior emissions if needed (for optimization)
    if ( aq_invent == 'epa' ) {    # only scale EDGAR emiss for PP based on EPA

        epa_co2_df = eco2_lst$epa_df 
        epa_no_df = eno_lst$epa_df 
        if (is.null(epa_co2_df) | is.null(epa_no_df)) 
            stop('p_emiss_prep(): PP emiss or EPA file is missing when trying to scale emissions based on EPA...\n')

        epa_df = full_join(epa_co2_df %>% rename(epa_co2_pp = epa_pp), 
                           epa_no_df %>% rename(epa_no_pp = epa_pp), 
                           by = c('timestr_lst', 'timestr_utc'))

        # re-locate the emission grid where PP falls into 
        # x, y coordinate in lower left corner now, DW, 03/22/2022
        names(eno_rt) = 'enox'
        eno_res = res(eno_rt)[1]; eco2_res = res(eco2_rt)[1]
        eno_xy = as.data.frame(eno_rt, xy = T) %>% filter(enox == eno_pp) %>%
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
    if (perturbTF) {
        
        cat('p_emiss_prep(): obtaining perturbation scaling factor...\n')
        
        # load perturbation scaling factor, DW, 04/15/2022
        sf_df = readRDS(perturb_fn) %>% filter(n == perturb_indx)
        sf_rt = rasterFromXYZ(sf_df[, c('x', 'y', 'sf')])
        pex$sf_perturb = raster::extract(sf_rt, p_sp, method = 'simple') 

        # if NA meaning no need to perturb the NOx emissions
        # also perturb CO2 and CO emissions if perturbTF = T, DW, 08/11/2022 
        pex = pex %>% mutate(sf_perturb = ifelse(is.na(sf_perturb), 
                                                 1, sf_perturb), 
                             #eno_init = eno,    # preserve original emission
                             eno = eno * sf_perturb, 
                             eco = eco * sf_perturb, 
                             eco2 = eco2 * sf_perturb) 
    }


    ### -----------------------------------------------------------------------
    # use initial footprint, AK, and a normalized PWF, DW, 08/25/2021
    # default PWF is weighted based on the entire column, modify it
    g = 9.8        	         # m s-2
    Mair = 29 / 1E3        	 # kg mol-1 

    # tropospheric column according to TROPOMI tropopause
    no2_all = get.wgt.tropomi.func(output, no2_fn, 'NO2')
    if ( is.null(no2_all) ) return()
    
    no2_df  = no2_all$combine.prof %>% filter(stiltTF)
    no2_info = no2_all$tropomi.info

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
    #no2_info$air_vcd_total = 

    if (!is.na(co_fn)) {      # also get total column dry air density 
        co_all = get.wgt.tropomi.func(output, co_fn, 'CO', tropo.no2TF = FALSE)
        #if ( is.null(co_all) ) return()    # no need to stop model if CO is NA

        if ( !is.null(co_all)) {
            co_df = co_all$combine.prof 
            co_info = co_all$tropomi.info 
            co_info$air_vcd_total = unique(co_df$xdry.tot)
            co_info$xco_ppb = co_info$co_vcd / co_info$air_vcd_total * 1e9

            # interpolate AK of XCO to TROPOMI NO2 levels
            wgt_df = wgt_df %>% 
                    mutate(ak_xco = approx(co_df$lower_pres, co_df$ak.norm, 
                                            xout = wgt_df$lower_pres)$y)
        } else co_info = NULL
    } else co_info = NULL # end if


    if (!is.na(co2_fn)) {   # get xCO2 ak and aprior profiles from OCO
        co2_info = get.oco.info(receptor = receptor, oco.fn = co2_fn, 
                                diff.td = 1e3) 
        
        # calculate dry air column density of OCO columns
        co2_info$air_vcd_total = co2_info$oco.psfc * 100 / g / Mair - co2_info$oco.xh2o
        print(co2_info$air_vcd_total)
        wgt_df = wgt_df %>% mutate(ak_xco2 = approx(co2_info$pres, 
                                                    co2_info$ak.norm, 
                                                    xout = wgt_df$lower_pres)$y,
                                   ap_xco2 = approx(co2_info$pres, 
                                                    co2_info$ap, 
                                                    xout = wgt_df$lower_pres)$y)
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
                 demiss_co2 = eco2 * foot_wgt * eco2_sf) %>% 
          dplyr::select(-c('pwf_norm'))
    

    # check delt - default is 1 min from trajec file
    times = unique(pex$time)
    delt = min(abs(unique(diff(times))))  # in min
    if (max(pex$time) != -delt) 
        stop(paste0('p_emiss_prep(): missing near-field trajec e.g., with time of -', delt, ' min.\nPLEASE check your xmin, xmax, ymin, ymax...DO NOT subset near-field trajec!!!\n'))

    return(list(pex = pex, no2_info = no2_info, co_info = co_info, 
                           co2_info = co2_info))
}




#debugging....
if (F) {

    # cut particles beyond emission grid
    # remove zero footprint to save time and space, DW, 01/29/2019 
    emiss = eco_rt

    # find emissions and calculate CO2
    p1 = p %>% filter(long >= extent(emiss)@xmin,
                        long <= extent(emiss)@xmax, 
                        lati >= extent(emiss)@ymin,
                        lati <= extent(emiss)@ymax, foot > 0) 
    p_sp = p1[, c('long', 'lati')]; coordinates(p_sp) = ~long + lati
    crs(p_sp) = crs(emiss)
    p1 = p1 %>% mutate(find.emiss = raster::extract(emiss, p_sp),
                        demiss = find.emiss * foot)

    # also, compute total dCO2 for each traj over all backwards hours
    sum.p1 = p1 %>% group_by(indx) %>% na.omit() %>%
                dplyr::summarize(ff = sum(demiss)) %>% ungroup()


    # current calculation ---------
    p_sp = p[, c('long', 'lati')]; coordinates(p_sp) = ~long + lati
    crs(p_sp) = crs(emiss)

    p2 = p %>% mutate(find.emiss = raster::extract(emiss, p_sp), 
                        find.emiss = ifelse(is.na(find.emiss),0, find.emiss), 
                        demiss = find.emiss * foot)
    sum.p2 = p2 %>% group_by(indx) %>% dplyr::summarise(ff = sum(demiss))

}
