
p_tau_icbc_prep = function(pex, tau_chem = NA, tau_fn = NA, perturb_fn = NA, 
                           perturb_tsTF, perturb_indx, bg_nox = NA, aux_path, 
                           xstilt_wd) {

    if ( !is.numeric(tau_chem) ) {  # using lifetime curves
        cat('load_tau_icbc_mix(): loading lifetime relationship...\n')
        
        if ( is.na(tau_fn) ) tau_fn = file.path(xstilt_wd, 'data/nox_curves_gapfill_v20230503.rds')
        
        tau_lst = readRDS(tau_fn) 
        tau_nox_df = tau_lst$nox %>% rename(n_nox = n)

        # perturb lifetime if needed, DW, 04/28/2022
        if ( !is.na(perturb_fn) & perturb_tsTF ) {
            tau_sf = read.csv(perturb_fn) %>% filter(n == perturb_indx)
            tau_nox_df = tau_nox_df %>% mutate(ts_nox = ts_nox * tau_sf$sf)
        }
        
    } else tau_nox_df = NULL  # end if tau

    ### ------------------------------------------------------------------------
    if ( is.na(bg_nox) ) {
        cat('xgas_solver(): assigning NO2 prior to trajec endpoint...\n')
        p_edpt = pex %>% group_by(indx) %>% filter(time == min(time)) %>% 
                 dplyr::select(indx, lati, long, time, run_time = date, temz, 
                               pres, one_of('sphu','rhfr')) %>% 
                 unique() %>% ungroup() 

        # assign TROPOMI TM5 prior to trajec endpoint as IC/BC or background
        p_no2i = get.tropomi.prior.no2(p_edpt, aux_path) 
        #plot(p_no2i$no2_dry * 1e9, p_no2i$pres, ylim = c(1013, 0)) 
        
        # IC/BC NO2 now in ppm 
        pex = pex %>% 
              left_join(p_no2i %>% dplyr::select(indx, no2_ic = no2_dry), 
                        by = 'indx') %>% mutate(no2_ic = no2_ic * 1e6) 
    }   # end if IC/BC

    load_list = list(pex = pex, tau_nox_df = tau_nox_df)
    return(load_list)
}