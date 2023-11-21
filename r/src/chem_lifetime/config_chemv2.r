
get_timestr_from_outpath = function(store_path, met) {

    out_path = list.files(store_path, 'out_20', full.names = T)
    out_path = out_path[grepl('TROPOMI', out_path) & grepl(met, out_path)]
    no2_dirs = list.files(out_path, 'NO2', full.names = T)
    all_timestr = unique(sort(strsplit.to.df(basename(dirname(no2_dirs)))$V2))
    return(all_timestr)
}


config_chemv2 = function(namelist) {

    site       = namelist$site 
    timestr    = namelist$timestr 
    aq_invent  = namelist$aq_invent 
    ghg_invent = namelist$ghg_invent 
    mx_res     = namelist$mx_res 
    mx_hr      = namelist$mx_hr
    ts_chem    = namelist$ts_chem
    lon_lat    = namelist$lon_lat[[1]]
    store_path = namelist$store_path 
    overwriteTF = namelist$overwriteTF
    tno2_path  = namelist$tno2_path 
    tno2x_path = namelist$tno2x_path
    xco_path  = namelist$xco_path 
    xco2_path = namelist$xco2_path 
    xch4_path = namelist$xch4_path
    perturb_emiTF = namelist$perturb_emiTF
    perturb_tsTF  = namelist$perturb_tsTF
    perturb_mixTF = namelist$perturb_mixTF
    n_perturb  = namelist$n_perturb
    transerrTF = namelist$transerrTF

    # locate correct EPA file  -------------------------
    epa_tz    = lon_lat$tz             # will convert to standard time
    epa_state = tolower(lon_lat$state_abb)
    epa_name  = namelist$epa_name
    epa_path  = namelist$epa_path
    if (ghg_invent == 'epa' | aq_invent == 'epa') {
        epa_fn = file.path(epa_path, paste0(substr(timestr, 1, 4), epa_state, 
                                            substr(timestr, 5, 6), '.csv'))
        if (length(epa_fn) == 0) stop('NO hourly EPA data found...\n')
    } else epa_fn = NA

    # check if trajec and satellite data are available before simulations
    out_path = list.files(store_path, 'out_20', full.names = T)
    if ( transerrTF ) 
        out_path = list.files(store_path, 'outerr_20', full.names = T)

    out_path = out_path[grepl('TROPOMI', out_path) & 
                        grepl(met, out_path) & 
                        grepl(timestr, out_path)]
    
    byid_path = file.path(out_path, 'NO2/by-id')
    fns_lst = prep_satelite_fns(timestr, byid_path, tno2_path, tno2x_path, 
                                xch4_path, xco_path, xco2_path, lon_lat)
    traj_info = all_info = fns_lst$traj_info
    tno2_fn = fns_lst$no2_fn
    xch4_fn = fns_lst$ch4_fn
    xco2_fn = fns_lst$co2_fn
    xco_fn  = fns_lst$co_fn
    
    # before simulation, create ENO, ECO for non-EDGAR inventories ----------
    cat('\n\npreparing emissions in advance...\n')
    emis_lst = preproc_emiss(site, timestr, traj_info, ghg_invent, aq_invent,
                             eco2_edgar_fn = namelist$eco2_fn, 
                             eco2_nonedgar_fn = namelist$eco2_fn2, 
                             eco_edgar_fn = namelist$eco_fn, 
                             eno_edgar_fn = namelist$eno_fn, 
                             ech4_edgar_fn = namelist$ech4_fn,
                             epa_fn, epa_name, epa_tz, out_path, 
                             overwriteTF = overwriteTF) 
    eno_fn  = emis_lst$eno_fn
    eco_fn  = emis_lst$eco_fn
    eco2_fn = emis_lst$eco2_fn
    ech4_fn = emis_lst$ech4_fn
    epa_lst = emis_lst$epa_lst


    # prep for perturbation runs -------------------------
    # use the same perturbation for all overpasses
    if ( perturb_emiTF ) {
        
        cat('creating perturbation sf for emiss...\n\n')
        inv_path = file.path(store_path, 'inversion')
        dir.create(inv_path, showWarnings = F)
        perturb_fn = file.path(inv_path, 
                               paste0('perturb_emiss_norm_', site, '.rds'))

        # default prior emiss uncert from EDGAR vs. ODIAC CO2
        # TO DO +++
        sf_df = create_prior_distri(site, timestr, eno_edgar_fn = eno_fn, 
                                    eco2_edgar_fn, eco2_odiac_fn, inv_path, lon_lat, n_perturb, F)
        
    } else if ( perturb_tsTF ) {
        
        # generate scaling factor for lifetime perturbation
        perturb_fn = file.path(store_path, 
                              paste0('perturb_ts_', site, '_', timestr, '.csv'))

        ts_sd = 0.41
        if ( !file.exists(perturb_fn) ) {
            check_sd = ts_sd + 0.5
            check_mean = 1.5
            i = 0
            while( abs(check_sd - ts_sd) > 0.01 | abs(check_mean - 1) > 0.01 ) {
                sf_df = data.frame(n = 1: n_perturb, 
                                   sf = sort(rnorm(n_perturb, mean = 1, 
                                                   sd = ts_sd)))
                check_sd = sd(sf_df$sf)
                check_mean = mean(sf_df$sf)
                i = i + 1
            }
            write.csv(sf_df, perturb_fn, row.names = F)   
        } else sf_df = read.csv(perturb_fn)

    } else if ( perturb_mixTF ) {
        
        # different mixing length scales and time scales
        perturb_fn = file.path(store_path, 
                          paste0('perturb_mixing_', site, '_', timestr, '.csv'))
        
        if ( !overwriteTF & file.exists(perturb_fn) ) {
            mx_df = read.csv(perturb_fn)

        } else {
            
            # test a range of mixing present in GMD2023
            mx_res = c(1, 3, 10)         # in km
            dffs = c(1e2, 1e3, 1e4)     # eddy diffusivity in m2 s-1
            mx_df  = expand.grid(mx_res = mx_res, mx_hr = mx_hr) %>% 
                     #mutate(mx_hr = signif((mx_res *1e3) ^2 / dff / 3600 * exp(-1), 1)) %>% 
                     filter(mx_res * mx_hr != 1) %>% 
                     tibble::rownames_to_column(var = 'n') %>% 
                     mutate(n = as.numeric(n))

            n_perturb = max(mx_df$n)     # number of perturbation 
            write.csv(mx_df, file = perturb_fn, row.names = F)
        }   # end if

    } else {
        # if no perturbation, set # of perturb param as 1
        perturb_fn = NA
        n_perturb = 1
    }   # end if

    
    # ------------------------------------------------------------------------
    # loop over each run, either best-estimation or perturbation for inversions
    # ------------------------------------------------------------------------
    for ( perturb_indx in 1: n_perturb ) {
        
        # double check for non-perturbation runs
        if ( !perturb_emiTF & !perturb_mixTF & !perturb_tsTF ) perturb_indx = NA
        
        # slurm job names
        jobname = paste0(aq_invent, '_', substr(timestr, 3, 8), '_', 
                         site, '_', met)
        if ( transerrTF ) jobname = paste0('transerr_', jobname)
        if ( !is.na(perturb_indx) ) 
            jobname = paste0(aq_invent, perturb_indx, '_', 
                             substr(timestr, 3, 8), '_', site, '_', met)

        # get output filename  -------------------------
        rds_patt = naming_fnsv2(aq_invent, ghg_invent, mx_res, namelist$eno_sf, 
                                ts_chem, perturb_indx, perturb_emiTF, 
                                perturb_tsTF, perturb_mixTF)
        
        # -------------------------
        slurm = namelist$slurm
        n_nodes = namelist$n_nodes
        n_cores = namelist$n_cores
        
        # see if simulations exist
        rds_fns = list.files(byid_path, rds_patt, full.names = T, recursive = T)
        if ( length(rds_fns) == 0 ) overwriteTF = T 
        
        if ( !overwriteTF ) {       # not re-run, try postprocessing first
            cat('start post-processing...\n')
            run_msTF = p_post_proc(site, timestr, byid_path, met, rds_patt,
                                   perturb_emiTF, perturb_tsTF, perturb_mixTF,
                                   transerrTF, overwriteTF, store_path)
            
            # modified time needs to be after this min time
            if (run_msTF) {
                mtime_min = '2023-04-01 00:00:00'        
                ms_info = subset_recp(byid_path, rds_patt, all_info, mtime_min) 
                traj_info = ms_info
                cat(paste('\n\n** Need to rerun', nrow(traj_info), 
                          'missing/outdated receptors...\n'))
                n_cores = ceiling(nrow(traj_info) / n_nodes)
                jobname = paste0(jobname, '_ms')
            }   # end if

        } else {    # if re-run for all receptors, delete previous versions
            file.remove(rds_fns)
        }   # end if 

        if (n_cores > 10) n_cores = 10
        if (!slurm) n_nodes = n_cores = 1       # for debugging purpose
        slurm_options = list(time = '06:00:00', account = 'pow', 
                             partition = 'any', 
                             mem = namelist$mem_per_node * n_nodes * 1024)

        setwd(namelist$slurm_path)
        xstilt_apply(FUN = xgas_solverv4, 
                     slurm = slurm, 
                     slurm_options = slurm_options, 
                     n_nodes = n_nodes, 
                     n_cores = n_cores, 
                     jobname = jobname, 
                     site    = site, 
                     traj_fn = traj_info$fn, 
                     timestr = traj_info$time, 
                     bg_nox  = namelist$bg_nox, 
                     eno_fn  = eno_fn, 
                     eco_fn  = eco_fn, 
                     ech4_fn = ech4_fn, 
                     eco2_fn = eco2_fn, 
                     epa_lst = epa_lst,
                     eno_sf  = namelist$eno_sf,
                     eco_sf  = namelist$eco_sf,
                     ech4_sf = namelist$ech4_sf, 
                     eco2_sf = namelist$eco2_sf,
                     perturb_emiTF = perturb_emiTF, 
                     perturb_tsTF  = perturb_tsTF,
                     perturb_mixTF = perturb_mixTF, 
                     perturb_fn    = perturb_fn, 
                     perturb_indx  = perturb_indx, 
                     tno2_fn    = tno2_fn,
                     tno2x_path = tno2x_path, 
                     xco_fn     = xco_fn, 
                     xco2_fn    = xco2_fn, 
                     xch4_fn    = xch4_fn, 
                     aq_invent  = aq_invent,
                     ghg_invent = ghg_invent, 
                     xmin = NA, 
                     xmax = NA, 
                     ymin = NA, 
                     ymax = NA, 
                     met = met, 
                     mx_hr   = mx_hr,
                     mx_res  = mx_res, 
                     ts_fn   = namelist$ts_fn, 
                     ts_chem = ts_chem, 
                     transerrTF = transerrTF, 
                     truncTF = T, 
                     xstilt_wd = namelist$xstilt_wd)
    }   # end for
    
}

if (F) {

    X = 758
    traj_fn = traj_info$fn[X]
    timestr = traj_info$time[X]
    aq_invent = ghg_invent = 'edgar'
    eno_sf = eco2_sf = eco_sf = ech4_sf = 1
    ts_chem = NA
    truncTF = T
    xmin = xmax = ymin = ymax = NA 
    bg_nox = NA 
}
