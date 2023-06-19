
### ----------------------------------------------------------------------
# start with forward calculation from the min time, from -1440 mins to -1 
# assuming the BC for the furthest particle is zero 
# 0. figure out the emission for all particles
# 1. calculate the N mixing ratio before chem 
# 2. estimate the lifetime according to [NOx]
# 3. calculate the scaling factor - i.e., net changes in [NOx], L - P 
# 4. scale [NOx] before the chem to update the [NOx]
# 5. move on to the next timestamp (delt << longest lifetime)

# traj_fn, timestr, recp_lon, recp_lat can be a vector, with same dim, DW
# add vertical dilution - if mixed layer height increases over time, 
# then add extra term of (Ci_aloft - Ci) / H(t) * (dH / delt), DW, 08/18/2021
# now remove chemTF, always calc C with lifetime and without, DW, 10/25/2021 
# add typical O3 of 50 ppb, update NO2/NOx if NO2 exceed 50 ppb
# nmins for how many mins back 
# add additional column for no-mixing concentration, DW, 11/09/2021
# tau_mix for typical mixing time scale in PBL, DW, 01/26/2022 
# further modulize the code, DW, 11/21/2022 
#'@param tau_chem either a constant number or NA for using NOx curves 

if (F) {

    X = 1   # 92
    traj_fn = traj_info$fn[X]
    timestr = traj_info$time[X]
    recp_lon = traj_info$lon[X] 
    recp_lat = traj_info$lat[X]
    aq_invent = ghg_invent = 'edgar'; epa_fn = NA
    #aq_invent = ghg_invent = 'epa'
    eno_sf = eco2_sf = eco_sf = 1
    tau_chem = NA
    tau_fn = '/home/dienwu/postdoc_proj/NOx/tau_rto_gapfill_v20211008.rds'
    nmins = -720
    mx_scheme = 'partial'
    perturb_emissTF = F
    perturb_tsTF = F 
    perturb_mixTF = F 
    perturb_indx = NA
    perturb_fn = NA
    truncTF = T
    xmin = xmax = ymin = ymax = NA 
    bg_nox = NA 
    
}

xgas_solver = function(site, traj_fn, timestr, recp_lon, recp_lat, bg_nox = NA,
                       eno_fn, eco_fn, eco2_edgar_fn, eco2_odiac_fn = NA,
                       eno_sf = 1, eco_sf = 1, eco2_sf = 1, perturb_emissTF = F,
                       perturb_tsTF = F, perturb_mixTF = F, perturb_fn = NA,
                       perturb_indx = NA, no2_fn, aux_path, co_fn = NA, 
                       co2_fn = NA, epa_fn = NA, epa_name = NA, epa_tz = 'PDT',
                       aq_invent = c('edgar', 'epa', 'odiac')[1], 
                       ghg_invent = c('edgar', 'epa', 'odiac')[1], nmins = -720,
                       xmin = NA, xmax = NA, ymin = NA, ymax = NA, mx_res = NA, 
                       mx_hr = 3, tau_fn = NA, tau_chem = NA, storeTF = F,
                       truncTF = T, xstilt_wd = '/central/home/dienwu/X-STILT'){

    start_time = Sys.time()
    setwd(xstilt_wd); source('r/dependencies.r')    # source all functions
    cat(traj_fn)

    ### ------------------------------------------------------------------------
    # 1. call functions for loading emissions and preparing particles
    # it will also calc dry-air column density and observed Mixing ratio
    # Rprofiling ----- this chunk takes ~20 sec at most, DW, 11/29/2022
    cat('\n\n# ---- xgas_solver(): Loading emissions and trajectories ---- #\n')
    try_error = try(readRDS(traj_fn), silent = T)
    if ( 'particle' %in% names(try_error) & 'qt_prof' %in% names(try_error) ) {
        output = try_error
    } else if (attributes(try_error)$condition$message == 'error reading from connection') {
        cat('xgas_solver(): trajec file is somehow broken, error reading it...\n\n')
        return()
    } else {
        cat('xgas_solver(): printing error messages ---\n')
        print(try_error)
    }   # end if

    #print(str(output$qt_prof))
    receptor = output$receptor
    p = output$particle %>% filter(!is.na(lati), !is.na(long)) %>%
        mutate(recp_lon = recp_lon, recp_lat = recp_lat, 
               date = receptor$run_time + time * 60) 
    
    prep1 = p_emiss_prep(site, output, p, timestr, eno_fn, eco_fn, 
                         eco2_edgar_fn, eco2_odiac_fn, epa_fn, epa_name, epa_tz,
                         no2_fn, co_fn, co2_fn, eno_sf, eco_sf, eco2_sf, 
                         perturb_emissTF, perturb_indx, perturb_fn, aq_invent, 
                         ghg_invent, xmin, xmax, ymin, ymax, nmins) 
    if (is.null(prep1)) return()
                        
    pex = prep1$pex %>% arrange(abs(time), indx) 
    no2_info = prep1$no2_info
    co_info  = prep1$co_info
    co2_info = prep1$co2_info 

    ### ------------------------------------------------------------------------
    # 2. grab NO2 and NOx net loss timescale [hr] and
    #    grab aux NO2 prior profile if needed, DW, 2022/01/13 
    # Rprofiling ----- this chunk takes 20 sec at most, DW, 11/29/2022
    cat('\n\n# ---- xgas_solver(): Loading NOx lifetime curves and [NO2] initial conditions from TM5 ---- #\n')
    prep2 = p_tau_icbc_prep(pex, tau_chem, tau_fn, perturb_fn, perturb_tsTF, 
                            perturb_indx, bg_nox, aux_path, xstilt_wd)
    
    pex = prep2$pex 
    tau_nox_df = prep2$tau_nox_df
    uni_sza = sort(unique(tau_nox_df$bin_sza))
    min_nox = ifelse( is.numeric(tau_chem), 1e-6, min(tau_nox_df$bin_nox) )
    
    # if multiple mixing length scales and time scales are given
    if (perturb_mixTF & !is.na(perturb_fn)) {        
        mx_df = read.csv(perturb_fn) %>% filter(n == perturb_indx)
        mx_res = mx_df$mx_res 
        mx_hr = mx_df$mx_hr
    }   # end if


    ### ------------------------------------------------------------------------
    # 3. loop over each backward time and simulate concentrations
    # Rprofiling ----- this chunk 6 to 20 mins depending on dt, DW, 11/29/2022
    cat('\n\n# ---- xgas_solver(): Calculating [C] per timestamp ---- #\n')
    p_chem = p_solver(pex, mx_res, mx_hr, uni_sza, bg_nox, tau_nox_df, tau_chem,
                      min_nox)

    cat('# ---- xgas_solver(): DONE with forward calculation ---- #\n\n')
    # pp = left_join(p_chem3 %>% dplyr::select(indx, time, p_nox_mix), 
    #                p_chem4 %>% dplyr::select(indx, time, p_nox_mix), 
    #                by = c('indx', 'time')) %>% 
    #     mutate(p_dff = p_nox_mix.x - p_nox_mix.y)
    # plot(pp$p_nox_mix.x, pp$p_nox_mix.y)

    ### ------------------------------------------------------------------------
    # 4. perform vertical AK PWF weighting to particle-level concentrations
    cat('\n\n# ---- xgas_solver(): Performing vertical weighting ---- #\n')
    p_recp = x_weighting(p_chem, no2_info, co_info, co2_info)
    
    if (truncTF) {  # select variables to save memory 
        colnms = c('indx', 'xhgt', 'time', 'date', 'psza', 'foot_wgt', 'temz', 
                   'dens', 'eno', 'eco2', 'eco', 'ak_tno2', 'no2_ic', 'p_nox_nochem', 'p_nox_mix', 'p_no2_mix', 'ts_nox_mix', 'p_rto_mix', 'p_co2_mix', 'p_co_mix', 'dchem_mix','p_OX_mix')
        p_chem = p_chem[, colnames(p_chem) %in% colnms] 
        p_recp = p_recp %>% dplyr::select(indx, lati, long, pres, psza, xhgt, 
                                          foot_wgt, date, no2_ic, 
                                          contains('ak_'), ends_with('_wgt'), 
                                          starts_with('p_no2_')) 

        if (perturb_emissTF | perturb_tsTF | perturb_mixTF) 
            perturb_info = list(perturb_indx = perturb_indx, 
                                perturb_fn = perturb_fn)
    } # end if


    # 5. store output --------------------------------------------------
    if (storeTF) {
        cat('\n\n# ---- xgas_solver(): Storing output to .rds files ---- #\n')
        p_list = list(p_chem = p_chem, p_recp = p_recp, no2_info = no2_info,
                      co_info = co_info, co2_info = co2_info)

        if (perturb_emissTF | perturb_tsTF | perturb_mixTF) 
            p_list = list(p_recp = p_recp, perturb_info = perturb_info, 
                          no2_info = no2_info)
        
        common_fn = paste0(timestr, '_', recp_lon, '_', recp_lat)
        perturb_emiss_indx = perturb_tau_indx = perturb_mix_indx = NA
        if (perturb_emissTF) perturb_emiss_indx = perturb_indx
        if (perturb_tsTF) perturb_tau_indx = perturb_indx
        if (perturb_mixTF) { perturb_mix_indx = perturb_indx; mx_res = NA }
        
        rds_patt = naming_fns(aq_invent, ghg_invent, mx_res, eno_sf, tau_chem, 
                              perturb_emiss_indx, perturb_tau_indx, 
                              perturb_mix_indx)        
        fn = file.path(dirname(traj_fn), paste0(common_fn, rds_patt))
        saveRDS(p_list, file = fn) 

        # create symbolic links for the simulation 
        
    }   # end if


    cat('xgas_solver(): Done with NOx calculation...;')
    end_time = Sys.time()
    cat(paste('Time used:', end_time - start_time, '\n\n\n'))
}   # end of function

