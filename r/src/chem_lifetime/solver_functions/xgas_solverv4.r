
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
# add additional column for no-mixing concentration, DW, 11/09/2021
# ts_mix for typical mixing time scale in PBL, DW, 01/26/2022 
# further modulize the code, DW, 11/21/2022 
# adding CH4 simulation, DW, 11/20/2023 

xgas_solverv4 = function(site, traj_fn, timestr, bg_nox = NA, eno_fn, eco_fn, 
                         ech4_fn, eco2_fn, epa_lst, eno_sf = 1, eco_sf = 1, 
                         ech4_sf = 1, eco2_sf = 1, perturb_emiTF = F, 
                         perturb_tsTF = F, perturb_mixTF = F, perturb_fn = NA,
                         perturb_indx = NA, tno2_fn, tno2x_path, xco_fn = NA,
                         xch4_fn = NA, xco2_fn = NA, 
                         aq_invent = c('edgar', 'epa', 'odiac', 'vulcan')[1], 
                         ghg_invent = c('edgar', 'epa', 'odiac', 'vulcan')[1], 
                         xmin = NA, xmax = NA, ymin = NA, ymax = NA, met, 
                         mx_hr = NA, mx_res = NA, ts_fn = NA, ts_chem = NA, 
                         transerrTF = F, truncTF = T,
                         xstilt_wd = '/central/home/dienwu/models/X-STILT'){

    start_time = Sys.time()
    setwd(xstilt_wd); source('r/dependencies.r')    # source all functions
    cat(traj_fn)
    
    ### ------------------------------------------------------------------------
    # 1. call functions for loading emissions and preparing particles
    # it will also calc dry-air column density and observed Mixing ratio
    # Rprofiling ----- this chunk takes ~20 sec at most, DW, 11/29/2022
    cat('\n\n# ---- xgas_solverv4(): Loading emissions and trajec ---- #\n')
    try_error = try(readRDS(traj_fn), silent = T)
    if ( 'particle' %in% names(try_error) & 'qt_prof' %in% names(try_error) ) {
        output = try_error
    } else if (attributes(try_error)$condition$message == 'error reading from connection') {
        cat('xgas_solverv4(): trajec file is somehow broken, error reading it...\n\n')
        return()
    } else {
        cat('xgas_solverv4(): printing error messages ---\n')
        print(try_error)
    }   # end if

    # double check for non-perturbation runs
    if ( !perturb_emiTF & !perturb_mixTF & !perturb_tsTF ) perturb_indx = NA 
    
    #print(str(output$qt_prof))
    receptor = output$receptor
    recp_lon = strsplit.to.df(basename(traj_fn))$V2
    recp_lat = strsplit.to.df(basename(traj_fn))$V3
    p = output$particle %>% filter(!is.na(lati), !is.na(long)) %>%
        mutate(recp_lon = recp_lon, recp_lat = recp_lat, 
               date = receptor$run_time + time * 60) 
    
    prp = p_emiss_prep(site, output, p, timestr, aq_invent, ghg_invent, eno_fn,
                       eco_fn, ech4_fn, eco2_fn, epa_lst, eno_sf, eco_sf, 
                       ech4_sf, eco2_sf, tno2_fn, xco_fn, xch4_fn, xco2_fn, 
                       perturbTF = perturb_emiTF, perturb_indx, perturb_fn)
    
    if (is.null(prp)) return()
    pex = prp$pex %>% arrange(abs(time), indx) 
    tno2_info = prp$no2_info
    xco_info  = prp$co_info
    xch4_info = prp$ch4_info
    xco2_info = prp$co2_info 


    ### ------------------------------------------------------------------------
    # 2. grab NO2 and NOx net loss timescale [hr] and
    #    grab aux NO2 prior profile if needed, DW, 2022/01/13 
    # Rprofiling ----- this chunk takes 20 sec at most, DW, 11/29/2022
    cat('\n\n# ---- xgas_solverv4(): Loading NOx lifetime curves and [NO2] initial conditions from TM5 ---- #\n')
    prp2 = p_ts_icbc_prep(pex, ts_chem, ts_fn, perturb_fn, perturb_tsTF, 
                          perturb_indx, bg_nox, tno2x_path, xstilt_wd)
    pex = prp2$pex 
    ts_nox_df = prp2$ts_nox_df
    min_nox = ifelse( is.numeric(ts_chem), 1e-6, min(ts_nox_df$bin_nox) )
    
    # if multiple mixing length scales and time scales are given
    if ( perturb_mixTF & !is.na(perturb_fn) ) {        
        mx_df = read.csv(perturb_fn) %>% filter(n == perturb_indx)
        mx_res = mx_df$mx_res 
        mx_hr = mx_df$mx_hr
    }   # end if


    ### ------------------------------------------------------------------------
    # 3. loop over each backward time and simulate concentrations
    # Rprofiling ----- this chunk 6 to 20 mins depending on dt, DW, 11/29/2022
    cat('\n\n# ---- xgas_solverv4(): Calculating [C] per timestamp ---- #\n')
    p_chem = p_solverv4(pex, mx_res, mx_hr, bg_nox, ts_nox_df, ts_chem, min_nox)
    cat('# ---- xgas_solverv4(): DONE with forward calculation ---- #\n\n')


    ### ------------------------------------------------------------------------
    # 4. perform vertical AK PWF weighting to particle-level concentrations
    cat('\n\n# ---- xgas_solverv4(): Performing vertical weighting ---- #\n')
    p_recp = x_weightingv4(p_chem, tno2_info, xco_info, xch4_info, xco2_info, 
                           output)
    
    # select variables to save memory 
    if (truncTF) {  
        colnms = c('indx', 'xhgt', 'time', 'date', 'psza', 'foot_wgt', 'temz', 
                   'dens', 'eno', 'eco2', 'eco', 'ech4', 'ak_tno2', 'no2_ic', 
                   'p_nox_nochem', 'p_nox_mix', 'p_nox_nomix', 'p_no2_mix', 
                   'ts_nox_mix', 'p_rto_mix', 'p_co2_mix', 'p_co_mix', 
                   'p_ch4_mix', 'dchem_mix','p_OX_mix')
        p_chem = p_chem[, colnames(p_chem) %in% colnms] 
        p_recp = p_recp %>% dplyr::select(indx, lati, long, pres, psza, xhgt, 
                                          foot_wgt, date, no2_ic, 
                                          contains('ak_'), ends_with('_wgt'), 
                                          starts_with('p_no2_')) 

        if (perturb_emiTF | perturb_tsTF | perturb_mixTF) 
            perturb_info = list(perturb_indx = perturb_indx, 
                                perturb_fn = perturb_fn)
    } # end if


    # 5. store output --------------------------------------------------
    cat('\n\n# ---- xgas_solverv4(): Storing output to .rds files ---- #\n')
    p_list = list(p_chem = p_chem, p_recp = p_recp, no2_info = tno2_info,
                  co_info = xco_info, ch4_info = xch4_info, 
                  co2_info = xco2_info)
    
    if (perturb_emiTF | perturb_tsTF | perturb_mixTF) 
        p_list = list(p_recp = p_recp, perturb_info = perturb_info, 
                      no2_info = tno2_info, co_info = xco_info, 
                      ch4_info = xch4_info, co2_info = xco2_info)

    rds_patt = naming_fnsv2(aq_invent, ghg_invent, mx_res, eno_sf, ts_chem, 
                            perturb_indx, perturb_emiTF, perturb_tsTF, 
                            perturb_mixTF) 
    
    # if not the last receptor that is working on
    # store the results temporarily
    common_fn = paste0(timestr, '_', recp_lon, '_', recp_lat)
    rds_fn = file.path(dirname(traj_fn), paste0(common_fn, rds_patt))
    saveRDS(p_list, file = rds_fn) 

    cat('xgas_solverv4(): Done with NOx calculation...;')
    end_time = Sys.time()
    cat(paste('Time used:', signif(end_time - start_time, 2), 'mins \n\n\n'))
}   # end of function

