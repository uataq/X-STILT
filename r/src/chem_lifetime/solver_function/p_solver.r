
# now calc NO2-to-NOx ratio using NOx-scaled OX and J(NO2) and K(OH+O3)
# DW, 02/07/2022

# remove [NOx] before and after, just stick to one data colume for [NOx]
# refactoring the loop for solving for particle-level concentration, 
# (avoid appending particle results to the data.frame; 
#  recycle the particle results of the current timestamp for time = t + 1)
# DW, 11/29/2022
p_solver = function(pex, mx_res, mx_hr, uni_sza, bg_nox, tau_nox_df, 
                    tau_chem = c('v20230503', 4, NA)[3], min_nox = 1e-6) {
    
    times = sort(unique(pex$time))  # increasing order
    delt = min(abs(unique(diff(times))))
    OX_df = readRDS('/home/dienwu/postdoc_proj/NOx/rds/ox_coef_bgint.rds')
    
    if (is.na(mx_res)) {
        mx_res = 5      # in km
        cat(paste('p_solver(): use TROPOMI grid spacing as mixing length of',   
                  mx_res, 'km...\n'))
    }
    
    p_fill = pex %>% tibble::add_column(p_nox_nochem = NA, p_nox_nomix = NA, 
                                        p_nox_mix = NA, p_co2_nomix = NA, 
                                        p_co2_mix = NA, p_co_nomix = NA, 
                                        p_co_mix = NA, p_rto_nomix = NA, 
                                        p_rto_mix = NA, p_no2_nomix = NA, 
                                        p_no2_mix = NA, p_OX_mix = NA, 
                                        dh = NA, dentrain_mix = NA, 
                                        dentrain_nomix = NA, 
                                        dentrain_nochem = NA,
                                        ts_nox_nomix = NA, dchem_nomix = NA, 
                                        ts_nox_mix = NA, dchem_mix = NA, 
                                        n_mix = NA) 
    cols_fill = colnames(p_fill)[!colnames(p_fill) %in% colnames(pex)]

    dt_prof = p_last = NULL  
    for (t in 1 : length(times)) {
        
        if (t %% 50 == 0) 
            cat(paste('# ---', signif(t / length(times), 3) * 1E2, '% --- #\n'))
        #time0 = Sys.time()
  
        # 0 - subset particle df for the current timestamp, 
        # as for the particle at last timestamp, we will recycle from last loop
        row_indx = which(pex$time == times[t])
        p_curr = pex %>% filter(time == times[t]) 
        #time1 = Sys.time()
        
        # 1 - initialize with info from last time stamp
        # p_last will be NULL for the first timestamp (i.e., endpoint)
        p_init = p_initial(p_curr, p_last, mx_res, uni_sza, bg_nox, a_rto = 1)
        #time2 = Sys.time()

        # 2 - if MLHT increase over time, perform vertical mixing --------------
        # min_nox as background concentration above mixed layer
        p_ent = p_entrain(p_init, delt, dilutionTF = T, min_nox)
        #time3 = Sys.time()

        # 3 - calculate changes in nox due to chem -----------------------------
        p_tau = p_chemicalv3(p_ent, delt, tau_nox_df, OX_df, tau_chem, 
                             num_chunk = 1)
        
        # double check if there is negative NOx 
        if ( length(p_tau$p_nox_nomix[p_tau$p_nox_nomix < 0]) > 1 ) 
            stop(paste('p_solver(): capture particles with negative [NOx] for backward time = ', times[t], 'min...\n') )
        #time4 = Sys.time()

        # 4. horizontal mixing of [NOx] among particles
        p_mix = p_mixing(p_tau, mx_hr, delt)
        #time5 = Sys.time()

        # 5. subset particle results and insert them back to the data frame
        p_sub = p_mix %>% dplyr::select(all_of(cols_fill))
        p_fill[row_indx, cols_fill] = p_sub

        # recycle the particle results at time = t for the next loop
        # `p_last` here will be used for the next loop when time = t + 1
        p_last = p_mix %>% dplyr::select(indx, mlht, starts_with('p_'))
        #time6 = Sys.time()

        # dt_prof = rbind(dt_prof, data.frame(timestamp = times[t], dt1 = time1 - time0, dt2 = time2 - time1, dt3 = time3 - time2, dt4 = time4 - time3, dt5 = time5 - time4, dt6 = time6 - time5))
    }   # end for t 

    # write.table(dt_prof, file = '/home/dienwu/postdoc_proj/NOx/run_scripts/Rprof_rev.txt', quote = F, row.names = F)
    
    # remove all *_mean columns to reduce size
    p_chem = p_fill %>% dplyr::select(-ends_with('_mean')) 
    return(p_chem)
}
