
# calculate changes in nox due to chem, DW, refactoring on 2021/10/25
# now remove chemTF, always calc C with lifetime and without, DW, 10/25/2021 

#' @param delt the time stamp for solving chemical tendency in minutes. 
#' @param num_chunk the number of sub-min chunks, e.g., if 5, meaning the total enhancements calculated within the default 1 min should be devided by 5 (i.e., 20 sec), in order to update the concentrations due to emissions, chem, etc. 
#' @param num_chunk = 1 for default value, no need to modify demiss_*

# try both ways - 
# 1. calc ratio at current time step using NO2 from last time step
# 2. calc ratio simply using the quadratic equation

# combine subscale calc with p_chemicalv3 --
# get ts for mix and nomix cases, and check if ts << delt = 1min
# if so, start subscale calc...
#' @param tau_chem either numeric for constant lifetime or NA using curves 
p_chemicalv3 = function(p_ent, delt, tau_nox_df, OX_df, 
                        tau_chem = c('v20230503', 4, NA)[3], num_chunk = 1) {
    
    # at this point - p_ent is the particle info for the current time step, 
    # but with all p_* from last time step, which require to be modified by 
    # all forcing/processes 

    # --------------------------------------------------------------------------
    # ------------------------- NON-MIXING cases -------------------------------
    # add calculations for non-mixing calculation first, DW, 11/16/2021
    # now merge lifetime table using non-mixing concentrations 
    # get dchem and remove columns, DW, 11/16/2021
    p_rename1 = p_ent %>% rename(tmp_nox = p_nox_nomix, 
                                 tmp_rto = p_rto_nomix, 
                                 tmp_no2 = p_no2_nomix, 
                                 tmp_dentrain = dentrain_nomix)

    # assign NOx NET change frequency to each particle
    p_tmp1 = lookup_nox_tendency(p_rename1, tau_nox_df, OX_df, delt, tau_chem)

    # check to see if delt >> NOx time scale; delt in min, ts in hr
    # at least delt / ts_nox should be smaller than 0.5
    # if delt << NOx timescale, reduce delt for solving chemistry only when NOx timescale is positive (implying net sink to NOx)
    p_eval = p_tmp1 %>% mutate(eval_t = delt/60/ts_nox) %>% filter(eval_t > 0.5)
    indx_eval = unique(p_eval$indx)

    # --------------------------------------------------------------------------
    # separate particle data.frame for regular 1-min calc and sub-minute calc
    # p_chemical_dtv3() will update 
    # 1) NOx, updated NOx with changes due to emission, chem loss, vert dilution
    # 2) NO2-to-NOx ratio, based on two approaches (1. NO2 from last time stamp; 2. quadratic equation using NOx at the current time stamp)
    # 3) NO2, based on updated ratio and updated NOx

    if ( length(indx_eval) > 0 ) {  
        
        # certain particles have ts much smaller than chemical timestamp = 1min
        # need to carry out sub-minute calculation by reducing chemical timestep
        # subset particles need subdt calc or normal dt calc
        p_subdt = p_tmp1 %>% filter(indx %in% indx_eval) 
        p_dt = p_tmp1 %>% filter(!indx %in% indx_eval)

        p_subdt = p_chemical_subdtv3(p_subdt, p_before = p_rename1, delt, 
                                     tau_nox_df, OX_df, tau_chem)
        p_dt = p_chemical_dtv3(p_df = p_dt, tau_nox_df, delt, num_chunk = 1)

        p_nomix = rbind(p_subdt, p_dt)

    } else  # no need to perform sub-minute calculation 
        p_nomix = p_chemical_dtv3(p_df = p_tmp1, tau_nox_df, 
                                  delt, num_chunk = 1)

    # rename column names back to represent NON-mixing calculations
    # after p_chemical_dtv3, we can remove data from previous time stamp
    p_nomix = p_nomix %>% 
              rename(p_nox_nomix = tmp_nox, p_rto_nomix = tmp_rto, 
                     p_no2_nomix = tmp_no2, p_OX_nomix  = tmp_OX, 
                     ts_nox_nomix = ts_nox, dchem_nomix = dchem, 
                     dentrain_nomix = tmp_dentrain)

    # ------------------------- END OF NON-MIXING cases ------------------------
    # --------------------------------------------------------------------------


    # --------------------------------------------------------------------------
    # ----------------------------- MIXING cases -------------------------------
    # now calculation NOx with chem and mixing within delt
    p_rename2 = p_nomix %>% rename(tmp_nox = p_nox_mix, 
                                   tmp_rto = p_rto_mix, 
                                   tmp_no2 = p_no2_mix, 
                                   tmp_dentrain = dentrain_mix)
    
    # assign NOx NET change frequency to each particle
    p_tmp2 = lookup_nox_tendency(p_rename2, tau_nox_df, OX_df, delt, tau_chem)

    # check to see if delt >> NOx time scale; delt in min, ts in hr
    # at least delt / ts_nox should be smaller than 0.5
    # if delt << NOx timescale, reduce delt for solving chemistry only when NOx timescale is positive (implying net sink to NOx)
    p_eval = p_tmp2 %>% mutate(eval_t = delt/60/ts_nox) %>% filter(eval_t > 0.5)
    indx_eval = unique(p_eval$indx)

    # --------------------------------------------------------------------------
    # SIMILARLY, separate particle data.frame for regular 1-min calc and sub-minute calc
    if ( length(indx_eval) > 0 ) {
        
        # subset particles need subdt calc or normal dt calc
        p_subdt = p_tmp2 %>% filter(indx %in% indx_eval) 
        p_dt = p_tmp2 %>% filter(!indx %in% indx_eval)

        p_subdt = p_chemical_subdtv3(p_subdt, p_before = p_rename2, delt, 
                                     tau_nox_df, OX_df, tau_chem)
        p_dt = p_chemical_dtv3(p_df = p_dt, tau_nox_df, delt, num_chunk = 1)
        p_mix = rbind(p_subdt, p_dt)
        
    } else {  # no need to perform sub-minute calculation 
        p_mix = p_chemical_dtv3(p_df = p_tmp2, tau_nox_df, delt, num_chunk = 1)
    }   # end if 

    p_all = p_mix %>% rename(p_nox_mix = tmp_nox, p_rto_mix = tmp_rto,    
                             p_no2_mix = tmp_no2, p_OX_mix = tmp_OX,
                             ts_nox_mix = ts_nox, dchem_mix = dchem, 
                             dentrain_mix = tmp_dentrain)

    # ------------------------- END OF MIXING cases ----------------------------
    # --------------------------------------------------------------------------


    # --------------------------------------------------------------------------
    # update CO, CO2 and NOx concentration without chemical transformation
    # update CO and CO2 enhancements due to FF emissions
    p_update = p_all %>% 
               mutate(p_co_mix = p_co_mix + demiss_co / num_chunk, 
                      p_co2_mix = p_co2_mix + demiss_co2 / num_chunk, 
                      p_co_nomix = p_co_nomix + demiss_co / num_chunk, 
                      p_co2_nomix = p_co2_nomix + demiss_co2 / num_chunk, 
                    
                      # add a column for excluding chem changes, 10/25/2021
                      # treat NOx as a passive tracer like CO2 for debugging
                      p_nox_nochem = p_nox_nochem + demiss_no / num_chunk + dentrain_nochem)

    return(p_update)
}



# ------------------------------------
lookup_nox_tendency = function(p_df, tau_nox_df, OX_df, delt, tau_chem) {

    # get NOx NET frequency in hr based on [NOx] bin first
    if ( !is.numeric(tau_chem) ) {    # if using curves for lifetimes

        nox_bin = sort(unique(tau_nox_df$bin_nox))  # nox bin in ppb
        sza_bin = sort(unique(OX_df$bin_sza))

        # find corresponding NOx bins based on concentrations
        p_tmp = p_df %>% mutate(nox_indx = findInterval(tmp_nox * 1E3, nox_bin),
                                nox_indx = ifelse(nox_indx == 0, 1, nox_indx), 
                                gnox = nox_bin[nox_indx]) %>% 

                # merge NOx look-up tables by SZA and concentration bins
                left_join(tau_nox_df, by = c('gsza' = 'bin_sza', 
                                             'gnox' = 'bin_nox')) %>% 
                left_join(OX_df, by = c('gsza' = 'bin_sza')) %>% 
                dplyr::select(-c(gnox, nox_indx, n_nox))
    
    } else p_tmp = p_df %>% mutate(ts_nox = tau_chem)  # using constant lifetime

    return(p_tmp)
}

