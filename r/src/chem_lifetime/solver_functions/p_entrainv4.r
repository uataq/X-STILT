
# add vertical dilution - if mixed layer height increases over time, 
# add extra term of (Ci_aloft - Ci) / H(t) * (dH / dt), DW, 08/18/2021

p_entrainv4 = function(p_curr, delt, dilutionTF = TRUE, min_nox = 1E-5) {
    
    library(dplyr)
    # preserve the p_* from current time stamp and grab MLH for the last time
    p_ent = p_curr %>%  
            mutate(
                dh = mlht - mlht_last,            # change in MLH in meters
                dh = ifelse(is.na(dh), 0, dh),
                dh_dt = dh / delt / 60,           # dH/dt [m /s]

                # only when dilutionTF is turned on and MLH increases
                # then mix [NOx] with air aloft
                # !!! only perform entrainment for particles within the 
                # mixed layer, where foot > 0 at the current time stamp
                dentrain_mix = ifelse(dh_dt > 0 & dilutionTF & foot > 0,   
                                      (min_nox - p_nox_mix) / mlht * dh_dt, 0), 
                dentrain_nomix = ifelse(dh_dt > 0 & dilutionTF & foot > 0,
                                     (min_nox - p_nox_nomix) / mlht * dh_dt, 0),
                dentrain_nochem = ifelse(dh_dt > 0 & dilutionTF & foot > 0, 
                                     (min_nox - p_nox_nochem) / mlht * dh_dt, 0)
            ) %>% dplyr::select(-c('mlht_last', 'dh_dt'))   
    return(p_ent)
}