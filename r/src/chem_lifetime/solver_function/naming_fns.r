
naming_fns = function(aq_invent, ghg_invent, mx_res = 1, eno_sf = 1, 
                      tau_chem = c(4, NA)[1], perturb_emiss_indx = NA, 
                      perturb_tau_indx = NA, perturb_mix_indx = NA) {
    
    fn = paste0('_pchem_', aq_invent, '_', ghg_invent, '.rds')
    
    if (!is.na(perturb_mix_indx)) {
        fn = gsub('.rds', paste0('_mix', perturb_mix_indx, '.rds'), fn)
        
    } else fn = paste0('_pchem_', aq_invent, '_', ghg_invent, 
                       '_', mx_res, 'mixing.rds')

    if (eno_sf != 1) fn = gsub('.rds', paste0('_sf', eno_sf, '.rds'), fn)

    if ( is.numeric(tau_chem) ) 
        fn = gsub('pchem', paste0('pchem_', tau_chem, 'hrs'), fn)

    if ( !is.na(perturb_emiss_indx) ) 
        fn = gsub('.rds', paste0('_emiss', perturb_emiss_indx, '.rds'), fn)
    
    if ( !is.na(perturb_tau_indx) ) 
        fn = gsub('.rds', paste0('_tau', perturb_tau_indx, '.rds'), fn)

    # if ( !is.na(tau_chem) ) {   # this is for debugging
    #     if ( tau_chem == 'v20230503' ) 
    #         fn = gsub('pchem', paste0('pchem_', tau_chem), fn)
    # }

    return(fn)
}

