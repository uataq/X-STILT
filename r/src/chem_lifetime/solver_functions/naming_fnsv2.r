
naming_fnsv2 = function(aq_invent, ghg_invent, mx_res = 1, eno_sf = 1, 
                        tau_chem = c(4, NA)[1], perturb_indx = NA, 
                        perturb_emiTF = F, perturb_tsTF = F, perturb_mixTF = F){
    
    fn = paste0('_pchem_', aq_invent, '_', ghg_invent, '.rds')
    
    if ( perturb_mixTF ) {
        fn = gsub('.rds', paste0('_mix', perturb_indx, '.rds'), fn)
        
    } else fn = paste0('_pchem_', aq_invent, '_', ghg_invent, 
                       '_', mx_res, 'mixing.rds')

    if (eno_sf != 1) fn = gsub('.rds', paste0('_sf', eno_sf, '.rds'), fn)

    if ( is.numeric(tau_chem) ) 
        fn = gsub('pchem', paste0('pchem_', tau_chem, 'hrs'), fn)

    if ( perturb_emiTF ) 
        fn = gsub('.rds', paste0('_emiss', perturb_indx, '.rds'), fn)

    if ( perturb_tsTF ) 
        fn = gsub('.rds', paste0('_tau', perturb_indx, '.rds'), fn)

    return(fn)
}

