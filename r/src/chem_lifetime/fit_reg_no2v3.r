
fit_reg_tno2v3 = function(pr_edgar = NULL, pr_epa = NULL, pr_vulcan = NULL, 
                          pr_odiac = NULL) {
    
    facs = c('01 - EDGAR (no chem, mix)',  
             '02 - EDGAR (chem, mix)', 
             '03 - EDGAR (chem, non-mix)',    
             '04 - EDGAR_FF (chem, mix)',  

             '05 - EPA hr (no chem, mix)', 
             '06 - EPA hr (chem, mix)', 
             '07 - EPA hr (chem, non-mix)',    
             '08 - EPA hr_FF (chem, mix)', 

             '09 - VULCAN (no chem, mix)',  
             '10 - VULCAN (chem, mix)', 
             '11 - VULCAN (chem, non-mix)', 

             '12 - ODIAC (no chem, mix)',  
             '13 - ODIAC (chem, mix)', 
             '14 - ODIAC (chem, non-mix)' )

    library(lmodel2)
    lm_df = NULL
    if ( !is.null(pr_edgar)) {
        lm1 = lmodel2(sim_tno2_nochem ~ obs_tno2, pr_edgar)$regression.results 
        lm2 = lmodel2(sim_tno2_mix ~ obs_tno2, pr_edgar)$regression.results 
        lm3 = lmodel2(sim_tno2_nomix ~ obs_tno2, pr_edgar)$regression.results 
        lm_df = rbind(lm_df, lm1 %>% mutate(fac = facs[1]), 
                             lm2 %>% mutate(fac = facs[2]), 
                             lm3 %>% mutate(fac = facs[3])) 
        
        if ('sim_tno2_ff' %in% colnames(pr_edgar)) {
            lm4 = lmodel2(sim_tno2_ff ~ obs_tno2_ff, pr_edgar)$regression.results 
            #lm4 = lmodel2(sim_tno2ff_mix ~ obs_tno2_ff, pr_edgar)$regression.results 
            lm_df = rbind(lm_df, lm4 %>% mutate(fac = facs[4]))
        }
    }

    if ( !is.null(pr_epa) ) {
        lm5 = lmodel2(sim_tno2_nochem ~ obs_tno2, pr_epa)$regression.results 
        lm6 = lmodel2(sim_tno2_mix ~ obs_tno2, pr_epa)$regression.results 
        lm7 = lmodel2(sim_tno2_nomix ~ obs_tno2, pr_epa)$regression.results 
        lm_df = rbind(lm_df, lm5 %>% mutate(fac = facs[5]), 
                             lm6 %>% mutate(fac = facs[6]),
                             lm7 %>% mutate(fac = facs[7]))
        if ('sim_tno2_ff' %in% colnames(pr_epa)) {
             lm8 = lmodel2(sim_tno2_ff ~ obs_tno2_ff, pr_epa)$regression.results 
            #lm8 = lmodel2(sim_tno2ff_mix ~ obs_tno2_ff, pr_epa)$regression.results 
            lm_df = rbind(lm_df, lm8 %>% mutate(fac = facs[8]))
        }
       
    }

    if ( !is.null(pr_vulcan) ) {
        lm9 = lmodel2(sim_tno2_nochem ~ obs_tno2, pr_vulcan)$regression.results 
        lm10 = lmodel2(sim_tno2_mix ~ obs_tno2, pr_vulcan)$regression.results 
        lm11 = lmodel2(sim_tno2_nomix ~ obs_tno2, pr_vulcan)$regression.results 
        lm_df = rbind(lm_df, lm9 %>% mutate(fac = facs[9]), 
                             lm10 %>% mutate(fac = facs[10]),
                             lm11 %>% mutate(fac = facs[11]))
    }

    if ( !is.null(pr_odiac) ) {
        lm12 = lmodel2(sim_tno2_nochem ~ obs_tno2, pr_odiac)$regression.results 
        lm13 = lmodel2(sim_tno2_mix ~ obs_tno2, pr_odiac)$regression.results 
        lm14 = lmodel2(sim_tno2_nomix ~ obs_tno2, pr_odiac)$regression.results 
        lm_df = rbind(lm_df, lm12 %>% mutate(fac = facs[12]), 
                             lm13 %>% mutate(fac = facs[13]),
                             lm14 %>% mutate(fac = facs[14]))
    }

    lm_sma = lm_df %>% filter(Method == 'SMA') %>% arrange(fac) %>% 
             mutate(sign = ifelse(Intercept < 0, '-', '+'), 
                    row = as.numeric(substr(fac, 1, 2)), 
                    text = paste0('sim', row, ' = ', signif(Slope, 3),' * obs ',
                                  sign, abs(signif(Intercept, 2)))) 
    return(lm_sma)
}