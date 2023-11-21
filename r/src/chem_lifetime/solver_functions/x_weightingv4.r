
# grab AK from different sensors/species
# approx the dry-air vertical column density [mol m-2] based on surface pressure

# for final NOx, need to extract by the assumed background
# for total column estimates, need to have either TROPOMI XCO or OCO XCO2 data
x_weightingv4 = function(p_chem, no2_info, co_info = NULL, ch4_info = NULL, 
                         co2_info = NULL, output = NULL) {

    ### ------------------------------------------------------------------------
    cat('x_weightingv4(): column weighting...\n')
    xdry_stilt = no2_info$air_vcd_stilt     # dry-air STILT column (e.g., 0-3km)
    xdry_tropo = no2_info$air_vcd_tropo     # dry-air tropo column
    
    # only select the particles at the receptor time
    pr0 = p_chem %>% filter(time == max(time)) 

    # vertical weighting using AK and dilute MR for stilt levels
    # to tropospheric mean MR using xdry_tropo` and `xdry_stilt`,
    # assuming 0 FF enhancement above model level
    # add NOx column, DW, 2021/12/16
    pr_tno2 = pr0 %>% mutate(
            p_nox_mix_wgt = p_nox_mix * ak_tno2 * xdry_stilt / xdry_tropo,
            p_nox_nomix_wgt = p_nox_nomix * ak_tno2 * xdry_stilt / xdry_tropo,
            p_nox_nochem_wgt = p_nox_nochem * ak_tno2 * xdry_stilt / xdry_tropo,
            #p_no2_nochem_wgt = p_nox_nochem * 0.74 * ak_tno2 * xdry_stilt /xdry_tropo,
            p_no2_nochem_wgt = p_nox_nochem_wgt * 0.74, 
            p_no2_mix_wgt = p_no2_mix * ak_tno2 * xdry_stilt / xdry_tropo, 
            p_no2_nomix_wgt = p_no2_nomix * ak_tno2 * xdry_stilt / xdry_tropo ) 
    
    ### ------------------------------------------------------------------------
    # vertical weighting using AK and convert MR for stilt levels
    # to total X for XCO using `xdry_tropo` and `xdry_stilt`, 
    # assuming zero enhancement above model levels
    if (!is.null(co_info)) {
        xdry = co_info$air_vcd_total
        pr_xco = pr_tno2 %>%   
                 mutate(p_co_mix_wgt = p_co_mix * ak_xco * xdry_stilt / xdry, 
                        p_co_nomix_wgt = p_co_nomix * ak_xco * xdry_stilt /xdry)

    } else pr_xco = pr_tno2 %>% mutate(p_co_mix_wgt = NA, p_co_nomix_wgt = NA)
    # if no TROPOMI XCO or AK CO is available, flag it as NA
    

    ### ------------------------------------------------------------------------
    # if for TROPOMI xCH4
    if (!is.null(ch4_info)) {
        xdry = ch4_info$air_vcd_total
        pr_xch4 = pr_xco %>%   
                  mutate(p_ch4_mix_wgt = p_ch4_mix * ak_xch4 * xdry_stilt /xdry,
                         p_ch4_nomix_wgt = p_ch4_nomix * ak_xch4 * xdry_stilt /xdry)

    } else pr_xch4 = pr_xco %>% mutate(p_ch4_mix_wgt = NA, p_ch4_nomix_wgt = NA)
    # if no TROPOMI XCO or AK CO is available, flag it as NA
    

    ### ------------------------------------------------------------------------
    # if no OCO info available, treat AK_xco2 as 1
    # use STILT modeled total column dry-air
    if ( is.null(co2_info) ) {
        xdry = unique(get.xdry.mod(output)$xdry.tot) 
        pr_xch4$ak_xco2 = 1
    } else xdry = co2_info$air_vcd_total
    
    pr = pr_xch4 %>% 
         mutate(p_co2_mix_wgt = p_co2_mix * ak_xco2 * xdry_stilt / xdry, 
                p_co2_nomix_wgt = p_co2_nomix * ak_xco2 * xdry_stilt / xdry)

    cat(paste('tNO2 (+ quadratic ratio + NOx chem + mixing):', 
              signif(mean(pr$p_no2_mix_wgt, na.rm = T) * 1E3, 3), 'ppb\n'))
    cat(paste('tNO2 (+ quadratic ratio + NOx chem - mixing):', 
              signif(mean(pr$p_no2_nomix_wgt, na.rm = T) * 1E3, 3), 'ppb\n'))
    cat(paste('tNO2 (+ const of 0.74 in Beirle2021 - NOx chem + mixing):', 
              signif(mean(pr$p_no2_nochem_wgt, na.rm = T) * 1E3, 3), 'ppb\n'))
    cat(paste('tNOx (+ NOx chem + mixing):', 
              signif(mean(pr$p_nox_mix_wgt, na.rm = T) * 1E3, 3), 'ppb\n'))
    cat(paste('tNOx (- NOx chem + mixing):', 
              signif(mean(pr$p_nox_nochem_wgt, na.rm = T) * 1E3, 3), 'ppb\n'))
    cat(paste('dXCO2 (+ mixing):', 
              signif(mean(pr$p_co2_mix_wgt, na.rm = T), 3), 'ppm\n'))
    cat(paste('dXCO (+ mixing):', 
              signif(mean(pr$p_co_mix_wgt, na.rm = T) * 1e3, 3), 'ppb\n'))
    cat(paste('dXCH4 (+ mixing):', 
              signif(mean(pr$p_ch4_mix_wgt, na.rm = T) * 1e3, 3), 'ppb\n'))

    return(pr)
}