
# grab AK from different sensors/species
# approx the dry-air vertical column density [mol m-2] based on surface pressure

# for final NOx, need to extract by the assumed background
# for total column estimates, need to have either TROPOMI XCO or OCO XCO2 data
x_weighting = function(p_chem, no2_info, co_info = NULL, co2_info = NULL) {

    ### ------------------------------------------------------------------------
    cat('x_weighting(): column weighting...\n')
    xdry_stilt = no2_info$air_vcd_stilt
    xdry_tropo = no2_info$air_vcd_tropo

    # only select the particles at the receptor time
    pr0 = p_chem %>% filter(time == max(time)) 

    # vertical weighting using AK and dilute MR for stilt levels
    # to tropospheric mean MR using xdry_tropo` and `xdry_stilt`,
    # assuming 0 FF enhancement above model level
    # add NOx column, DW, 2021/12/16
    pr_tno2 = pr0 %>%
              mutate(
                p_nox_mix_wgt = p_nox_mix * ak_tno2 * xdry_stilt / xdry_tropo,

                p_nox_nomix_wgt = p_nox_nomix * ak_tno2 * xdry_stilt / xdry_tropo,

                p_nox_nochem_wgt = p_nox_nochem * ak_tno2 * xdry_stilt / xdry_tropo,

                p_no2_nochem_wgt = p_nox_nochem * 0.74 * ak_tno2 * xdry_stilt / xdry_tropo, 

                p_no2_mix_wgt = p_no2_mix * ak_tno2 * xdry_stilt / xdry_tropo, 

                p_no2_nomix_wgt = p_no2_nomix * ak_tno2 * xdry_stilt / xdry_tropo ) 

    ### ------------------------------------------------------------------------
    # # if neither exists, use modeled dry air columns 
    # g = 9.8        	      	# m s-2
    # Mdry = 29 / 1E3        	# kg mol-1
    # psfc = max(pr0$pres)
    # # according to Odell 2012, Xdry = integrate( (1 - q ) / g / Mdry * dp )
    # # assuming q comes mainly from PBL to approx total H2ov column
    # xdrym = (1 - mean(pr0$sphu)) / g / Mdry * psfc * 100  # in mol m2
    
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
    

    # same for OCO-2/3
    if (!is.null(co2_info)) {
        xdry = co2_info$air_vcd_total
        pr = pr_xco %>% 
             mutate(p_co2_mix_wgt = p_co2_mix * ak_xco2 * xdry_stilt / xdry, 
                    p_co2_nomix_wgt = p_co2_nomix * ak_xco2 * xdry_stilt / xdry)
    } else pr = pr_xco %>% mutate(p_co2_mix_wgt = NA, p_co2_nomix_wgt = NA)
    
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
    cat(paste('XCO2 (+ mixing):', 
              signif(mean(pr$p_co2_mix_wgt, na.rm = T), 3), 'ppm\n'))
    cat(paste('XCO2 (- mixing):', 
              signif(mean(pr$p_co2_nomix_wgt, na.rm = T), 3), 'ppm\n'))

    return(pr)
}