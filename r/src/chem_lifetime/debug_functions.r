
grab_p = function(rds_list, fac_list) {
       
       library(lubridate)
       pp = pp_avg = NULL
       for (l in 1 : length(rds_list)) {
              
              sel_fn = rds_list[[l]]
              if (length(sel_fn) == 0) next

              dat = readRDS(sel_fn)
              p_chem = dat$p_chem
              p_recp = dat$p_recp  

              emiss_max = as.numeric(quantile(p_chem$eno, 0.99))
              print(emiss_max)

              count_indx = length(unique((p_chem %>% filter(eno >= emiss_max))$indx))
              frac_indx = count_indx / max(p_chem$indx) * 100

              ### --------------------------------------------------------------
              tmp_pp = p_chem %>% #filter(indx == 1) %>% 
                     mutate(p_nox_mix = p_nox_mix * 1e3, 
                            p_nox_nomix = p_nox_nomix * 1e3, 
                            p_nox_nochem = p_nox_nochem * 1e3, 

                            p_co_mix = p_co_mix * 1E3, 
                            crnxc2 = p_nox_mix / p_co2_mix, 
                            crnxc2 = ifelse(crnxc2 == Inf, NA, crnxc2),
                            crn2c2 = p_no2_mix / p_co2_mix, 
                            crn2c2 = ifelse(crn2c2 == Inf, NA, crn2c2),
                            ern = eno / eco2 * 1E3, 
                            ern = ifelse(ern == Inf, NA, ern), 
                            crcc2 = p_co_mix / p_co2_mix, 
                            crcc2 = ifelse(crcc2 == Inf, NA, crcc2), 
                            tsf = air_vcd_stilt / dat$no2_info$air_vcd_tropo, 
                            xsf = air_vcd_stilt / dat$no2_info$air_vcd, 
                            fac = fac_list[[l]], dayTF = psza < 90) %>% 
                     dplyr::select(fac, time, lati, long, dayTF, xhgt, date, 
                                   temz, pres, p_rto_mix, p_rto_nomix, 
                                   p_nox_mix, p_no2_mix, p_no2_nomix, 
                                   p_nox_nomix, p_nox_nochem, ts_nox_mix, 
                                   ts_nox_nomix, eno, contains('demiss_no'), 
                                   ern, mlht, foot, 
                                   contains('p_co2_mix'), 
                                   contains('p_co2_nomix'), crnxc2, 
                                   crn2c2, crcc2, psza, xsf, tsf)
                            
              tmp_avg = tmp_pp %>% group_by(fac, time) %>% 
                        summarise_all(mean) %>% 
                        mutate(lt = with_tz(date, tzone = lon_lat$tz), 
                               frac_indx = frac_indx) %>% as.data.frame()
              
              tmp_avg$emiss_t = min((tmp_pp %>% filter(eno >= as.numeric(quantile(tmp_pp$eno, 0.99))))$time)

              pp = rbind(pp, tmp_pp)
              pp_avg = rbind(pp_avg, tmp_avg)
       }      # end for l

    return(list(pp = pp, pp_avg = pp_avg))
}


# for plotting trajec-level info
plot_pchem = function(pp) {

       emiss_sec = min(unique(pp$emiss_t) * 60)         # now in second
       frac_indx = max(unique(pp$frac_indx))
       min_day = min(pp[pp$dayTF, 'date'])

       r0 = ggplot(data = pp, aes(x = date, colour = fac)) + theme_bw() + 
            scale_color_manual(name = NULL, values = c('blue', 'purple', 
                                                       'orange', 'red')) + 
            scale_linetype_manual(name = NULL, values = c(1, 2)) +
            labs(y = NULL, x = NULL) + theme(legend.position = 'bottom', 
                                        legend.key.width = unit(1.2, 'cm')) + 
            geom_vline(xintercept = c(min_day, max(pp$date) + emiss_sec), 
                       linetype = 3, color = 'brown') +
            guides(linetype = guide_legend(nrow = 2), 
                   color = guide_legend(nrow = 1))  

       r1 = r0 + ggtitle('NO2-NOx ratio') + 
              geom_path(aes(y = p_rto_mix, linetype = 'mixing')) + 
              geom_path(aes(y = p_rto_nomix, linetype = 'no mixing')) 

       r2 = r0 + ggtitle('[NOx] in ppb') + 
              geom_path(aes(y = p_nox_mix, linetype = 'chem + mixing')) + 
              geom_path(aes(y = p_nox_nomix, linetype = 'chem + no mixing')) 

       r3 = r0 + ggtitle('[NO2] in ppb') + 
              geom_path(aes(y = p_no2_mix, linetype = 'chem + mixing')) + 
              geom_path(aes(y = p_no2_nomix, linetype = 'chem + no mixing')) 

       r4 = r0 + ggtitle('NOx "net loss" timescale [hr]')+
              geom_hline(yintercept = 4, linetype = 2, color = 'gray50') + 
              geom_path(aes(y = ts_nox_mix, linetype = 'mixing')) + 
              geom_path(aes(y = ts_nox_nomix, linetype = 'no mixing')) 

       r6 = r0 + geom_path(aes(y = eno)) + ggtitle('ENOx in umol m-2 s-1')
            annotate('text', mean(pp$date), max(pp$eno) * 0.85, 
                     label = paste(frac_indx, '% of trajec\naffected by max emiss')) 

       if ('p_co2_mix' %in% colnames(pp) & 'p_co2_nomix' %in% colnames(pp)) {
              r7 = r0 + ggtitle('[CO2] in ppm') + 
                   geom_path(aes(y = p_co2_mix, linetype = 'mixing')) +
                   geom_path(aes(y = p_co2_nomix, linetype = 'non mixing'))

       } else r7 = r0 + ggtitle('footprint [ppm /flux]') + 
                        geom_path(aes(y = foot)) 

       r5 = ggplot(data = pp, aes(long, lati, group = time, color = psza)) + 
            theme_bw() + geom_point(size = 0.2) + 
            scale_color_gradientn(colors = rev(def.col())) + 
            geom_point(data = pp %>% filter(time == -1), 
                       aes(long, lati), shape = 8, size = 5) + 
            ggtitle('Mean trajectories') + theme(legend.position = 'bottom') 

       r8 = r0 + geom_path(aes(y = psza)) + ggtitle('SZA')       
       rr = ggarrange(r6, r7, r8, r1, r2, r3, r4, r5, 
                      ncol = 4, nrow = 2, common.legend = T)

       return(rr)
}




if (F) {

r5 = r0 + geom_path(aes(y = ern)) + ylim(c(0, 20)) + 
          ggtitle('Flux ratio ENOx / ECO2 [mmol / mol]') 
r6 = r0 + geom_path(aes(y = mlht)) + ggtitle('Mixed layer hgt [m]') 
r7 = r0 + geom_path(aes(y = foot))+ ggtitle('Footprint [ppm / (umol m-2 s-1)]')+
          geom_hline(yintercept = mean(pp$foot), linetype = 2, color = 'gray50')
r9 = r0 + geom_path(aes(y = crnxc2)) + ylim(c(0, 20)) +
          ggtitle('ppb-[NOx] / ppm-[CO2]')
r10 = r0 + geom_path(aes(y = crn2c2)) + ylim(c(0, 20)) +
           ggtitle('ppb-[NO2] / ppm-[CO2]')
r11 = r0 + geom_path(aes(y = crcc2)) + ggtitle('ppb-[CO] / ppm-[CO2]')


}