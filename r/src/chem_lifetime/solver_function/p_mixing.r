
# 4. horizontal mixing of [NOx] among particles, DW, 04/21/2021
# if they fall within the same footprint grid (1km)
# calculate mean [NOx] among particles falling within the same grid
# select particles within the mixed layer, if foot > 0 (zagl < mlht)

# always perform mixing and no-mixing, using two columns, DW, 11/09/2021
# allow for partial mixing, DW, 01/26/2022
p_mixing = function(p_tau, mx_hr, delt) {

    # add mixing of [NOx] among particles, DW, 04/21/2021
    p_tly = p_tau %>% group_by(glon, glat) %>% 
            tally() %>% ungroup() %>% rename(n_mix = n)
    
    p_mix = p_tau %>% group_by(glon, glat) %>% filter(foot > 0) %>% 
            dplyr::summarise(p_nox_mix_mean = mean(p_nox_mix, na.rm = T),
                             p_nox_nochem_mean = mean(p_nox_nochem, na.rm = T),
                             p_co2_mix_mean = mean(p_co2_mix, na.rm = T), 
                             p_co_mix_mean = mean(p_co_mix, na.rm = T), 
                             .groups = 'drop') %>% ungroup() %>% 
            left_join(p_tly, by = c('glon', 'glat'))


    # introduce mixing timescale
    # C = C(t) * exp(-delt/tmix) + C_mean(t) * [1 - exp(-delt/tmix)]
    # define diff rate, dff = exp(-delt/tmix)
    # delC = C * (dff - 1) + C_mean * (1 - dff)
    # DW, 01/26/2022
    dffs = exp(-delt / 60 / mx_hr)        # delt in min, mx_hr in hour
    p_fn2 = p_tau %>% left_join(p_mix, by = c('glon', 'glat')) %>% 

            # !!! only modify p_* due to mixing if its foot > 0, DW, 11/09/2021
            mutate(dmix_nox = ifelse(is.na(p_nox_mix_mean) | foot == 0, 0, 
                                     p_nox_mix * (dffs - 1) + p_nox_mix_mean * (1 - dffs)),

                   dmix_nox_nochem = ifelse(is.na(p_nox_nochem_mean)| foot == 0,
                                            0, p_nox_nochem * (dffs - 1) + p_nox_nochem_mean * (1 - dffs)), 

                   dmix_co2 = ifelse(is.na(p_co2_mix_mean) | foot == 0, 0, 
                                     p_co2_mix * (dffs - 1) + p_co2_mix_mean * (1 - dffs)),

                   dmix_co = ifelse(is.na(p_co_mix_mean) | foot == 0, 0, 
                                    p_co_mix * (dffs - 1) + p_co_mix_mean * (1 - dffs)),

                   n_mix = ifelse(is.na(n_mix) | foot == 0, 0, n_mix), 
                   
                   # now dmix_* is stored, only update p_* with mixing
                   p_co_mix = p_co_mix + dmix_co,
                   p_co2_mix = p_co2_mix + dmix_co2,
                   p_nox_mix = p_nox_mix + dmix_nox,
                   p_nox_nochem = p_nox_nochem + dmix_nox
                   ) %>% dplyr::select(-ends_with('_mean'))

    return(p_fn2)
}






if (F) {

sz = 2

p0 = ggplot() + theme_bw() + theme(legend.position = 'bottom') + 
     geom_tile(data = p_mix, aes(glon + 0.045, glat + 0.045, 
               fill = p_nox_mean * 1e3), color = 'gray', size = 0.5) + 
     scale_fill_gradient2(low = 'white', high = 'black', name = 'C_mean') + 
     xlim(c(-113.2, -112.8)) + ylim(c(39.9, 40.3))

p1 = p0 + ggtitle('C for each parcel in ppb before mixing') + 
     geom_point(data = p_tau %>% filter(foot > 0), 
                aes(long, lati, color = p_noxa * 1e3), size = sz) + 
     scale_color_gradient2(low = 'yellow', high = 'orange', name = 'C_par')
     
p2 = p0 + ggtitle('dC in ppb due to PARTIAL mixing in 1 min') + 
     geom_point(data = p_fn2 %>% filter(foot > 0), 
                aes(long, lati, color = dmix_nox * 1e3), size = sz) + 
     scale_color_gradient2(low = 'blue', high = 'orange', name = 'dC_mix')

p2a = p0 + ggtitle('C in ppb after PARTIAL mixing in 1 min') + 
     geom_point(data = p_fn2 %>% filter(foot > 0), 
                aes(long, lati, color = p_noxa * 1e3), size = sz) + 
     scale_color_gradient2(low = 'yellow', high = 'orange', name = 'C_par')

p3 = p0 + ggtitle('dC in ppb due to COMPLETE mixing in 1 min (wrong)') + 
     geom_point(data = p_fn1 %>% filter(foot > 0), 
                aes(long, lati, color = dmix_nox * 1e3), size = sz) + 
     scale_color_gradient2(low = 'blue', high = 'orange', name = 'dC_mix') 

p3a = p0 + ggtitle('C in ppb after COMPLETE mixing in 1 min (wrong)') + 
     geom_point(data = p_fn1 %>% filter(foot > 0), 
                aes(long, lati, color = p_noxa * 1e3), size = sz) + 
     scale_color_gradient2(low = 'yellow', high = 'orange', name = 'C_par')

# as a ref, plot dC due to chemistry
p4 = p0 + ggtitle('dC in ppb due to CHEM in 1 min') + 
     geom_point(data = p_fn2 %>% filter(foot > 0), 
                aes(long, lati, color = -dchem * 1e3), size = sz) + 
     scale_color_gradient2(low = 'blue', high = 'yellow', name = 'dC_chem') 

# double check lifetime
c1 = ggplot(data = p_fn2) + geom_point(aes(p_noxb, ts_nox_nomix)) + 
     scale_x_continuous(trans = 'log10')

pp = ggarrange(p1, p4, p2, p3, p2a, p3a, ncol = 2, nrow = 3)
ggsave(pp, filename = '/home/dienwu/postdoc_proj/NOx/mixing_t1_zoom.png', 
           width = 10, height = 15)


}