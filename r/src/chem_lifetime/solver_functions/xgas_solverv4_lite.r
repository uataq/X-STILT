
### ----------------------------------------------------------------------
# lite version of xgas_solverv4() for debugging, DW, 07/26/2023

if (F) {

    X = 151
    traj_fn = traj_info$fn[X]
    timestr = traj_info$time[X]
    #aq_invent = ghg_invent = 'vulcan'; epa_fn = NA
    aq_invent = ghg_invent = 'epa'
    eno_sf = eco2_sf = eco_sf = 1
    ts_chem = NA
    perturb_emissTF = perturb_tsTF = perturb_mixTF = F
    perturb_indx = perturb_fn = NA
    xmin = xmax = ymin = ymax = NA 
    bg_nox = NA 

}

xgas_solverv4_lite = function(site, traj_fn, timestr, bg_nox = NA, eno_fn, 
                              eco_fn, eco2_fn, eco2_fn2 = NA, eno_sf = 1, 
                              eco_sf = 1, eco2_sf = 1, perturb_emissTF = F, 
                              perturb_tsTF = F, perturb_mixTF = F, 
                              perturb_fn = NA, perturb_indx = NA, tno2_fn, 
                              tno2_aux_path, xco_fn = NA, xco2_fn = NA, 
                              epa_fn = NA, epa_name = NA, epa_tz = 'PDT',
                              aq_invent = c('edgar', 'epa','odiac','vulcan')[1],
                              ghg_invent = c('edgar','epa','odiac','vulcan')[1],
                              mx_res = NA, ts_fn = NA, ts_chem = NA, xstilt_wd = '/central/home/dienwu/models/X-STILT'){

    start_time = Sys.time()
    setwd(xstilt_wd); source('r/dependencies.r')    # source all functions
    cat(traj_fn)
    
    ### ------------------------------------------------------------------------
    # 1. call functions for loading emissions and preparing particles
    # it will also calc dry-air column density and observed Mixing ratio
    # Rprofiling ----- this chunk takes ~20 sec at most, DW, 11/29/2022
    cat('\n\n# ---- xgas_solverv4_lite(): Loading emissions and trajectories ---- #\n')
    try_error = try(readRDS(traj_fn), silent = T)
    if ( 'particle' %in% names(try_error) & 'qt_prof' %in% names(try_error) ) {
        output = try_error
    } else if (attributes(try_error)$condition$message == 
               'error reading from connection') {
        cat('xgas_solverv4_lite(): trajec file is somehow broken, error reading it...\n\n')
        return()
    } else {
        cat('xgas_solverv4_lite(): printing error messages ---\n')
        print(try_error)
    }   # end if

    #print(str(output$qt_prof))
    receptor = output$receptor
    recp_lon = strsplit.to.df(basename(traj_fn))$V2
    recp_lat = strsplit.to.df(basename(traj_fn))$V3
    p = output$particle %>% filter(!is.na(lati), !is.na(long)) %>%
        mutate(recp_lon = recp_lon, recp_lat = recp_lat, 
               date = receptor$run_time + time * 60) 
    
    prp = p_emiss_prep(site, output, p, timestr, aq_invent, ghg_invent, eno_fn,
                       eco_fn, eco2_fn, eco2_fn2, epa_fn, epa_name, epa_tz, 
                       tno2_fn, xco_fn, xco2_fn, eno_sf, eco_sf, eco2_sf, 
                       perturb_emissTF, perturb_indx, perturb_fn) 

    if (is.null(prp)) return()
    pex = prp$pex %>% arrange(abs(time), indx) 
    tno2_info = prp$no2_info
    xco_info  = prp$co_info
    xco2_info = prp$co2_info 


    ### ------------------------------------------------------------------------
    # 2. grab NO2 and NOx net loss timescale [hr] and
    #    grab aux NO2 prior profile if needed, DW, 2022/01/13 
    # Rprofiling ----- this chunk takes 20 sec at most, DW, 11/29/2022
    cat('\n\n# ---- xgas_solverv4_lite(): Loading NOx lifetime curves and [NO2] initial conditions from TM5 ---- #\n')
    prp2 = p_ts_icbc_prep(pex, ts_chem, ts_fn, perturb_fn, perturb_tsTF, 
                          perturb_indx, bg_nox, tno2_aux_path, xstilt_wd)
    pex = prp2$pex 
    ts_nox_df = prp2$ts_nox_df
    min_nox = ifelse( is.numeric(ts_chem), 1e-6, min(ts_nox_df$bin_nox) )
    
    # if multiple mixing length scales and time scales are given
    if ( perturb_mixTF & !is.na(perturb_fn) ) {        
        mx_df = read.csv(perturb_fn) %>% filter(n == perturb_indx)
        mx_res = mx_df$mx_res 
        mx_hr = mx_df$mx_hr
    }   # end if

    p_chem = p_solverv4(pex, mx_res, mx_hr, bg_nox, ts_nox_df, ts_chem, min_nox)

    ### -------------------------------------------------------------
    # 4. perform vertical AK PWF weighting to particle-level concentrations
    cat('Performing vertical weighting...\n')
    p_recp = x_weightingv4(p_chem, tno2_info, xco_info, xco2_info)
    
}   # end of function


if (F) {


    ### ------------------------------------------------------------------------
    # 3. loop over each backward time and simulate concentrations
    # Rprofiling ----- this chunk 6 to 20 mins depending on dt, DW, 11/29/2022
    cat('\n\n# --- xgas_solverv4_lite(): Calculating [C] per timestamp --- #\n')
    x_df = NULL
    dffs = c(10, 50, 100, 500, 1500, 3000, 5000, 1e4, 5e4, 1e5) # m2 s-1
    for (mx_res in c(1, 3, 10)) {
        print(mx_res)
        
        for (dff in dffs) {
            cat('\n\nStarting a new diffusivity...\n'); print(dff); print(mx_hr)
            
            # mixing time scale in hr, calc as ts_mix = mx_res^2 / Dh
            # initial value of 3 hr for v1
            mx_hr = (mx_res * 1e3) ^2 / dff / 3600      
            p_chem = p_solverv4(pex, mx_res, mx_hr, bg_nox, ts_nox_df, 
                                ts_chem, min_nox)

            ### -------------------------------------------------------------
            # 4. perform vertical AK PWF weighting to particle-level concentrations
            cat('Performing vertical weighting...\n')
            p_recp = x_weightingv4(p_chem, tno2_info, xco_info, xco2_info)
            tmp_df = data.frame(
                    tno2_mix = mean(p_recp$p_no2_mix_wgt, na.rm = T) * 1e3, 
                    tno2_nomix = mean(p_recp$p_no2_nomix_wgt, na.rm = T) *1e3,
                    tno2_nochem = mean(p_recp$p_no2_nochem_wgt, na.rm = T) *1E3,
                    tnox_mix = mean(p_recp$p_nox_mix_wgt, na.rm = T) * 1E3, 
                    tnox_nomix = mean(p_recp$p_nox_nomix_wgt, na.rm = T) * 1E3, 
                    tnox_nochem = mean(p_recp$p_nox_nochem_wgt, na.rm = T) *1E3,
                    xco2_mix = mean(p_recp$p_co2_mix_wgt, na.rm = T), 
                    xco2_nomix = mean(p_recp$p_co2_nomix_wgt, na.rm = T)
                    ) %>% mutate(mx_res = mx_res, mx_hr = mx_hr, dff = dff)
            x_df = rbind(x_df, tmp_df)
        }   # end for
    }   # end for

    write.csv(x_df, file = 'test_tno2_tnox_mixing_202006151842_-89.16393_36.8939.csv', quote = F)
    
    csv_fns = list.files('./', 'test_tno2_tnox', full.names = T)
    x_df = do.call(rbind, lapply(csv_fns, function(x) {
        read.csv(x) %>% mutate(lon = as.numeric(strsplit.to.df(x)$V6), 
                               lat = as.numeric(strsplit.to.df(gsub('.csv', '', x))$V7)) }))

    #sf = 0.65
    plt_df = x_df %>% dplyr::select(c('lon', 'lat', 'tno2_mix', 'tno2_nomix', 
                                      'tnox_mix', 'tnox_nomix', 'mx_res', 
                                      'mx_hr', 'dff')) %>% 
             mutate(norm_frac = tno2_mix / tno2_nomix, 
                    dff_frac = (tno2_mix - tno2_nomix) / tno2_nomix) %>%
             filter(dff >= 50, dff <= 1e4) 

    # pp = pex %>% filter(ppTF) %>% group_by(indx) %>% add_count() %>% 
    #      ungroup() #%>% filter(n >= 1e4 / 10 / 60)
    # frac = sort(length(unique(pp$indx))) / 2000
    # if (!'ppTF' %in% colnames(pp)) frac = 0 
    lab_df = data.frame(lon = unique(plt_df$lon), lat = unique(plt_df$lat), 
                        y = rep(-0.04, 3), x = rep(50, 3), 
                        lab = c('17%', '61%', '0%'), mx_res = NA)

    m0 = ggplot(data = plt_df, aes(mx_hr, dff_frac, color = as.factor(mx_res)))+
         theme_bw() + geom_hline(yintercept = 0, linetype = 2) + 
         labs(x = 'Mixing e-folding timescale [hr]', 
              y = 'Normalized difference [unitless]', 
              title = 'Sensitivity test of mixing impact on modeled tNO2 over three receptors around New Madrid Power Plant on Jun 15th', 
              subtitle = 'Normalized tNO2 difference: (MIX - NOMIX) / NOMIX [ppb ppb-1]') +
         scale_x_log10(breaks = 10^seq(-3, 3), labels = 10^seq(-3, 3)) + 
         facet_wrap(~paste0(-lon,'W, ', lat,'N'), nrow = 1) +
         scale_shape_manual(name = 'SIM', values = c(19, 1)) +
         scale_linetype_manual(name = 'SIM', values = c(2, 1)) + 
         scale_color_manual(name = 'Mixing Box [km] ', 
                            values = RColorBrewer::brewer.pal(5, 'Dark2')) + 
         theme(legend.position = 'bottom', legend.key.width = unit(1.2, 'cm'), 
               legend.key.height = unit(0.3, 'cm')) +
         geom_point(size = 2) + geom_path() + 
         geom_text(data = lab_df, aes(x, y, label = paste(lab, 'trajec encountered\npower plant emission')))
    ggsave(m0, filename = 'FigS_mix.png', width = 13, height = 4.2, bg = 'white')


}