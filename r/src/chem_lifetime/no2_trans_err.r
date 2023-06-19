
# rindx can be a vector if parallel computing is activated
calc_no2_trans_err = function(site, timestr, emiss, met, rindx = 1, 
                              out_init_path, out_err_path, out_fn = NULL, 
                              xstilt_wd) {
       
       setwd(xstilt_wd); source('r/dependencies.r')

       # ----------------------------------------------------------
       # load sim without or with transport errors
       fns1 = list.files(out_init_path, paste0(emiss, '_1mixing.rds'), 
                         recursive = T, full.names = T)
       fns2 = list.files(out_err_path, paste0(emiss, '_1mixing.rds'), 
                         recursive = T, full.names = T)

       info = strsplit.to.df(basename(fns1)) %>% 
              dplyr::select(timestr = V1, lon = V2, lat = V3) %>% 
              mutate_if(is.character, as.numeric)

       # work on one receptor --------------------------------
       tmp_info = info[rindx, ]
       tmp_fn1 = fns1[grepl(tmp_info$lat, fns1) & grepl(tmp_info$lon, fns1)]
       tmp_fn2 = fns2[grepl(tmp_info$lat, fns2) & grepl(tmp_info$lon, fns2)]
       if (length(tmp_fn1) * length(tmp_fn2) == 0) return()
       

       # load particle-level tropospheric NO2 -----------------------------
       pinfo = readRDS(tmp_fn1)$p_recp %>% 
               dplyr::select(indx, long, lati, pres, psza, xhgt)
       pr1 = readRDS(tmp_fn1)$p_chem %>% filter(time == max(time)) %>% 
             mutate(pno2_ini = p_nox_mix * p_rto_mix) %>%
             dplyr::select(indx, pno2_ini)
       pr2 = readRDS(tmp_fn2)$p_chem %>% filter(time == max(time)) %>% 
             mutate(pno2_err = p_nox_mix * p_rto_mix) %>%
             dplyr::select(indx, pno2_err)
       pr = left_join(pr1, pr2, by = 'indx') %>% left_join(pinfo, by = 'indx')


       # now load TROPOMI profile and AK -----------------------------
       trp_prof = as.data.frame(readRDS(tmp_fn1)$no2_info[c('lower_pres', 'ak_tropo')])
       tropause = min(trp_prof[trp_prof$ak_tropo > 0, 'lower_pres'])
       
       # instead of using TROPOMI levels, create more arbitrary levels 
       #p_lower = trp_prof$lower_pres[-1]

       # ----------------------------------------------------------
       # create a while loop to determine the proper pressure levels
       # initialization: 
       nbins = 6 
       npar_min = 1    # numbers of particles per levels
       nlev_pos = 1    # numbers of pres levels with positive diff in variances
       check_nlev_pos = 2
       check_npar_min = 50

       # ** sufficient # of levels and sufficient # of particles per level
       while ( nlev_pos <= check_nlev_pos | npar_min < check_npar_min ) {
              
              p_lower = seq(floor(min(pr$pres)), 
                            floor(max(pr$pres)), length = nbins)

              ak_int = approx(x = trp_prof$lower_pres, 
                              y = trp_prof$ak_tropo, 
                              xout = p_lower, rule = 2)$y
              int_prof = data.frame(lower_pres = p_lower, ak_tropo = ak_int)
              dp = diff(p_lower)

              # use the variance of two simulations as uncertainties
              pr_rev = pr %>% 
                       mutate(pres_indx = findInterval(pres, p_lower) + 1, 
                              pres_indx = ifelse(pres_indx > length(p_lower), 
                                                 length(p_lower), pres_indx), 
                              lower_pres = p_lower[pres_indx]) %>%
                       left_join(int_prof, by = 'lower_pres') %>% 

                       # convert to ppb
                       mutate(pno2_ini = pno2_ini * 1e3, 
                              pno2_err = pno2_err * 1e3)

              # for each layer/level, calc the difference in variance in NO2 
              var_df = NULL
              for (ll in unique(pr_rev$pres_indx)) {

                     tmp = pr_rev %>% filter(pres_indx == ll) 
                     tmp_p_err = sort(tmp$pno2_err)
                     tmp_p_ini = sort(tmp$pno2_ini)
                     n = length(tmp_p_ini)
                     tmp_var_df = data.frame(pres_indx = ll, n = n, 
                                             lower_pres = p_lower[ll], 
                                             var_err = sd(tmp_p_err)^2, 
                                             var_ini = sd(tmp_p_ini)^2, 
                                             avg_pno2 = mean(tmp$pno2_ini)) %>% 
                                  mutate(dvar = var_err - var_ini)
                     var_df = rbind(var_df, tmp_var_df)
              }      # end for

              npar_min = min(var_df$n)
              nlev_pos = nrow(var_df[var_df$dvar > 0, ])
              nbins = nbins + 1

              # to prevent infinite loop, 
              # reduce threshold for # of pos level as 1; 
              # if no positive levels for all loop, drop this receptor? 
              if (nbins > 20 & check_nlev_pos == 2) check_nlev_pos = 1 
              if (nbins > 20 & check_nlev_pos == 1) break 
       }      # end while
       

       # ----------------------------------------------------------
       # correct for negative diff in variances if there is
       if ( TRUE %in% (var_df$dvar < 0) ) {
              stat = var_df %>% rename(dVAR = dvar, var.err = var_err, 
                                   var.orig = var_ini)
              rev_info = scale.dvar(stat)
              fit_df = rev_info[[2]]
              rev_df = rev_info[[1]] %>% rename(dvar = scaled.dvar) %>% 
                       left_join(int_prof, by = 'lower_pres')
              
              # if linear reg intercept is negative, it creates negative dVAR, 
              # drop levels (assuming zero errors, since var.orig is small)
              if (fit_df[fit_df$validTF, 'i'] < 0) rev_df = rev_df %>% na.omit()

       } else rev_df = var_df %>% left_join(int_prof, by = 'lower_pres')


       # ----------------------------------------------------------
       # variance of column = sum(wgt^2 * variance), wgt includes PW and AK
       # since tNO2 = sum(wgt * PNO2), PWF = dp_level / dp_tropo
       pwf = dp[unique(pr_rev$pres_indx) - 1] / (max(p_lower) - tropause)
       tno2_sd = sqrt(sum(rev_df$dvar * rev_df$ak_tropo^2 * pwf^2))
       print(tno2_sd)

       # ----------------------------------------------------------
       # store all back
       tmp_df = tmp_info %>% mutate(trans_err_ppb = tno2_sd, 
                                    nlev_tot = nbins, 
                                    npar_min = npar_min, 
                                    nlev_pos = nlev_pos)
       
       if (is.null(out_fn)) 
              out_fn = file.path(dirname(out_init_path), 
                                 paste0('transerr_stat_', site, '_', timestr, 
                                        '_', emiss, '_', met, '.csv'))
       
       write.table(tmp_df, out_fn, sep = ',', row.names = F, 
                   col.names = !file.exists(out_fn), 
                   append = file.exists(out_fn))
       cat(paste('Done', rindx, '\n'))
}


load_trans_err = function(site, timestr, emiss, met, out_init_path, rds_path) {

       # calculate trans error
       out_fn = paste0('../error_inversion/tmp/transerr_stat_', site, 
                       '_', timestr, '_', emiss, '_', met, '.csv')

       err_df = read.csv(out_fn) %>% 
                rename(var.orig = var_ini, var.err = var_err) %>% 
                mutate(dVAR = var.err - var.orig)
       lr_stat_info = scale.dvar(err_df)
       err_df$sd_trans = lr_stat_info[[1]]$sd.trans * 1e3  # ppm to ppb

       # load initial results 
       out_init = obs_sim_pairv3(site, timestr, out_init_path, rds_path, 
                                 paste0(emiss, '_1mixing'), met, F)$x_df

       all_df = out_init %>% left_join(err_df %>% dplyr::select(-timestr),  
                                       by = c('recp_lon' = 'lon', 
                                              'recp_lat' = 'lat'))
       return(all_df)
}