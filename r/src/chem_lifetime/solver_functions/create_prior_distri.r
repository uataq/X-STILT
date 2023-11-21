
# requires ENO from at least two emission inventories including one from EDGAR
create_prior_distri = function(site, timestr, eno_edgar_fn, eco2_edgar_fn, 
                               eco2_odiac_fn, inv_path, lon_lat, n_perturb, 
                               overwriteTF = F) {

    perturb_fn = file.path(inv_path, 
                           paste0('perturb_emiss_norm_', site, '.rds'))
    
    if ( !file.exists(perturb_fn) | overwriteTF ) {   
        err_stk = load_eno_uncert_norm(site, timestr, eno_edgar_fn, 
                                       eco2_edgar_fn, eco2_odiac_fn, 
                                       aq_invent = 'edgar', 
                                       xmin = lon_lat$minlon, 
                                       xmax = lon_lat$maxlon,
                                       ymin = lon_lat$minlat, 
                                       ymax = lon_lat$maxlat, plotTF = T)

        err_df = err_stk %>% as.data.frame(xy = T) %>% arrange(x, y) %>% 
                 tibble::rownames_to_column(var = 'indx') %>% 
                 mutate(indx = as.numeric(indx), 
                        CV_ENO = ifelse(CV_ENO > 10, 10, CV_ENO))
        
        # work on each emission grid and minimize the diff in emissions
        sf_df = NULL
        for (ii in err_df$indx) {
            
            print(ii)
            tmp_df = err_df %>% filter(indx == ii)

            # generate scaling factor for emission perturbation assuming no bias
            # minimize the mismatch in mean between perturbation and EDGAR
            # minimize the number of negative SF. 
            check_mu = check_sd = 1  # initialization
            check_ng = n_perturb
            
            while ( check_mu > 0.01 | check_sd > 0.05 | 
                    check_ng >= ifelse(tmp_df$CV_ENO < 1, 2, n_perturb / 2) ) {
                loop_df = data.frame(x = rep(tmp_df$x, n_perturb), 
                                     y = rep(tmp_df$y, n_perturb),
                                     ens = 1: n_perturb, 
                                     mu = tmp_df$EDGAR_ENO, 
                                     cv = tmp_df$CV_ENO) %>%
                        
                        # prior stat in normal distribution
                        mutate(sf = rnorm(n_perturb, mean = 1, sd = cv))
                check_mu = abs(mean(loop_df$sf) - 1)
                check_sd = abs(sd(loop_df$sf) - unique(loop_df$cv))
                check_ng = length(loop_df$sf[loop_df$sf < 0])
            }   # end while
            #hist(loop_df$sf)

            sf_df = rbind(sf_df, loop_df %>% mutate(emis_indx = ii))
        }   # end for

        # SF for emission scaling factor, 
        # emis for perturbed emissions, 
        # mu for original inventory emissions, 
        # indx for emission grid cell
        sf_df = sf_df %>% mutate(sf = ifelse(sf < 0, 1e-6, sf), emis = sf * mu) 
        
        saveRDS(sf_df, file = perturb_fn)
    } else sf_df = readRDS(perturb_fn)

    #sfm = sf_df %>% group_by(x, y) %>% summarise_all(mean)# %>% filter(sf < 2)
    #s1 = ggplot(sfm) + geom_raster(aes(x, y, fill = emis - mu))

    return(sf_df)
}