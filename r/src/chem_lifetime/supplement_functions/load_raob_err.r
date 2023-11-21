
load_raob_err = function(site, met, err_path) {

    raob_df = data.frame(site = c('Flagstaff', 'Yuma', 'Phoenix', 'Tucson', 
                                  'SaltLakeCity', 'Dugway', 'Springfield', 
                                  'Lincoln'), 
                         lat = c(32.23, 32.50, 33.45, 35.23, 40.77, 
                                 40.18, 37.23, 40.15), 
                         lon = c(-110.96, -114.00, -111.95, -111.82, -111.97,
                                 -112.92, -93.40, -89.33))

    err_fns = list.files(err_path, paste0(met, '_raob'), full.names = T)

    # load raob results
    if ( length(err_fns) > 0 ) {

        err_df = do.call(rbind, lapply(err_fns, function(x){
            #print(x)
            tmp_tt = gsub('.txt', '', strsplit.to.df(basename(x))$V4)
            readr::read_csv(x, show_col_types = FALSE) %>% 
            mutate(tropomi = tmp_tt)
        })) %>% left_join(raob_df, by = c('lat', 'lon')) 

        # no need to average over all vertical levels below 2km
        comb_df = err_df %>% mutate(site = 'ALL') 
        all_df = rbind(comb_df, err_df) %>% group_by(site, tropomi) %>% 
                  mutate(date = as.POSIXct(tropomi, 'UTC', format = '%Y%m%d%H'),
                         yr = substr(tropomi, 1, 4), 
                         siguverr = sqrt(mean(c(u.err^2, v.err^2))), 
                         bias_uv = mean(c(u.err, v.err)), 
                         rmse_ws = sqrt(mean(ws.err^2)), 
                         rmse_wd = sqrt(mean(wd.err^2)), 
                         bias_ws = mean(ws.err), 
                         bias_wd = mean(wd.err), 
                         avguv = mean(ws.met), frac = siguverr / avguv, 
                         met = met) %>% ungroup() %>% unique()
                  #dplyr::select(site, tropomi, met, frac, siguverr, 
                  #              contains('rmse_'), contains('bias_')) 
                  
    } else all_df = NULL

    return(all_df)
}