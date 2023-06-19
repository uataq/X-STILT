
gapfill = function(df, uni_species, species = c('NOx', 'NO2')[1]) {


    if (species == 'NOx') {  # fix extremely high lifetime for low [NOx]
        fix_df = df %>% filter(bin_species <= 0.01)
        q0p95  = quantile(fix_df$ts_species, 0.95)

        if (q0p95 >= 5 & nrow(fix_df %>% filter(ts_species >= q0p95)) > 0) {
            fix_avg = mean((fix_df %>% filter(ts_species <= q0p95))$ts_species, 
                        na.rm = T)
            df = df %>% mutate(ts_species = ifelse(bin_species <= 0.01 & 
                                               ts_species >= q0p95, 
                                               fix_avg, ts_species)) 
        }
    }    # end if

    # fix few outliers, extremely high lifetime for high [NO2] during daytime
    if (species == 'NO2' & unique(df$bin_sza) < 90) {
        fix_avg = mean((df %>% filter(bin_species >= 40, ts_species < 15))$ts_species, na.rm = T)

        df = df %>% mutate(ts_species = ifelse(bin_species >= 40 & 
                                               ts_species > 15, 
                                               fix_avg, ts_species))
    }


    na = uni_species[!uni_species %in% df$bin_species]
    if ( length(na) > 0 ) {          # now gap fill

        # average timescales from the highest or lowest three [NOx] bins 
        nn = 3#; if (species == 'NO2') nn = 5
        high_df  = df[(nrow(df) - nn + 1) : nrow(df), ]
        high_avg = mean(high_df$ts_species)
        if (species == 'NOx' && high_avg < 72) high_avg = 72
        
        low_df  = df[1 : nn, ]
        low_avg = mean(low_df$ts_species)
        if (species == 'NOx' && low_avg > 3) low_avg = 3

        # select high or low [NOx] bins
        na_high = na[na >= min(high_df$bin_species)]
        na_low  = na[na <= min(high_df$bin_species)]

        if (length(na_high) > 0) {
            high_gf = data.frame(bin_species = na_high, bin_sza = x, 
                                 n = 0, ts_species = high_avg)
        } else high_gf = NULL
    
        if (length(na_low) > 0) {
            low_gf = data.frame(bin_species = na_low, bin_sza = x, 
                                n = 0, ts_species = low_avg)
        } else low_gf = NULL
        
        # merge low df, initial df, high_df with modified tau and NA df
        tmp_gf = rbind(low_gf, #df, 
                       df[1 : (nrow(df) - nn), ], 
                       high_df %>% mutate(ts_species = high_avg), 
                       high_gf)

    } else tmp_gf = df   # end if

    return(tmp_gf)
}




gapfillv2 = function(df, uni_species, species = c('NOx', 'NO2')[1]) {


    if (species == 'NOx') {  # fix extremely high lifetime for low [NOx]
        fix_df = df %>% filter(bin_species <= 0.01)
        q0p95  = quantile(fix_df$ts_species, 0.95)

        if (q0p95 >= 5 & nrow(fix_df %>% filter(ts_species >= q0p95)) > 0) {
            fix_avg = mean((fix_df %>% filter(ts_species <= q0p95))$ts_species, 
                        na.rm = T)
            df = df %>% mutate(ts_species = ifelse(bin_species <= 0.01 & 
                                               ts_species >= q0p95, 
                                               fix_avg, ts_species)) 
        }
    }    # end if

    # fix few outliers, extremely high lifetime for high [NO2] during daytime
    if (species == 'NO2' & unique(df$bin_sza) < 90) {
        fix_avg = mean((df %>% filter(bin_species >= 40, ts_species < 15))$ts_species, na.rm = T)

        df = df %>% mutate(ts_species = ifelse(bin_species >= 40 & 
                                               ts_species > 15, 
                                               fix_avg, ts_species))
    }


    na = uni_species[!uni_species %in% df$bin_species]
    if ( length(na) > 0 ) {          # now gap fill

        # linear interpolate timescales from highest or lowest three [NOx] bins 
        nn = 5
        high_df = df[(nrow(df) - nn + 1) : nrow(df), ]
        low_df  = df[1 : nn, ]

        # select high or low [NOx] bins
        na_high = na[na >= min(high_df$bin_species)]
        na_low  = na[na <= min(high_df$bin_species)]

        # predict ts for the two ends
        high_coef = coef(lm(ts_species ~ log(bin_species), data = high_df))
        high_pred = log(na_high) * high_coef[2] + high_coef[1]
        high_pred = ifelse(high_pred > 72, 72, high_pred)
        
        low_coef = coef(lm(ts_species ~ log(bin_species), data = low_df))
        low_pred = log(na_low) * low_coef[2] + low_coef[1]
        low_pred = ifelse(low_pred < -72, -72, low_pred)
    
        if (length(na_high) > 0) {
            high_gf = data.frame(bin_species = na_high, bin_sza = x, 
                                 n = 0, ts_species = high_pred)
        } else high_gf = NULL

        if (length(na_low) > 0) {
            low_gf = data.frame(bin_species = na_low, bin_sza = x, 
                                n = 0, ts_species = low_pred)
        } else low_gf = NULL
        
        # merge low df, initial df, high_df with modified tau and NA df
        tmp_gf = rbind(low_gf, df, high_gf)

    } else tmp_gf = df   # end if

    return(tmp_gf)
}