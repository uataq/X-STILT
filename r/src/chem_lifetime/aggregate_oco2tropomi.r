
agg_oco2tropomi = function(oco_df, tropomi_df, bg_xco2 = NULL) {

    if (!is.null(bg_xco2)) oco_df$xco2_ff = oco_df$xco2 - bg_xco2 

    # loop over each TROPOMI polygon and find the corresponding TROPOMI polygon 
    tropomi_add = NULL
    for (p in unique(tropomi_df$polygon)) {

        tmp = tropomi_df %>% filter(polygon == p)
        pipTF = point.in.polygon(oco_df$lon, oco_df$lat, tmp$lons, tmp$lats)
        pindx = which(pipTF > 0)
        
        if (length(pindx) == 0) next 
        tmp_oco = oco_df[pindx, ] %>% dplyr::select(contains('xco2'))
        #tmp_oco = oco_df[pindx, ] %>% dplyr::select(-c('lon', 'lat'))
        tmp_avg = as.data.frame.list(colMeans(tmp_oco, na.rm = T))
        tmp_add = cbind(tmp, tmp_avg)

        if (!is.null(bg_xco2)) 
            tmp_add$xco2_ff = as.numeric(tmp_avg[names(tmp_avg) == 'xco2_ff'])

        tropomi_add = rbind(tropomi_add, tmp_add)
    }

    return(tropomi_add)
}


