# function to obtain TROPOMI NO2 prior profile, DW, 01/11/2022
# requires particle-level humidity data!!! 
# requires auxiliary NO2 data from TROPOMI, not the L2 
#' if @param aux_fn is not specified, please provide @param aux_path 

if (F) {
    receptor = p_edpt
    aux_fn = aux_no2_fn
}

get.tropomi.prior.no2 = function(receptor, aux_path = NULL, aux_fn = NULL) {
    

    if (is.null(aux_fn)) {
        timestr = strftime(min(receptor$run_time), 'UTC', format = '%Y%m%d%H')
        fns = list.files(aux_path, 'AUX_CTMANA', full.names = T)

        aux_fn = do.call(c, lapply(timestr, function(x) 
            fns[grepl(paste0('CTMANA_', substr(x, 1, 8)), fns)] ))

        if (length(aux_fn) == 0) 
            stop('No auxiliary TROPOMI file found for this timestr...\n')
    } # end if


    # ------------- load all TROPOMI aux grid center lat/lon ----------------
    dat = nc_open(aux_fn[1])

    # move centered lat/lon to lower left lat/lon
    lat = ncvar_get(dat, 'lat') - 0.5   
    lon = ncvar_get(dat, 'lon') - 0.5   
    lev = ncvar_get(dat, 'lev')     
    day = ncvar_get(dat, 'time')   # days since 1950-01-01 00:00:00
    date = as.POSIXct(day * 24 * 3600, 'UTC', origin = '1950-01-01 00:00:00')

    # coefficient at layer interfaces and midpoints
    # hyam hybm (mlev=hyam+hybm*ps)
    hyai = ncvar_get(dat, 'hyai')       # Pa, vec of 34
    hybi = ncvar_get(dat, 'hybi')       # unitless
    hyam = ncvar_get(dat, 'hyam')       # Pa
    hybm = ncvar_get(dat, 'hybm')       # unitless
    ps = ncvar_get(dat, 'ps')           # surface pres in Pa, [lon, lat, day]
    
    # locate lat, long, pressure indics before loading NO2 
    receptor$lat_indx = findInterval(receptor$lati, lat)
    receptor$lon_indx = findInterval(receptor$long, lon)
    receptor$day_indx = findInterval(receptor$run_time, date)
    uni_indx = receptor %>% dplyr::select(lon_indx, lat_indx, day_indx) %>% 
               unique() %>% arrange(lon_indx, lat_indx)

    # [lon, lat, lev, day], 1x1, 34 hybrid level, layer midpoints, 30 mins
    # NO2 VMR in humid air in ppb
    receptor = do.call(rbind, lapply(1 : nrow(uni_indx), function(x){
        
        if (x %% 20 == 0) 
        cat(paste('# ---', signif(x/nrow(uni_indx), 3) * 1e2, '% DONE --- #\n'))

        tmp_indx = uni_indx[x, ]
        tmp_recp = receptor %>% filter(receptor$lon_indx == tmp_indx$lon_indx, 
                                       receptor$lat_indx == tmp_indx$lat_indx, 
                                       receptor$day_indx == tmp_indx$day_indx)

        # surface pressure and midpoint of pressure layer
        tmp_psfc = ps[tmp_indx$lon_indx, tmp_indx$lat_indx, tmp_indx$day_indx] 
        tmp_pmid = (hyam + hybm * tmp_psfc) / 100
        tmp_psfc = tmp_psfc / 100 

        # lower pressure boundaries 
        tmp_plwr = c(tmp_psfc, rep(NA, length(tmp_pmid)))
        for (y in 1 : length(tmp_pmid)) {
            tmp_plwr[y + 1] = tmp_pmid[y] * 2 - tmp_plwr[y]
            if (tmp_plwr[y + 1] < 0) tmp_plwr[y + 1] = 0
        }

        # sanity check for the lowest 4 layers, 5 boundaries
        # n = 5; plot(rep(1, n), tmp_plwr[1:n]); points(rep(1, n - 1), tmp_pmid[1: n - 1], col = 'red')
        
        # select NO2 profile for the nearest sounding to trajec endpoint
        tmp_no2 = ncvar_get(dat, 'no2')[tmp_indx$lon_indx, tmp_indx$lat_indx, , tmp_indx$day_indx]

        tmp_df = data.frame(no2_vmr_moist = tmp_no2, pres_mid = tmp_pmid, 
                            pres_lower = tmp_plwr[-length(tmp_plwr)]) %>% 
                            arrange(pres_lower)
        
        # locate a particular parcel onto TM5 pressure grid
        # use lower pressure boundary as a reference
        pres_indx = findInterval(tmp_recp$pres, tmp_df$pres_lower) + 1
        pres_indx[pres_indx > nrow(tmp_df)] = nrow(tmp_df)
        tmp_recp$no2_moist = tmp_df[pres_indx, 'no2_vmr_moist'] 
        return(tmp_recp)
    }))


    # --------- need to convert moist MR to dry MR -----------
    if ( 'sphu' %in% colnames(receptor) ) {     
        
        # use specific humidity to calc water vapor mixing ratio
        receptor = receptor %>% mutate(w = sphu / (1 - sphu))   

    } else if (length(unique(receptor$sphu)) == 1 | 
               'rhfr' %in% colnames(receptor)) { 
        
        # if specific humidity is zero, meaning no SH is modeled in met fields, 
        # use RH and saturation vapor pressure to calculate specific humidity 
        # Clausius–Clapeyron equation for saturation vapor pressure
        receptor = receptor %>% mutate(
                    temc = temz - 273.15, 
                    es = 6.112 * exp( (17.67 * temc) / (temc + 243.5) ), 
                    e = rhfr * es,               # vapor pressure in air
                    w = 0.622 * e / (pres - e)   # H2Ov mixing ratio
        )   

    } else stop('get.tropomi.no2.prior(): require trajec-level H2Ov info\n either from specific humidity ("sphu") or relative humidity ("rhfr"). Consider regenerating trajec...\n')  

    # moist MR * moist air / dry air, now in ppb
    # since MR of moist_air > dry_air, NO2 MR in moist air < NO2 in dry air
    # still a very small adjustmnt
    receptor$no2_dry = receptor$no2_moist * 1 / (1 - receptor$w)
    # plot( (receptor$no2_moist - receptor$no2_dry) * 1e9, receptor$pres, ylim = c(1013, 0), xlab = 'ppb')  # * 1e9 for ppb

    return(receptor)
} # end of subroutine
