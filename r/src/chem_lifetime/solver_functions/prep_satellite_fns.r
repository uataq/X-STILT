
# always using TROPOMI v2
prep_satelite_fns = function(timestr, byid_path, tno2_path, tno2x_path, 
                             xch4_path, xco_path, xco2_path, lon_lat) {
    
    # check to see if trajectories are existed -------------------------------
    traj_fns = list.files(byid_path, 'X_traj.rds', full.names = T,recursive = T)
    if (length(traj_fns) == 0) { stop('NO trajec created yet...\n'); return() }

    traj_info = strsplit.to.df(basename(traj_fns)) %>%
                dplyr::select(time = V1, lon = V2, lat = V3) %>% 
                mutate_if(is.character, as.numeric) %>% mutate(fn = traj_fns)
    cat('prep_fns4no2(): Found', nrow(traj_info), 'receptors!\n')

    # check to see if satellite data are existed for simulations -------------
    no2_info = find.tropomi(tno2_path, timestr, lon_lat, bufferTF = F)
    co_info  = find.tropomi(xco_path, timestr, lon_lat, bufferTF = F)
    ch4_info = find.tropomi(xch4_path, timestr, lon_lat, bufferTF = F)
    co2_fn = list.files(xco2_path, paste0('_', substr(timestr, 3, 8), '_'), 
                        full.names = T, recursive = T)
    
    if (length(co2_fn) == 0) co2_fn = NA
    if (!is.null(co_info)) {
        co_fn = co_info[co_info$tot.count == max(co_info$tot.count), 'fn']
    } else co_fn = NA 
    if (!is.null(ch4_info)) {
        ch4_fn = ch4_info[ch4_info$tot.count == max(ch4_info$tot.count), 'fn']
    } else ch4_fn = NA 

    if (is.null(no2_info)) { stop('NO TROPOMI NO2 files found...NO2 AKs are required for simulations\n\n'); return() }

    no2_fn = no2_info[no2_info$tot.count == max(no2_info$tot.count), 'fn']
    if (is.na(no2_fn)) { stop('NO TROPOMI NO2 files found...\n\n'); return() }

    # check aux no2 files -------
    tmp_receptor = readRDS(traj_info$fn[1])$receptor
    tmp_edpt = strftime(min(tmp_receptor$run_time), 'UTC', format = '%Y%m%d%H')
    aux_fns = list.files(tno2x_path, 'AUX_CTMANA', full.names = T, recursive =T)
    aux_fn = do.call(c, lapply(tmp_edpt, function(x) 
                aux_fns[grepl(paste0('CTMANA_', substr(x, 1, 8)), aux_fns)] ))

    if (length(aux_fn) == 0) {
        stop(paste('No TROPOMI auxiliary NO2 file found on', tmp_edpt, 
                   '...required for NO2 boundary conditions\n'))
        return()
    }

    # if multiple versions of CO are found, choose the reprocessed one
    # Jan 10, 2023, DW
    if ( length(co_fn) > 1) {
        cat('prep_fns4no2(): found multiple CO files, choose reprocessed version with filename of "RPRO"...\n')
        co_fn = co_fn[grepl('RPRO', co_fn)]
    }

    prep_list = list(traj_info = traj_info, no2_fn = no2_fn, ch4_fn = ch4_fn,
                     co_fn = co_fn, co2_fn = co2_fn)
    return(prep_list)
}


