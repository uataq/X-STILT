
# find out receptors that have not been updated since last simulations according to a desired modified time @param mtime_string '%Y-%m-%d %M-%H-%S'
subset_recp = function(byid_path, rds_patt, traj_info, mtime_min = NULL) {
    
    # for running missing simulations - use this section 
    # if no simulation rds files found, simulate all receptors
    rds_fns = list.files(byid_path, rds_patt, full.names = T, recursive = T)
    if (length(rds_fns) == 0) return(traj_info)

    # get existing simulation info
    rds_info = strsplit.to.df(basename(rds_fns)) %>% 
               dplyr::select(time = V1, lon = V2, lat = V3) %>%
               mutate_if(is.character, as.numeric) %>% 
               mutate(mtime = file.info(rds_fns)$mtime)
    
    # check missing simulations
    df1 = rds_info[, c('lon', 'lat')]
    df2 = traj_info[, c('lon', 'lat')]
    colsToUse = intersect(colnames(df1), colnames(df2))
    common_indx = match(do.call("paste", df1[, colsToUse]), 
                        do.call("paste", df2[, colsToUse]))
    all_indx = 1:nrow(traj_info)
    ms_indx = all_indx[!all_indx %in% common_indx]

    if (length(ms_indx) > 0) {
        ms_info = traj_info[ms_indx, ]
    } else ms_info = NULL


    # check outdated simulations with incorrect modified time string 
    if (is.null(mtime_min)) mtime_min = max(unique(as.Date(rds_info$mtime)))
    outd_indx = which(as.Date(rds_info$mtime) < mtime_min)

    if (length(outd_indx) > 0) {
        outd_traj = list.files( dirname(rds_fns[outd_indx]), 'traj.rds', 
                                full.names = T, recursive = T )
        outd_info = rds_info[outd_indx, ] %>% mutate(fn = outd_traj)
        outd_info$mtime = NULL 
        ms_info = rbind(ms_info, outd_info)
    }

    return(ms_info)
}