
#' @param out_path for path that stores out_YYYYMMDD
#' @param pproc_dir for path that store rds results 
if (F) {
    perturb_indx = 1
    perturb_type = 'emiss'
    timestr = all_timestr[1]
    overwriteTF = T
    transerrTF = F
    mx_res = 1
    ts_chem = NA
    eno_sf = 1
    pproc_dir = rds_path
}

postproc_chemv4 = function(site, timestr, met, store_path, pproc_dir, 
                           overwriteTF = T, 
                           proc_emis = list(c('edgar', 'epa', 'vulcan')), 
                           eno_sf = 1, ts_chem = NA, mx_res = 1, perturb_type, 
                           perturb_indx, transerrTF = FALSE, xstilt_wd) {
    
    setwd(xstilt_wd); source('r/dependencies.r') 
    #if (nchar(timestr) > 8) timestr = substr(timestr, 1, 8)
    print(timestr)
    
    perturb_emiTF = perturb_tsTF = perturb_mixTF = FALSE
    if ( !is.na(perturb_type)) {
        perturb_emiTF = ifelse(perturb_type == 'emiss', T, F)
        perturb_tsTF  = ifelse(perturb_type == 'chem', T, F)
        perturb_mixTF = ifelse(perturb_type == 'mixing', T, F)
    } 

    # get rds output files
    out_path = list.files(store_path, 'out', full.names = T)
    out_path = out_path[grepl('TROPOMI', out_path) & 
                        grepl(met, out_path) & 
                        grepl(timestr, out_path)]
    
    out_path = dirname(list.files(out_path, 'NO2', full.names = T))
    if (transerrTF) {
        out_path = out_path[grepl('outerr_', out_path)]
    } else out_path = out_path[grepl('out_', out_path)]
    met = strsplit.to.df(basename(out_path))$V3
    byid_path = file.path(out_path, 'NO2/by-id')

    # work on EDGAR-based results if there are
    if ( 'edgar' %in% unlist(proc_emis) ) {
        rds_nms = basename(list.files(byid_path, 'edgar_', recursive = T))
        if (length(rds_nms) > 0) {
            cat('working on EDGAR-based simulations...\n')
            rds_patt = naming_fnsv2(aq_invent = 'edgar', ghg_invent = 'edgar', 
                                    mx_res, eno_sf, ts_chem, perturb_indx, 
                                    perturb_emiTF, perturb_tsTF, perturb_mixTF)
            print(rds_patt)
            edgar_fn = obs_sim_pairv6(site, timestr, byid_path, pproc_dir, 
                                      rds_patt, met, overwriteTF, transerrTF)
        }
    }
    
    # work on EPA-based results if there are
    if ( 'epa' %in% unlist(proc_emis) ) {
        rds_nms = basename(list.files(byid_path, 'epa_', recursive = T))
        if (length(rds_nms) > 0) {
            cat('working on EPA-based simulations...\n')
            rds_patt = naming_fnsv2(aq_invent = 'epa', ghg_invent = 'epa', 
                                    mx_res, eno_sf, ts_chem, perturb_indx, 
                                    perturb_emiTF, perturb_tsTF, perturb_mixTF)
            print(rds_patt)
            epa_fn = obs_sim_pairv6(site, timestr, byid_path, pproc_dir, 
                                    rds_patt, met, overwriteTF, transerrTF)
        } 
    }


    # work on Vulcan-based results if there are
    if ( 'vulcan' %in% unlist(proc_emis) ) {
        rds_nms = basename(list.files(byid_path, 'vulcan_', recursive = T))
        if (length(rds_nms) > 0) {
            cat('working on Vulcan-based simulations...\n')
            rds_patt = naming_fnsv2(aq_invent = 'vulcan', ghg_invent = 'vulcan',
                                    mx_res, eno_sf, ts_chem, perturb_indx, 
                                    perturb_emiTF, perturb_tsTF, perturb_mixTF)
            print(rds_patt)
            vulcan_fn = obs_sim_pairv6(site, timestr, byid_path, pproc_dir, 
                                       rds_patt, met, overwriteTF, transerrTF)
        } 
    }

    cat('DONE with post processing STILT-NOx...\n')
}