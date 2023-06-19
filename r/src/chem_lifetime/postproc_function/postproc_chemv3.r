
#' @param out_path for path that stores out_YYYYMMDD
#' @param rds_path for path that store rds results 
if (F) {
    perturb_indx = 1
    perturb_type = 'emiss'
    timestr = all_timestr[1]
    overwriteTF = T
    site_lon = lon_lat$site_lon 
    site_lat = lon_lat$site_lat
    transerrTF = F
    mx_res = 1
}

postproc_chemv3 = function(site, timestr, met, store_path, rds_path, 
                           overwriteTF = T, 
                           proc_emis = list(c('edgar', 'epa')), 
                           mx_res = 1, 
                           perturb_indx = NA, 
                           perturb_type = c(NA, 'emiss', 'tau', 'mixing')[1],
                           transerrTF = FALSE, 
                           site_lon, site_lat, xstilt_wd) {

    setwd(xstilt_wd); source('r/dependencies.r')    # source all functions
    print(timestr)

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
    
    perturb_emiss_indx = perturb_tau_indx = perturb_mix_indx = NA
    if (!is.na(perturb_type)) {
        if (perturb_type == 'emiss') perturb_emiss_indx = perturb_indx 
        if (perturb_type == 'tau') perturb_tau_indx = perturb_indx 
        if (perturb_type == 'mixing') {
        perturb_mix_indx = perturb_indx; mx_res = NA }
    } 

    # work on EDGAR-based results if there are
    if ( 'edgar' %in% unlist(proc_emis) ) {
        rds_nms = basename(list.files(out_path, 'edgar_', recursive = T))
        if (length(rds_nms) > 0) {
            cat('working on EDGAR...\n')
            rds_patt = naming_fns(aq_invent = 'edgar', ghg_invent = 'edgar', 
                                  mx_res, eno_sf = 1, tau_chem = NA, 
                                  perturb_emiss_indx, perturb_tau_indx, 
                                  perturb_mix_indx)
            rds_patt = gsub('.rds', '', rds_patt)
            rds_patt = gsub('_pchem_edgar_', '', rds_patt); print(rds_patt)
        
            edgar_df = obs_sim_pairv5(site, timestr, out_path, rds_path, 
                                      rds_patt, met, overwriteTF, transerrTF, 
                                      site_lon, site_lat)
        }
    }
    
    # work on EPA-based results if there are
    if ( 'epa' %in% unlist(proc_emis) ) {
        rds_nms = basename(list.files(out_path, 'epa_', recursive = T))
        if (length(rds_nms) > 0) {
            cat('working on EPA...\n')
            rds_patt = naming_fns(aq_invent = 'epa', ghg_invent = 'epa',        
                                  mx_res, eno_sf = 1, tau_chem = NA, 
                                  perturb_emiss_indx, perturb_tau_indx, 
                                  perturb_mix_indx)
            rds_patt = gsub('.rds', '', rds_patt)
            rds_patt = gsub('_pchem_epa_', '', rds_patt); print(rds_patt)

            epa_df = obs_sim_pairv5(site, timestr, out_path, rds_path, 
                                    rds_patt, met, overwriteTF, transerrTF, 
                                    site_lon, site_lat)
        } 
    }

    cat('DONE with post processing STILT-NOx...\n')
}