

load_chem_out = function(site, timestr, lon_lat, met, emiss, out_path, rds_path,
                         bg_quad, xbuf = 0.3, ybuf = 0.3, obs_tno2_df = NULL) {
    
    # GFS EDGAR, GFS EPA, HRRR EPA
    
    #pr_edgar = obs_sim_pairv4(site, timestr, out_path, rds_path, 'edgar_1mixing', met, F)$x_df

    if (emiss == 'epa')
        pr_epa = obs_sim_pairv4(site, timestr, out_path, rds_path, 
                                'epa_1mixing', met, F)$x_df
    #if (is.null(pr_edgar)) return()

    # define background using upwind soundings and fit linear regression -----
    pr_list = bg_def(pr_edgar, pr_epa, bg_quad, lon_lat, xbuf, ybuf,obs_tno2_df)
    
    # prep obs-sim -----------------------------------------------
    pr_df = prep_sim_obsv3(pr_list$pr_edgar, pr_list$pr_epa) %>% 
            dplyr::select(-ends_with('xco2'), -ends_with('ratio'))
    
    if (!is.null(pr_epa))
        colnames(pr_df)[grepl('epa_tno2', colnames(pr_df))] = 
            paste0(toupper(met), '_', colnames(pr_df)[grepl('epa_tno2', colnames(pr_df))])
    
    if (!is.null(pr_edgar))
        colnames(pr_df)[grepl('edgar_tno2', colnames(pr_df))] = 
            paste0(toupper(met), '_', colnames(pr_df)[grepl('edgar_tno2', colnames(pr_df))])

    return(list(pr_df = pr_df, xmn = pr_list$xmn, xmx = pr_list$xmx, 
                               ymn = pr_list$ymn, ymx = pr_list$ymx))
}   