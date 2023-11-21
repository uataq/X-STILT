
# only when simulations for all receptors have completed !!!
# if for perturbation, postprocess results from individual .rds
# delete original .rds files from the by-id folder to save space
# --------------------------------------------------

p_post_proc = function(site, timestr, byid_path, met, rds_patt, perturb_emiTF, 
                       perturb_tsTF, perturb_mixTF, transerrTF, overwriteTF,
                       store_path = NA) {
    
    if (is.na(store_path)) store_path = dirname(dirname(dirname(byid_path)))
    check_trj = list.files(byid_path, 'X_traj.rds', recursive = T)
    check_sim = list.files(byid_path, rds_patt, recursive = T)

    # postprocessing location --------
    pproc_dir = file.path(store_path, 'rds')
    if ( perturb_emiTF | perturb_mixTF | perturb_tsTF ) {
        if ( perturb_emiTF ) pstr = 'emiss'
        if ( perturb_mixTF ) pstr = 'mixing'
        if ( perturb_tsTF )  pstr = 'chem'
        pproc_dir = file.path(pproc_dir, paste0('perturb_', pstr))
    }
    if (!dir.exists(pproc_dir)) dir.create(pproc_dir, recursive = T)
    
    # --------------------------------------------------
    # store postprocessed results when all sims are done to save space ---
    if ( length(check_sim) == length(check_trj) ) {
        
        cat('p_post_proc(): post-processing all simulations...\n')
        pproc_fn = obs_sim_pairv5(site, timestr = substr(timestr, 1, 8), 
                                  byid_path, pproc_dir, rds_patt, met, 
                                  overwriteTF, transerrTF)
        
        # if for perturbation runs, delete original rds
        # if for best-estimation runs, sym link to simulation folder
        if ( length(pproc_fn) == 0 ) {
            if ( perturb_emiTF | perturb_tsTF | perturb_mixTF ) {
                file.remove(file.path(byid_path, check_sim))
            } else {
            
                # create symbolic links for the simulation 
                sim_path = file.path(dirname(byid_path), 'simulation')
                dir.create(sim_path, showWarnings = F)
                file.symlink(from = check_sim, to = sim_path)
            }   # end if delete or sym link
        }   # end if

        run_missingTF = FALSE
        cat('p_post_proc(): DONE with post-processing...\n')

    } else {
        
        run_missingTF = TRUE 
        
        # check to see if post-processing file exists
        pproc_patt = gsub('.rds', '', rds_patt)
        pproc_patt = gsub('_pchem_', '', pproc_patt)
        pproc_patt = paste0(met, '_', pproc_patt)
        if ( transerrTF ) pproc_patt = paste0(pproc_patt, '_transerr')
        pproc_fn = file.path(pproc_dir, paste0(site, '_', timestr, '_', 
                                               pproc_patt, '.rds'))
        if (file.exists(pproc_fn)) run_missingTF = FALSE
    }   # end if end of all simulations

    return(run_missingTF)
}
