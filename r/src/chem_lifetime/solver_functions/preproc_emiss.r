
# prepare and subset emissions over trajec domains to reduce computation cost
# CO2 and NO2 emissions are required, CO and CH4 are optional
preproc_emiss = function(site, timestr, traj_info, ghg_invent, aq_invent, 
                         eco2_edgar_fn, eco2_nonedgar_fn = NULL, 
                         eco_edgar_fn = NA, eno_edgar_fn, ech4_edgar_fn = NA, 
                         epa_fn, epa_name, epa_tz, out_path, overwriteTF = T) {
    
    # ----------------------------------------------------------------------
    tmp_eno_fn = file.path(out_path, paste0(aq_invent, '_NOx.tif'))
    tmp_eco_fn = file.path(out_path, paste0(aq_invent, '_CO.tif'))
    tmp_eco2_fn = file.path(out_path, paste0(aq_invent, '_CO2.tif'))
    tmp_ech4_fn = file.path(out_path, paste0(aq_invent, '_CH4.tif'))

    # ----------------------------------------------------------------------
    tmp_p = readRDS(traj_info$fn[1])$particle
    xmin = min(tmp_p$long) - 3
    xmax = max(tmp_p$long) + 3
    ymin = min(tmp_p$lati) - 3
    ymax = max(tmp_p$lati) + 3
    
    # load EDGAR CO2 emissions first --------------------------------------
    if (is.null(eco2_edgar_fn)) stop('missing EDGAR CO2 emissions...\n')
    eco2_lst = load_eco2(timestr, invent = ghg_invent, eco2_edgar_fn, 
                         eco2_fn2 = eco2_nonedgar_fn, epa_fn, epa_name, epa_tz, 
                         xmin, xmax, ymin, ymax)
    eco2_edgar = eco2_lst$edgar_rt     # EDGAR eCO2
    
    # load EDGAR NO emissions ----------------------------------------------
    if (is.null(eno_edgar_fn)) stop('missing EDGAR NO emissions...\n')
    eno_lst = load_eno(invent = aq_invent, emiss_fn = eno_edgar_fn, epa_fn, 
                       epa_name, epa_tz, xmin, xmax, ymin, ymax)
    eno_edgar = eno_lst$emiss_rt    # EDGAR eNO
    
    # load EDGAR CO and CH4 emissions ------------------------------------------
    if (!is.na(eco_edgar_fn)) {
        eco_edgar = load_eco(aq_invent, eco_edgar_fn, xmin, xmax, ymin, ymax)
    } else eco_edgar = NA
        
    if (!is.na(ech4_edgar_fn)) {
        ech4_edgar = load_ech4(aq_invent, ech4_edgar_fn, xmin, xmax, ymin, ymax)
    } else ch4_edgar = NA   
    
    # ----------------------------------------------------------------------
    # create CO and NOx emissions for non-EDGAR inventories
    if (aq_invent == 'edgar' | aq_invent == 'epa') 
        emis_stk = stack(eco2_edgar, eco_edgar, ech4_edgar, eno_edgar)

    if (aq_invent == 'odiac') 
        emis_stk = create_aq_odiac(eco2_non_edgar = eco2_lst$odiac_rt, 
                                   eco2_edgar, eco_edgar, ech4_edgar, eno_edgar)
    
    if (aq_invent == 'vulcan')
        emis_stk = create_aq_non_edgar(eco2_non_edgar = eco2_lst$vulcan_rt, 
                                       eco2_edgar, eco_edgar, eno_edgar)
    names(emis_stk) = c('eco2', 'eco', 'ech4', 'eno')

    # store the tmp file of emissions 
    if (overwriteTF) {
        writeRaster(emis_stk$eco2, tmp_eco2_fn, format = 'GTiff', overwrite = T)
        writeRaster(emis_stk$ech4, tmp_ech4_fn, format = 'GTiff', overwrite = T)
        writeRaster(emis_stk$eco, tmp_eco_fn, format = 'GTiff', overwrite = T)
        writeRaster(emis_stk$eno, tmp_eno_fn, format = 'GTiff', overwrite = T)
    }
    
    if (ghg_invent == 'epa' | aq_invent == 'epa') {
        epa_co2_df = eco2_lst$epa_df 
        epa_no_df = eno_lst$epa_df 
        eno_pp = eno_lst$emiss_pp    # EDGAR eNO for PP if needs to be scaled 
        epa_lst = list(epa_co2_df = epa_co2_df, epa_no_df = epa_no_df, 
                       eno_pp = eno_pp)
    } else epa_lst = NULL

    return(list(eco2_fn = tmp_eco2_fn, ech4_fn = tmp_ech4_fn, 
                eco_fn = tmp_eco_fn, eno_fn = tmp_eno_fn, epa_lst = epa_lst))
}
