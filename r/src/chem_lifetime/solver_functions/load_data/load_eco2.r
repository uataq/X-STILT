
load_eco2 = function(timestr, invent = c('edgar', 'epa', 'odiac', 'vulcan'), 
                     eco2_edgar_fn, eco2_fn2 = NA, epa_fn = NA, epa_name = NA,
                     epa_tz = NA, xmin, xmax, ymin, ymax) {
    
    # loading EDGAR emission ------------------------------------------
    edgar_rt = load_eco2_edgar(eco2_edgar_fn)
    
    # if using CO2 emissions from ODIAC 
    if ( invent == 'odiac' ) {
        odiac_rt = load_eco2_odiac(eco2_fn2, timestr, xmin, xmax, ymin, ymax)
    } else odiac_rt = NULL  # end if
    
    # if using CO2 emissions from Vulcan 
    if ( invent == 'vulcan' ) {     # only for annual total emissions
        vulcan_rt = load_eco2_vulcan(eco2_fn2, xmin, xmax, ymin, ymax)
    } else vulcan_rt = NULL 
    
    # if scaled by EPA -------
    # ad hoc fix for Hunter Power plant, simply shift the EDGAR emission grid eastward by 0.1degree, DW, 03/27/2022
    if (epa_name == 'Hunter' & invent == 'epa') {
        edgar_rt = crop(edgar_rt, extent(xmin - 0.1, xmax + 0.1, ymin, ymax))  
        extent(edgar_rt)[1:2] = extent(edgar_rt)[1:2] + 0.1

    } else edgar_rt = crop(edgar_rt, extent(xmin, xmax, ymin, ymax))  

    if ( !is.na(epa_fn) & invent == 'epa' ) {
        epa_list = correct_emiss_epav2(edgar_rt, epa_name, epa_tz, epa_fn,'CO2')
        epa_df = epa_list$epa_df 
        edgar_pp = epa_list$emiss_pp 
    } else epa_df = edgar_pp = NULL
    

    # returning ------
    emiss_list = list(edgar_rt = edgar_rt, odiac_rt = odiac_rt, 
                      vulcan_rt = vulcan_rt, epa_df = epa_df, 
                      edgar_pp = edgar_pp)

    return(emiss_list)
}
