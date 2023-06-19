
if (F) {

    eno_edgar_fn = eno_fn
    xmin = lon_lat$minlon
    xmax = lon_lat$maxlon
    ymin = lon_lat$minlat
    ymax = lon_lat$maxlat

}


# calculate emission uncertainties as a log-normal distribution rather than normal given large emission magnitude
load_eno_uncert_lognorm = function(site, timestr, eno_edgar_fn, eco2_edgar_fn, 
                             eco2_odiac_fn, aq_invent, xmin, xmax, ymin, ymax) {
    
    cat('load_eno_uncert(): use SD between EDGAR and ODIAC-scaled ENO to create grid-level fractional uncertainty...\n\n')
    edgar_eno = load_eno(timestr, invent = 'edgar', emiss_fn = eno_edgar_fn, 
                         epa_fn = NA, epa_name = site, epa_tz = NA, 
                         xmin, xmax, ymin, ymax)$emiss_rt
                         
    edgar_eco2 = load_eco2(timestr, invent = 'edgar', eco2_edgar_fn, 
                           eco2_odiac_fn, epa_fn = NA, epa_name = site, 
                           epa_tz = NA, xmin, xmax, ymin, ymax)$edgar_rt

    odiac_eco2 = load_eco2(timestr, invent = 'odiac', eco2_edgar_fn, 
                           eco2_odiac_fn, epa_fn = NA, epa_name = site, 
                           epa_tz = NA, xmin = extent(edgar_eno)[1], 
                           xmax = extent(edgar_eno)[2], 
                           ymin = extent(edgar_eno)[3], 
                           ymax = extent(edgar_eno)[4])$odiac_rt

    # EDGAR based emission ratio
    erno = edgar_eno / edgar_eco2         # EDGAR based ER
    odiac_eco2_agg = aggregate(odiac_eco2, fact = 12, fun = mean)
    odiac_eno = crop(odiac_eco2_agg * erno, extent(edgar_eno))    
    odiac_eno[odiac_eno == 0] = min(getValues(edgar_eno))
    
    # simple SD between two emissions
    eno_stk = stack(edgar_eno, odiac_eno, log10(edgar_eno), log10(odiac_eno))
    sd_eno  = calc(eno_stk[[1:2]], fun = sd)
    sdlog_eno  = calc(eno_stk[[3:4]], fun = sd)
    out_stk = stack(eno_stk, sd_eno, sdlog_eno)
    names(out_stk) = c('EDGAR_ENO', 'ODIACscaled_ENO', 'EDGAR_ENO_log10',
                       'ODIACscaled_ENO_log10', 'SD_ENO', 'SD_ENO_log10')

    if (F) {
        # for plotting
        library(RColorBrewer)
        ecol = rev(RColorBrewer::brewer.pal(9, 'Spectral'))
        zlim = c(-4.5, 0)
        cv2 = cv_eno * 100; cv2[cv2 > 500] = 500 
         dff_eno = edgar_eno - odiac_eno

        png('eno_prior.png', width = 2400, height = 1800, res = 300)
        par(mfrow = c(2, 3), family = 'Arial', mar = c(3, 3, 4, 5))
        o1 = plot(log10(odiac_eco2), font.main = 1, main = '1km ODIAC ECO2\n[umol m-2 s-1] in log10 space', col = ecol)
        e1 = plot(log10(edgar_eno), zlim = zlim, col = ecol, font.main = 1, 
                main = '0.1deg EDGAR ENO\n[umol m-2 s-1] in log10 space')
        o1p = plot(log10(odiac_eno), zlim = zlim, col = ecol, font.main = 1, 
                main = '0.1deg ODIAC-based ENO\n[umol m-2 s-1] in log10 space')
        d1 = plot(dff_eno, zlim = c(-0.04, 0.04), font.main = 1, 
                col = rev(brewer.pal(11, 'RdBu')), 
                main = 'diff (EDGAR - ODIAC)\nin ENO [umol m-2 s-1]', )
        f1 = plot(cv_eno * 100, font.main = 1, col = ecol, main = '0.1deg %uncertainty [%] in ENO\n= sd(EDGAR, ODIAC) / EDGAR')
        f2 = plot(cv2, font.main = 1, col = ecol, main = '0.1deg %uncertainty [%] in ENO\n= sd(EDGAR, ODIAC) / EDGAR')
        dev.off()
    }
    
    return(out_stk)
}