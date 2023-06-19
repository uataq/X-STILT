
if (F) {
    
    edgar_eno_fn = eno_fn
    xmin = lon_lat$minlon
    xmax = lon_lat$maxlon
    ymin = lon_lat$minlat
    ymax = lon_lat$maxlat

}

load_eno_uncert_norm = function(site, timestr, edgar_eno_fn, eco2_edgar_fn, 
                                eco2_odiac_fn, aq_invent, xmin, xmax, ymin, 
                                ymax, plotTF = F) {
    
    cat('load_eno_uncert(): use SD between EDGAR and ODIAC-scaled ENO to create grid-level fractional uncertainty...\n\n')
    edgar_eno = load_eno(timestr, invent = 'edgar', emiss_fn = edgar_eno_fn, 
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
    eno_stk = stack(edgar_eno, odiac_eno)
    names(eno_stk) = c('EDGAR_ENO', 'ODIACscaled_ENO')

    # simple SD between two emissions
    sd_eno  = calc(eno_stk, fun = sd)
    dff_eno = edgar_eno - odiac_eno

    if (aq_invent == 'edgar') cv_eno = sd_eno / edgar_eno 
    if (aq_invent == 'odiac') cv_eno = sd_eno / odiac_eno 
    #cv_eno[cv_eno > 100] = 100         # previous version

    if (plotTF) {
        # for plotting
        library(RColorBrewer)
        ecol = rev(RColorBrewer::brewer.pal(9, 'Spectral'))
        zlim = c(-4.5, 0)

        png(paste0('eno_prior_norm_', site, '.png'), width = 2400, 
            height = 1800, res = 300)
        par(mfrow = c(2, 3), family = 'Arial', mar = c(3, 3, 4, 5))
        o1 = plot(log10(odiac_eco2), font.main = 1, main = '1km ODIAC ECO2\n[umol m-2 s-1] in log10 space', col = ecol)
        e1 = plot(log10(edgar_eno), zlim = zlim, col = ecol, font.main = 1, 
                main = '0.1deg EDGAR ENO\n[umol m-2 s-1] in log10 space')
        o1p = plot(log10(odiac_eno), zlim = zlim, col = ecol, font.main = 1, 
                main = '0.1deg ODIAC-based ENO\n[umol m-2 s-1] in log10 space')
        d1 = plot(dff_eno, zlim = c(-0.5, 0.5), font.main = 1, 
                  col = rev(brewer.pal(11, 'RdBu')), 
                  main = 'diff (EDGAR - ODIAC)\nin ENO [umol m-2 s-1]', )
        s1 = plot(sd_eno, font.main = 1, col = ecol, main = '0.1deg SD in ENO')
        f1 = plot(cv_eno, font.main = 1, col = ecol, main = '0.1deg coefficient of variation [unitless] in ENO\n= sd(EDGAR, ODIAC) / EDGAR')
        dev.off()
    }
    
    out_stk = stack(eno_stk, sd_eno, dff_eno, cv_eno)
    names(out_stk) = c('EDGAR_ENO', 'ODIACscaled_ENO', 'SD_ENO', 'DIFF_ENO', 'CV_ENO')

    return(out_stk)
}