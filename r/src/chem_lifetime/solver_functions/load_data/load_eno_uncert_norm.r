
if (F) {
    
    eno_edgar_fn = eno_fn
    xmin = lon_lat$minlon
    xmax = lon_lat$maxlon
    ymin = lon_lat$minlat
    ymax = lon_lat$maxlat

}

load_eno_uncert_norm = function(site, timestr, eno_edgar_fn, eco2_edgar_fn, 
                                eco2_odiac_fn, aq_invent, xmin, xmax, ymin, 
                                ymax, plotTF = F) {
    
    cat('load_eno_uncert(): use SD between EDGAR and ODIAC-scaled ENO to create grid-level fractional uncertainty...\n\n')
    edgar_eno = load_eno(invent = 'edgar', emiss_fn = eno_edgar_fn, epa_fn = NA,
                         epa_name = site, epa_tz = NA, xmin, xmax, ymin, ymax)$emiss_rt

    edgar_eco2 = load_eco2(timestr, invent = 'edgar', eco2_edgar_fn, 
                           eco2_odiac_fn, epa_fn = NA, epa_name = site, 
                           epa_tz = NA, xmin, xmax, ymin, ymax)$edgar_rt
    
    odiac_eco2 = load_eco2(timestr, invent = 'odiac', eco2_edgar_fn, 
                           eco2_odiac_fn, xmin = extent(edgar_eno)[1], 
                           xmax = extent(edgar_eno)[2], 
                           ymin = extent(edgar_eno)[3], 
                           ymax = extent(edgar_eno)[4])$odiac_rt

    # EDGAR based emission ratio, only for high ENO grids
    erno = edgar_eno / edgar_eco2         # EDGAR based spatially-exp ER
    odiac_eco2_agg = aggregate(odiac_eco2, fact = 12, fun = mean)
    odiac_eno = crop(odiac_eco2_agg * erno, extent(edgar_eno))     
    eno_stk = stack(edgar_eno, odiac_eno)
    names(eno_stk) = c('EDGAR_ENO', 'ODIACscaled_ENO')

    # simple SD between two NOx and CO2 emissions
    sd_eno  = calc(eno_stk, fun = sd)
    dff_eno = edgar_eno - odiac_eno

    eco2_stk = stack(edgar_eco2, odiac_eco2_agg)
    names(eco2_stk) = c('EDGAR_ECO2', 'ODIACagg_ECO2')
    sd_eco2 = calc(eco2_stk, fun = sd)
    dff_eco2 = edgar_eco2 - odiac_eco2_agg

    if (aq_invent == 'edgar') { cv_eno = sd_eno / edgar_eno ; cv_eco2 = sd_eco2 / edgar_eco2 } 
    if (aq_invent == 'odiac') { cv_eno = sd_eno / odiac_eno ; cv_eco2 = sd_eco2 / odiac_eco2_agg } 
    #cv_eno[cv_eno > 100] = 100         # previous version

    if (plotTF) {
        # for plotting
        library(RColorBrewer)
        ecol = rev(RColorBrewer::brewer.pal(9, 'Spectral'))
        ccol = RColorBrewer::brewer.pal(9, 'PuRd')
        dcol = rev(brewer.pal(11, 'PRGn'))

        nn = getValues(edgar_eno * 1e3); nn = nn[nn > 0]
        zlim = c(min(log10(nn)), max(log10(nn))); print(zlim)
        axag = list(at = seq(-3, 4), labels = 10 ^ seq(-3, 4))

        rr = getValues(odiac_eco2); rr = rr[rr > 0]
        zlim2 = c(min(log10(rr)), max(log10(rr))); print(zlim2)
        axag2 = list(at = seq(-2, 3), labels = 10 ^ seq(-2, 3))

        png(paste0('eno_prior_norm_', site, '_18', substr(timestr, 5, 6), 
                   '.png'), height = 2500, width = 1500, res = 300)
        par(mfrow = c(4, 2), family = 'Arial', mar = c(3, 2.5, 3, 4.))
        lgwd = 1.5

        e1 = plot(log10(edgar_eco2), zlim = zlim2, col = ecol, font.main = 1,
                  asp = 'fill', legend.width = lgwd, axis.args = axag2, 
                  main = 'a) 0.1deg EDGAR ECO2\n[umol m-2 s-1]')
        o1 = plot(log10(odiac_eco2), zlim = zlim2, font.main = 1, col = ecol,
                  asp = 'fill', legend.width = lgwd, axis.args = axag2, 
                  main = 'b) 1km ODIAC ECO2\n[umol m-2 s-1]')

        e2 = plot(log10(edgar_eno * 1e3), zlim = zlim, col = ecol, asp = 'fill',
                  font.main = 1, legend.width = lgwd, axis.args = axag,
                  main = 'c) 0.1deg EDGAR ENO\n[mmol m-2 s-1]')
        o2 = plot(log10(odiac_eno * 1e3), zlim = zlim, col = ecol, asp = 'fill',
                  font.main = 1, legend.width = lgwd, axis.args = axag,
                  main = 'd) 0.1deg ODIAC-based ENO\n[mmol m-2 s-1]')
        d1 = plot(dff_eco2, font.main = 1, col = dcol, legend.width = lgwd, 
                  zlim = c(-max(abs(getValues(dff_eco2))), 
                            max(abs(getValues(dff_eco2)))), asp = 'fill',
                  main = 'e) ECO2 (EDGAR - ODIAC)\n[umol m-2 s-1]')
        d2 = plot(dff_eno * 1e3, font.main = 1, col = dcol, legend.width = lgwd,
                  zlim = c(-max(abs(getValues(dff_eno))) * 1e3, 
                            max(abs(getValues(dff_eno))) * 1e3), asp = 'fill',
                  main = 'f) ENO (EDGAR - ODIAC)\n[mmol m-2 s-1]')
        # s2 = plot(sd_eco2, font.main = 1, col = ecol, legend.width = 2, 
        #           asp = 'fill', main = '0.1deg SD in ECO2\n[umol m-2 s-1]')
        # s1 = plot(sd_eno * 1e3, font.main = 1, col = ecol, asp = 'fill', 
        #           legend.width = 2, main = '0.1deg SD in ENO\n[mmol m-2 s-1]')
        
        
        # fraction uncert
        cv = cv_eno * 100
        cmax = trunc(quantile(getValues(cv), 0.9, na.rm = T))
        dc = 10; if (cmax > 150) dc = 50; if (cmax > 400) dc = 100
        brks = seq(0, cmax, dc)
        axc = list(at = brks, labels = c(brks[-length(brks)], 
                                         paste0('>', brks[length(brks)])))
        
        # only plot fraction uncertainty when ENO relatively large
        cv[cv > cmax] = cmax
        plt_stk = stack(cv, edgar_eno * 1e3); names(plt_stk) = c('cv', 'eno')
        plt_stk$cv[plt_stk$eno < 1] = NA 

        f1 = plot(plt_stk$cv, font.main = 1, col = ccol, asp = 'fill', 
                  axis.args = axc, zlim = c(0, cmax), legend.width = lgwd, 
                  main = paste('g) Fract uncert [%] for high ENO\nover', site))
        dev.off()
    }
    
    out_stk = stack(eno_stk, sd_eno, dff_eno, cv_eno, 
                    eco2_stk, sd_eco2, dff_eco2, cv_eco2)
    names(out_stk) = c('EDGAR_ENO', 'ODIACscaled_ENO', 'SD_ENO', 'DIFF_ENO', 'CV_ENO', 'EDGAR_ECO2', 'ODIACagg_ECO2', 'SD_ECO2', 'DIFF_ECO2', 'CV_ECO2')

    return(out_stk)
}