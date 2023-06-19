# subroutine to read OH fields from chemical transport model or reanalysis
# including CAMS, WACCM, and TCR, DW, 10/04/2020

### ----------------------------------------------------------------------------
grab_oh_ctm <- function(ctm.fn, ctm = c('tcr', 'cams', 'waccm')[1], 
                        xmn, xmx, ymn, ymx, tmn, tmx, xstilt_wd) {

    ### -----------------------------------------------------------------------
    cat(paste('grab_oh_ctm(): loading prescribed OH fields from', toupper(ctm), '...\n'))
    ctm.df <- NULL
    
    if (ctm == 'tcr') {     # grab TCR, OH field in mol cm-3
        for (c in ctm.fn) ctm.df <- rbind(ctm.df, grab_tcr(ctm.fn, xmn, xmx, ymn, ymx, tmn, tmx)) 
    
    } else if (ctm == 'cams') {     # grab CAMS

        for (c in ctm.fn) 
            ctm.df <- rbind(ctm.df, grab_cams(c, xmn, xmx, ymn, ymx, tmn, tmx, 
                                              cams.speci = 'oh', xstilt_wd))
        ctm.df$pres <- ctm.df$lower.pres # always use the lower pressure boundary 

    } else if (ctm == 'waccm') {     # grab WACCM
        for (c in ctm.fn) ctm.df <- rbind(ctm.df, grab_waccm(c, xmn, xmx, ymn, ymx, tmn, tmx))
    
    } else {
        cat('grab_oh_ctm(): NO functions to read your CTM files...Check ctm_name\n')
        return()
    }  # end if

    return(ctm.df)
}


### ----------------------------------------------------------------------------
# subroutine to find the corresponding OH to particle location and time

# instead of using the starting hour of each hour interval, find the closest 
# hour to the particle time
wgt_foot_oh <- function(p, r_run_time, ctm, ctm.fn, xstilt_wd) {
    
    cat(paste('wgt_foot_oh(): loading', ctm.fn), '\n')
    xmn = min(p$long) - 1.5; xmx = max(p$long) + 1.5 
    ymn = min(p$lati) - 1.5; ymx = max(p$lati) + 1.5
    all_run_times = r_run_time + p$time * 60
    tmn = min(all_run_times) - 6 * 3600
    tmx = max(all_run_times) + 6 * 3600

    # grab OH fields from CTMs
    ctm.df   <- grab_oh_ctm(ctm.fn, ctm, xmn, xmx, ymn, ymx, tmn, tmx, xstilt_wd)
    ctm.long <- unique(ctm.df$lon)
    ctm.lati <- unique(ctm.df$lat)
    cat(paste('*** debugging: for', as.numeric(p[p$indx == 1 & p$time == -1, 'lati']), 
              '; CTM.df', nrow(ctm.df), 'rows\n\n'))

    # find the cloest CTM hour to the particle time 
    # so we calculate the mid-point as the interval coordinate
    ctm.time  <- unique(ctm.df$time)
    ctm.mtime <- sort(zoo::rollmean(ctm.time, 2))

    # then grab [OH] based on lon, lat, and pressure of particles with non-zero foot
    # add actual UTC time of each particle, need to convert 'time' from mins to sec
    # locate lon/lon grid and hours first and then look for vertical level
    p.prep <- p %>% mutate(run_time_par = r_run_time + time * 60, 
                           ctm.lati = ctm.lati[findInterval(lati, ctm.lati)], 
                           ctm.long = ctm.long[findInterval(long, ctm.long)], 

                           # use the middle time as a coordinate, need to plus 1
                           time.indx = findInterval(run_time_par, ctm.mtime) + 1, 
                           ctm.time = ctm.time[time.indx])

    cat(paste('wgt_foot_oh(): Figuring out unique', toupper(ctm), 'grids that particles fall into...\n'))
    uni.grid <- p.prep[, c('ctm.long', 'ctm.lati', 'ctm.time')] %>% unique()
    cat(paste('*** debugging: for', as.numeric(p[p$indx == 1 & p$time == -1, 'lati']), 
              ';', nrow(uni.grid), '\n\n'))

    ### -----------------------------------------------------------------------
    # then loop over each unique WACCM grid and locate the pressure layer 
    # where particles fall into
    cat('wgt_foot_oh(): Looking for OH at each particle at each timestamp...this takes a while\n')
    
    p.chem <- NULL
    for (r in 1 : nrow(uni.grid)) {

        if (r %% 10 == 0) cat(paste(signif(r / nrow(uni.grid) * 100, 3), '% done...\n'))
        tmp.ctm <- ctm.df %>% filter(lon == uni.grid$ctm.long[r], 
                                     lat == uni.grid$ctm.lati[r], 
                                     time == uni.grid$ctm.time[r]) 
        
        if (ctm == 'waccm') {   

            # 88 layers, 89 levels, always use lower press border as vertical coordinates
            lower.pres <- c(zoo::rollmean(tmp.ctm$pmid, 2), tmp.ctm$ps[1])
            tmp.ctm <- tmp.ctm %>% mutate(pres = lower.pres) %>% arrange(pres) %>% 
                       dplyr::select(lon, lat, pres, time, oh, ps)
            tmp.pres <- tmp.ctm$pres

        } else {

            # SOMETIMES, TCR gives NA OH for the layers near the surface
            # overwrite to use the OH for the lowest available model layer
            tmp.ctm <- tmp.ctm %>% arrange(pres) %>% filter(!is.na(oh)) %>% 
                                    dplyr::select(lon, lat, pres, time, oh)
            tmp.pres <- tmp.ctm$pres    # pressure for lower pressure boundary
            #if (r %% 10 == 0) cat(paste(signif(r / nrow(uni.grid) * 100, 3), '% done2..\n'))
        }   # end if
        
        ## use middle pressure to locate which vertical box particle falls into
        p.tmp <- p.prep %>% filter(ctm.long == uni.grid$ctm.long[r], 
                                   ctm.lati == uni.grid$ctm.lati[r], 
                                   ctm.time == uni.grid$ctm.time[r])
        
        ## locate level indx
        p.tmp <- p.tmp %>% 
                 mutate(pres.indx = findInterval(pres, tmp.pres) + 1, 
                        
                        # if pres.indx is NA or larger than # of levs, 
                        # meaning particle is below CTM-based Psfc, 
                        # so we simply replace NA pressure with sfc pressure 
                        pres.indx = ifelse(pres.indx > length(tmp.pres), 
                                           length(tmp.pres), pres.indx), 
                        
                        ctm.pres = tmp.pres[pres.indx]) %>% 
                 left_join(tmp.ctm, by = c('ctm.lati' = 'lat', 'ctm.long' = 'lon', 
                                           'ctm.time' = 'time', 'ctm.pres' = 'pres')) 

        p.chem <- rbind(p.chem, p.tmp)
        #if (r %% 10 == 0) cat(paste(signif(r / nrow(uni.grid) * 100, 3), '% done3..\n'))
    }     # end for r

    # sanity check, diff < 0, abs(diff) should be smaller than horizonal resolution
    #range(p.chem$ctm.lati - p.chem$lati)
    #range(p.chem$ctm.long - p.chem$long)
    #range(p.chem$ctm.pres - p.chem$pres)

    ### ------------------------------------------------------------------------
    # convert OH mass mixing ratio to mol mol-1 and 
    # then to number density in molec cm-3 using Av * p / (R * T)
    # need to multiple OH mixing ratio with the number density of dry air
    M_air = 29              # molar mass of air in 28.97 g mol-1
    Av <- 6.02214 * 1E23    # molec per mol
    R  <- 8.31              # universal gas const in J / mol / K -- Av * Boltzmann const

    # convert WACCM OH from mol mol-1 to molec cm-3; CAMS OH from kg kg-1 to molec cm-3
    cat('wgt_foot_oh(): calculating instant rate coefficient and scaling footprint...\n')
    if (ctm == 'waccm') p.chem$oh <- p.chem$oh * Av * p.chem$pres * 100 / R / p.chem$temz / 1E6
    if (ctm == 'cams') p.chem$oh <- p.chem$oh / 17 * M_air * Av * p.chem$pres * 100 / R / p.chem$temz / 1E6

    # now calculate the rate constant K(T) using 2.8 * 1E-11 * (T_K / 300) ^ (-1.3)      
    # K should have a unit of cm3 / molec / s
    p.out <- p.chem %>% 
             mutate(k = 2.8 * 1E-11 * (temz / 300) ^ (-1.3), 

                    # because both K and n_OH vary with time backward, 
                    # we need to perform integral by calculating the
                    # product of k and [OH] at each timestamp   
                    kn = k * oh ) %>% 
             dplyr::select(-c('ctm.lati', 'ctm.long', 'ctm.time', 'ctm.pres', 'pres.indx'))

    ### ------------------------------------------------------------------------
    # since Tau relies on time backward, we perform integral by 
    # calculating the accumulative sum of tau and dt from 
    # receptor time to the time backward
    delt <- max(abs(unique(diff(p$time)))) * 60  # 60 sec

    # for cumsum(), we must order the time backward, DW, 10/05/2020
    # convert delt from mins to second
    # calculate the integral of instant K * instant [OH] * delt, which is unitless
    p.cumsum <- p.out %>% mutate(abs_time = abs(time)) %>% arrange(abs_time) %>% 
                group_by(indx) %>% mutate(tot_knt = cumsum(kn * delt)) %>% ungroup()            

    # finally weight footprint, this takes a while
    p.chem.wgt <- p.cumsum %>% mutate(sf_lifetime = exp(-tot_knt), 
                                      foot_chem = foot * sf_lifetime) %>%
                                      # convert UTC time to local time 
                                      #local_time = format(run_time_par, tz = tz)) %>% 
                  rename(foot_before_chem = foot, foot = foot_chem) %>% 

                  # after the renaming, foot_before_wgt is the foot before AKPWF weighting 
                  # foot_before_chem is the foot before NOx lifetime after AKPWF weighting 
                  # foot is the weighting foot after considering the NOx lifetime
                  # DW, 10/06/2020
                  dplyr::select(time, indx, long, lati, xhgt, pres, zagl, zsfc, 
                                temz, mlht, k, oh, tot_knt, sf_lifetime, wgt_no2, 
                                foot_before_wgt, foot_before_chem, foot)

    return(p.chem.wgt)

}

