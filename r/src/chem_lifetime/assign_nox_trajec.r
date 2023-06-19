# assign NOx to trajectories, DW, 01/13/2020
# only WACCM has been implemented

assign_nox_trajec <- function(p, r_run_time, ctm.fn, xstilt_wd) {

    cat(paste('assign_nox_trajec(): loading', ctm.fn), '\n')
    xmn = min(p$long) - 1.5; xmx = max(p$long) + 1.5 
    ymn = min(p$lati) - 1.5; ymx = max(p$lati) + 1.5
    all_run_times = r_run_time + p$time * 60
    tmn = min(all_run_times) - 6 * 3600
    tmx = max(all_run_times) + 6 * 3600

    # grab NOx fields from CTMs
    ctm.df = NULL
    for (c in ctm.fn) 
        ctm.df <- rbind(ctm.df, grab_waccm_nox(c, xmn, xmx, ymn, ymx, tmn, tmx))
    
    ctm.long <- unique(ctm.df$lon)
    ctm.lati <- unique(ctm.df$lat)
    cat(paste('*** debugging: for', as.numeric(p[p$indx == 1 & p$time == -1, 'lati']), 
              '; CTM.df', nrow(ctm.df), 'rows\n\n'))

    # find the cloest CTM hour to the particle time 
    # so we calculate the mid-point as the interval coordinate
    ctm.time  <- unique(ctm.df$time)
    ctm.mtime <- sort(zoo::rollmean(ctm.time, 2))

    # then grab [NOx] based on lon, lat, and pressure of particles with non-zero foot
    # add actual UTC time of each particle, need to convert 'time' from mins to sec
    # locate lon/lon grid and hours first and then look for vertical level
    p.prep <- p %>% mutate(run_time_par = r_run_time + time * 60, 
                           ctm.lati = ctm.lati[findInterval(lati, ctm.lati)], 
                           ctm.long = ctm.long[findInterval(long, ctm.long)], 

                           # use the middle time as a coordinate, need to plus 1
                           time.indx = findInterval(run_time_par, ctm.mtime) + 1, 
                           ctm.time = ctm.time[time.indx])

    cat('assign_nox_trajec(): Figuring out unique WACCM grids that particles fall into...\n')
    uni.grid <- p.prep[, c('ctm.long', 'ctm.lati', 'ctm.time')] %>% unique()
    cat(paste('*** debugging: for', as.numeric(p[p$indx == 1 & p$time == -1, 'lati']), 
              ';', nrow(uni.grid), '\n\n'))

    ### -----------------------------------------------------------------------
    # then loop over each unique WACCM grid and locate the pressure layer 
    # where particles fall into
    cat('assign_nox_trajec(): Looking for NOx at each particle at each timestamp...this takes a while\n')
    
    p.chem <- NULL
    for (r in 1 : nrow(uni.grid)) {

        if (r %% 10 == 0) cat(paste(signif(r / nrow(uni.grid) * 100, 3), '% done...\n'))
        tmp.ctm <- ctm.df %>% filter(lon == uni.grid$ctm.long[r], 
                                     lat == uni.grid$ctm.lati[r], 
                                     time == uni.grid$ctm.time[r]) 
        
        # 88 layers, 89 levels, always use lower press border as vertical coordinates
        lower.pres <- c(zoo::rollmean(tmp.ctm$pmid, 2), tmp.ctm$ps[1])
        tmp.ctm <- tmp.ctm %>% mutate(pres = lower.pres) %>% arrange(pres) %>% 
                    dplyr::select(lon, lat, pres, time, nox, ps)
        tmp.pres <- tmp.ctm$pres
        
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

    return(p.chem)
}