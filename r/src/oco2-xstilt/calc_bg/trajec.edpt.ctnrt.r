#### subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
# use CO2_total mole fractions 3D fields to get the CO2 concentration
# for STILT particle endpoints --> trajec-endpoint boundary condition
# ALSO, NEED TO WEIGHT THE ENDPOINT CONCENTRATION USING AK & PWF !!!

# for current purpose (Riyadh, Middle East), only use global 3x2 files,
# may want to add a nested NAM 1x1 for US sites

# Written by Dien Wu, 09/12/2016

########################################################################
# updates --
# Before, we only consider CO2 boundary condition for model levels
# using trajec endpoints, so AK aloft is ZERO
# HOWEVER, what we should do is to use CT for upper atmos
#  + keep original OCO-2 AK aloft
# last fix on 04/06/2017, DW
#
# separate STILT level and CT-NRT levels, and output two XCO2, 04/24/2017, DW
# clean codes, DW, 06/12/2018
# optimize codes, DW, 09/15/2018
# store results back to weighted trajec, DW, 09/16/2018
#########################################################################

# input variables needed--
#adjusted ak*pw profiles for STILT levels, .RData traj file, receptor time,

### list variable names for CO2_components,
# bg (background, component due to initial condition),
# bio (terrestrial biopshere exchange), ff (fossil fuel burning),
# fires (direct fire emissions), ocean (air-sea exchange),
# latitude, longitude, level, time (since 2000-01-01 00:00:00 UTC)

### list variable names for CO2_total,
# air_mass, blh (PBL thickness), co2, decimal_date,
# gph (geoponential_height), orography (sfc geopotential), pressure,
# specific_humidity, temperature, time_components

trajec.edpt.ctnrt <- function(traj.file, foot.ext, ct.ver, mf.path, timestr, 
  txtfile = NULL, writeTF = T){

    # from traj.file, get receptor info
    receptor <- unlist(strsplit(gsub('_X_wgttraj.rds', '', basename(traj.file)), '_'))
    receptor <- as.data.frame(matrix(receptor, byrow = T, ncol = 3),
        stringsAsFactors = F) %>% mutate_all(funs(as.numeric), colnames(receptor))
        # mutate_all() convert character to numberic
    colnames(receptor) <- list('timestr', 'lon', 'lat')

    order.index <- order(receptor$lat)
    receptor    <- receptor[order.index, ]
    traj.file   <- traj.file[order.index]
    receptor$xco2.bound <- NA
    receptor$xco2.ap <- NA 

    for (t in 1:length(traj.file)) {

        cat('Reading trajec..#', t, '..\n')
        output <- readRDS(traj.file[t])
        particle <- output$particle
        combine.prof <- output$wgt.prof

        storeTF <- T; if ('ctnrt.bound' %in% colnames(combine.prof)) storeTF <- F

        recp.yr  <- substr(timestr, 1, 4); recp.mon <- substr(timestr, 5, 6)
        recp.day <- substr(timestr, 7, 8); recp.hr  <- substr(timestr, 9, 10)

        # subset particle by footprint domain, to yield fair comparisons
        p <- particle %>% filter(long >= foot.ext[1], long <= foot.ext[2], 
                                 lati >= foot.ext[3], lati <= foot.ext[4])

        # --------- PART 1: DEALING WITH ENDPOINTS OF TRAJ for model levels -------- #
        # for grabbing the CO2 concentration at endpoints,
        # we need the lon, lat, asl height and absolute time of STILT particles
   
        # find out the min time according to index number
        edpt.time <- p %>% group_by(indx) %>% 
          dplyr::summarize(min.time = min(time)) %>% ungroup()

        # NOW, select endtraj rows
        edpt.p <- p %>% left_join(edpt.time, by = 'indx') %>% 
          group_by(indx) %>% filter(time == min.time) %>% ungroup() %>% arrange(indx)
  
        # THUS, endtraj data frame should always have numpar rows
        # compute the ASL heights for STILT particles
        edpt.p <- edpt.p %>% mutate(zasl = zagl + zsfc)
    
        # calculate backtime
        # add recp.time to edpt.time for local CT at upper levels
        recp.time <- as.POSIXct(timestr, format = '%Y%m%d%H', tz = 'UTC')
        edpt.time <- recp.time + c(0, edpt.p$time) * 60  
        edpt.timestr <- format(edpt.time, format = '%Y-%m-%d %H')
        edpt.times   <- data.frame(
            yr  = as.numeric(substr(edpt.timestr, 1, 4)), 
            mon = as.numeric(substr(edpt.timestr, 6, 7)),
            day = as.numeric(substr(edpt.timestr, 9, 10)),
            hr  = as.numeric(substr(edpt.timestr, 12, 13))) %>% mutate(ct.hr = hr)

        # fix actual hours to CT hours, if back.hours %% !=0, 08/29/2017, DW
        if (as.numeric(substr(timestr, 9, 10)) %% 3 != 0) {
            std.hr <- seq(0, 21, 3)  # standard hour from CT, 00, 03,..., 21 UTC
            hr.index <- edpt.times$ct.hr %% 3 != 0
            find.index <- findInterval(edpt.times[hr.index, 'ct.hr'], std.hr)
            edpt.times[hr.index, 'ct.hr'] <- std.hr[find.index]
        }
        
        # get unique day, year, mon for endpoints, and order them
        uni.edpt.day  <- sort(unique(edpt.times$day))	 # e.g., day 26, 27, 28
        uni.edpt.yr   <- sort(unique(edpt.times$yr))
        uni.edpt.mon  <- sort(unique(edpt.times$mon))
        uni.edpt.timestr <- sort(unique(substr(edpt.timestr, 1, 10)))

        # since CT files contain 3D CO2 fields every 3 hours in UTC,
        # assign a endpoint time index to each particle
        # Open all CT files into a list before assigning background CO2
        all.ct.gph <- rep(list(NULL), length(uni.edpt.day))
        all.ct.co2 <- all.ct.gph
        all.ct.pres <- all.ct.gph

        for (f in 1:length(uni.edpt.day)) {
            # find and open the CT file
            ctfile <- list.files(path = mf.path, pattern = uni.edpt.timestr[f], 
              full.names = T)
            ctdat <- nc_open(ctfile)

            # ---------------- readin CT-NRT for getting molefractions --------------- #
            # 2 deg in N-S, 3 deg in E-W, move centered lat, lon to lower left
            ct.lat  <- ncvar_get(ctdat, 'latitude') - 1
            ct.lon  <- ncvar_get(ctdat, 'longitude') - 1.5

            # CT time: centered time on 1:30, 4:30, 7:30, 10:30, 13:30, 16:30, 19:30 & 22:30
            # time range from 00-03, 03-06, 06-09, 09-12, 12-15, 15-18, 18-21 and 21-24.
            ct.time <- ncvar_get(ctdat, 'time_components')	# UTC time components
            rownames(ct.time) <- c('year', 'mon', 'day', 'hour', 'min', 'sec')
            ct.time <- data.frame(t(ct.time))
            ct.hr   <- ct.time$hour - 1 
            ct.min  <- ct.time$min - 30
            ct.level <- ncvar_get(ctdat, 'level')	       # 25 levels
            ct.bound <- ncvar_get(ctdat, 'boundary')	   # 26 boundaries

            # grab geopotentail height, 26 BOUNDS FOR HGT
            ct.gph <- ncvar_get(ctdat, 'gph')	           # [LON, LAT, BOUND, TIME]
            dimnames(ct.gph) <- list(ct.lon, ct.lat, ct.bound, ct.hr)

            # grab CO2 fields, 25 LEVELS FOR CO2
            ct.co2 <- ncvar_get(ctdat, 'co2')	           # [LON, LAT, LEVEL, TIME]
            dimnames(ct.co2) <- list(ct.lon, ct.lat, ct.level, ct.hr)

            # grab pressure levels, 
            ct.pres <- ncvar_get(ctdat, 'pressure') / 100 # convert to mb
            dimnames(ct.pres) <- list(ct.lon, ct.lat, ct.bound, ct.hr)

            # combine all arrays into list
            all.ct.gph[[f]] <- ct.gph
            all.ct.co2[[f]] <- ct.co2
            all.ct.pres[[f]] <- ct.pres

            nc_close(ctdat)
        } # end f for reading in ct files

        names(all.ct.gph) <- uni.edpt.day
        names(all.ct.co2) <- uni.edpt.day
        names(all.ct.pres) <- uni.edpt.day

        # convert to data.frame 
        ct.gph.df <- melt(all.ct.gph)
        ct.co2.df <- melt(all.ct.co2)
        ct.pres.df <- melt(all.ct.pres)
        colnames(ct.gph.df) <- list('lon', 'lat', 'level', 'hr', 'gph', 'date')
        colnames(ct.co2.df) <- list('lon', 'lat', 'level', 'hr', 'co2', 'date')
        colnames(ct.pres.df) <- list('lon', 'lat', 'level', 'hr', 'pres', 'date')
        sel.gph <- ct.gph.df %>% filter(lon >= foot.ext[1] - 3, lon <= foot.ext[2] + 3, 
                                        lat >= foot.ext[3] - 2, lat <= foot.ext[4] + 2)
        sel.co2 <- ct.co2.df %>% filter(lon >= foot.ext[1] - 3, lon <= foot.ext[2] + 3, 
                                        lat >= foot.ext[3] - 2, lat <= foot.ext[4] + 2)
        sel.pres <- ct.pres.df %>% filter(lon >= foot.ext[1] - 3, lon <= foot.ext[2] + 3, 
                                          lat >= foot.ext[3] - 2, lat <= foot.ext[4] + 2)

        ## use lat/lon/asl to locate the CT grid (all CT files have same lat lon grid)
        edpt.p <- cbind(edpt.p, edpt.times[-1, ])
        edpt.p <- edpt.p %>% mutate(
            ct.lati = as.numeric(ct.lat[findInterval(lati, ct.lat)]), 
            ct.long = as.numeric(ct.lon[findInterval(long, ct.lon)]), 
            ct.date = as.character(day), ct.co2 = NA) 
      
        # ---------- DEALING WITH CT 26 BOUNDARY AND 25 LEVEL for hgt.index -------- #
        # assign 'hgt.index' to each particles, referring which CT levels to grab,
        # by comparing ASL from STILT and GPH from CT
        # our goal is to relate CO2 concentration with the hgts, not levels
        cat('trajec.edpt.ctnrt(): Grabbing CT CO2 for each particle...it takes time...\n')

        # loop over each endpoint
        for (pp in 1:nrow(edpt.p)) {

            tmp.p <- edpt.p %>% filter(indx == pp) 
            tmp.gph <- sel.gph %>% filter(lon == tmp.p$ct.long, lat == tmp.p$ct.lati, 
                date == tmp.p$ct.date, hr == tmp.p$ct.hr)

            #ct.gph <- tmp.gph$gph[findInterval(tmp.p$zasl, tmp.gph$gph)]
            #if (length(ct.gph) == 0) ct.gph <- 1
            ct.level <- findInterval(tmp.p$zasl, tmp.gph$gph)
            if (ct.level == 0) ct.level <- 1
            
            # assign CO2 
            tmp.co2 <- sel.co2 %>% filter(lon == tmp.p$ct.long, lat == tmp.p$ct.lati, 
                       date == tmp.p$ct.date, hr == tmp.p$ct.hr, level == ct.level)
            edpt.p$ct.co2[pp] <- tmp.co2$co2
        }  # end for pp

        # take the mean of ensemble at each STILT releasing levels
        co2.bound.lower <- tapply(edpt.p$ct.co2, edpt.p$xhgt, mean)

        # because ct.bg.co2 is not normal distributed,
        # try 50th percentile, instead of mean, DW, 11/16/2017
        #co2.bound.lower <- tapply(edpt.p$ct.co2, edpt.p$xhgt, median)

        # --------- PART 2: DEALING WITH ENDPOINTS OF TRAJ for OCO-2 levels -------- #
        # grab CT pressure and interpolate CT at local receptors for high up
        cat('trajec.edpt.ctnrt(): Grabbing CT CO2 above max model level...\n')

        # NOW, grab the CT CO2 profiles at the recptor
        tmp.pres <- sel.pres %>% 
            filter(lon == ct.lon[findInterval(receptor$lon[t], ct.lon)], 
                   lat == ct.lat[findInterval(receptor$lat[t], ct.lat)], 
                   date == recp.day, 
                   hr == ct.hr[findInterval(recp.hr, ct.hr)])
        tmp.co2 <- sel.co2 %>% 
            filter(lon == ct.lon[findInterval(receptor$lon[t], ct.lon)], 
                   lat == ct.lat[findInterval(receptor$lat[t], ct.lat)], 
                   date == recp.day, 
                   hr == ct.hr[findInterval(recp.hr, ct.hr)])

        # interpolating CT CO2 profiles onto combine levels, only above STILT levels
        model.pres <- combine.prof$pres
        max.model.pres <- model.pres[max(which(combine.prof$stiltTF))]
        upper.pres <- model.pres[model.pres < max.model.pres]

        # assigning CT CO2 values to each model pressure, go through each
        co2.bound.upper  <- NULL
        for (l in 1:length(upper.pres)) {

            # locate which layer that model level falls in
            new.pres <- tmp.pres %>% mutate(diff.pres = pres - upper.pres[l]) %>% 
              filter(diff.pres > 0)  

            # --> lower boundary index, if no press, use sfc level
            if (nrow(new.pres) == 0) {
                level.index <- 1
            } else {
                level.index <- max(new.pres$level)
            }  # end if

            co2.bound.upper <- c(co2.bound.upper , 
                tmp.co2[tmp.co2$level == level.index, 'co2'])
        }   # end for l, assigning CO2 background for upper levels


        # ------------ PART 3: COMBINING UPPER AND LOWER bg CO2 ------------------- #
        combine.prof$ctnrt.bound <- c(co2.bound.lower, co2.bound.upper)

        # apply ak pw weighting for each level mean and add prior portion
        xco2.bound <- sum(combine.prof$ak.pwf * combine.prof$ctnrt.bound)
        xco2.ap <- sum((1 - combine.prof$ak.norm) * combine.prof$pwf * combine.prof$ap)
        
        receptor$xco2.bound[t] <- xco2.bound
        receptor$xco2.ap[t] <- xco2.ap

        # store results back to weighted trajec, DW, 09/16/2018
        output$wgt.prof <- combine.prof 
        if (storeTF) saveRDS(object = output, file = traj.file[t])

        gc()
    }  # end for t
 
    # finally, write in a txt file
    if (writeTF) 
      write.table(x = receptor, file = txtfile, sep = ',', row.names = F, quote = F)

    # write results in txtfile
    return(receptor)
}

# end of subroutine
