# script to interpolate the winds from GDAS, using trajwind() for each receptor based on asosiosonde data
# read output model-data comparisons and make plots
# written by DW, 05/15/2017

# Add temporal interpolation of GDAS winds, DW, 05/26/2017
# optimize and convert original script to subroutine, DW, 08/29/2018
#
# 'asos.wind' is a date frame that is generated from grab.asos()
# switch to Ben's calc_trajectory, via get.ground.hgt(), DW
#
# fix wd err, if wd.err is closed to -360 or 360, cycle it, DW, 08/31/2018
# add surface wind from ASOS, DW, 09/19/2018 

cal.met.wind.asos <- function(filename, met, met.path, met.format, workdir, 
                              site, timestr, overwrite = F, asos.file, nhrs = -120){

  # loop over each time period
  if (overwrite == T | !file.exists(filename)){

    # read ASOS data 
    asos <- read.table(asos.file, header = T, sep = ',', stringsAsFactors = F) 
    asos <- asos %>% 
        mutate(date = as.POSIXct(valid, format = '%Y-%m-%d %H:%M', tz = 'UTC'), 
                timestr = as.numeric(format(date, '%Y%m%d%H')), 
                temp.asos = as.numeric(tmpc),
                wd.asos = as.numeric(drct),  # degree from true north
                ws.asos = as.numeric(sknt) * 0.514, # convert knots to m/s
                u.asos = sin((wd.asos - 180) * pi/180) * ws.asos,
                v.asos = cos((wd.asos - 180) * pi/180) * ws.asos) %>% 
        dplyr::select(station, lon, lat, date, timestr, temp.asos, 
                      wd.asos, ws.asos, u.asos, v.asos) %>% na.omit()

    # only interpolate model wind for OCO-2 overpasses 
    recp.date <- as.POSIXct(timestr, format = '%Y%m%d%H', tz = 'UTC')
    end.date <- recp.date + nhrs * 60 * 60
    asos <- asos %>% filter(date >= min(end.date, recp.date), date <= max(end.date, recp.date))
    asos.hgt <- 612 # 612 m for two ASOS surface stations
    asos <- asos %>% mutate(met.grdhgt = ifelse(station == 'OERK', 605.1306, 537.8216), 
                            agl = asos.hgt - met.grdhgt)

    # compute unique receptors (time, lat, lon)
    recpstr <- paste(asos$timestr, asos$lat, asos$lon, asos$agl)
    sort.recpstr <- sort(unique(recpstr))
    uni.recp <- matrix(unlist(strsplit(sort.recpstr,' ')), ncol = 4, byrow = T)
    colnames(uni.recp) <- list('time', 'lat', 'lon', 'agl')

    # compute 'receptor'
    receptor <- data.frame(
      lati = as.numeric(uni.recp[, 'lat']),
      long = as.numeric(uni.recp[, 'lon']), 
      zagl = as.numeric(uni.recp[, 'agl'])) %>% mutate(
      run_time = as.POSIXct(uni.recp[, 'time'], '%Y%m%d%H', tz = 'UTC'))

    rundir <- file.path(workdir, paste0('out_wind_', timestr, '_', met, '_asos'))
    #var1 <- c('time', 'index', 'lon', 'lat', 'agl', 'grdht', 'zi', 'temp', 'pres')
    var2 <- c('time', 'indx', 'long', 'lati', 'zagl', 'zsfc', 'mlht', 'temp', 'pres')

    # Ensure necessary files and directory structure are established in the
    # current rundir
    dir.create(rundir, showWarnings = FALSE)

    # change the path of hymodelc executable, DW, 07/31/2018
    # since we added another version of hymodelc (AER_NOAA_branch) under exe
    link_files <- dir(file.path(workdir, 'exe', 'master'))
    if (!file.exists(file.path(rundir, link_files))[1])
      file.symlink(file.path(workdir, 'exe', 'master', link_files), rundir)

    # loop over each unique location + time
    cat(paste(nrow(receptor), 'unique asos station + asos time\n'))
    err.info <- NULL

    for (i in 1 : nrow(receptor)){
      cat(paste('working on', i, 'unique asos\n'))

      sel.asos <- asos[i, ]

      # if no postive AGL
      int.info <- get.ground.hgt(varsiwant = var2, met_loc = met.path,
                                 met_file_format = met.format, n_hours = 1, 
                                 receptor = receptor[i,], rundir = rundir, 
                                 timeout = 20 * 60, r_zagl = receptor$zagl[i], 
                                 run_trajec = T) 
      if (is.null(int.info)) next

      # merge obs and sim
      merge.info <- cbind(sel.asos, int.info) %>%
        rename(u.met = ubar, v.met = vbar, w.met = wbar, temp.met = temp)

      # calculate ws and wd (FROM which dir, degree from true North)
      merge.info <- merge.info %>% mutate(
          ws.met = sqrt(u.met^2 + v.met^2),
          wd.met = atan2(u.met/ws.met, v.met/ws.met) * 180 / pi + 180)

      # when ws == 0, wd is NA, thus, replace NA with 0
      merge.info[is.na(merge.info$wd.met), 'wd.met'] <- 0

      # calculate wind errors
      tmp.err.info <- merge.info %>% mutate(temp.err = temp.met - temp.asos,
                                            u.err  = u.met - u.asos,  
                                            v.err  = v.met - v.asos,
                                            ws.err = ws.met - ws.asos, 
                                            wd.err = wd.met - wd.asos)

      # if wd.err is closed to -360 or 360, cycle it, DW, 08/31/2018
      tmp.err.info[abs(tmp.err.info$wd.err) > 180, 'wd.err'] <-
        abs(tmp.err.info[abs(tmp.err.info$wd.err) > 180, 'wd.err']) - 360

      print(tmp.err.info$wd.err)

      colTF <- F; if (i == 1) colTF <- T
      appTF <- T; if (i == 1) appTF <- F
      write.table(tmp.err.info, file = filename, append = appTF, sep = ',',
                  quote = F, row.names = F, col.names = colTF)
      err.info <- rbind(err.info, tmp.err.info)
    } # end for i

  } else {
    err.info <- read.table(file = filename, sep = ',', header = T)

  } # end if reading txtfile

  return(err.info)
}
# end of function
