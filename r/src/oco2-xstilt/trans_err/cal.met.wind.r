# script to interpolate the winds from GDAS, using trajwind() for each receptor based on raobiosonde data
# read output model-data comparisons and make plots
# written by DW, 05/15/2017

# Add temporal interpolation of GDAS winds, DW, 05/26/2017
# optimize and convert original script to subroutine, DW, 08/29/2018
#
# 'raob.wind' is a date frame that is generated from grab.raob()
# switch to Ben's calc_trajectory, via get.ground.hgt(), DW
#
# fix wd err, if wd.err is closed to -360 or 360, cycle it, DW, 08/31/2018
# add surface wind from ASOS, DW, 09/19/2018 

cal.met.wind <- function(filename, met, met.path, met.format, met.files = NULL,
  workdir, site, timestr, overwrite = F, raob.path = NULL, nhrs = -120,
  raob.format = 'fsl'){

  # loop over each time period
  if (overwrite == T | !file.exists(filename)){

    # grab RAOB
    raob <- grab.raob(raob.path, timestr, workdir, nhrs, format = raob.format,
      overwrite) %>% dplyr::rename(temp.raob = temp, u.raob = u, v.raob = v,
      ws.raob = ws, wd.raob = wd)

    # compute unique receptors (time, lat, lon)
    recpstr <- paste(raob$timestr, raob$lat, raob$lon)
    sort.recpstr <- sort(unique(recpstr))
    uni.recp <- matrix(unlist(strsplit(sort.recpstr,' ')), ncol = 3, byrow = T)
    colnames(uni.recp) <- list('time', 'lat', 'lon')

    # compute 'receptor'
    receptor <- data.frame(
      lati = as.numeric(uni.recp[, 'lat']),
      long = as.numeric(uni.recp[, 'lon'])) %>% mutate(
      run_time = as.POSIXct(uni.recp[, 'time'], '%Y%m%d%H', tz = 'UTC'))

    rundir <- file.path(workdir, paste0('out_wind_', timestr, '_', met))
    var1 <- c('time', 'index', 'lon', 'lat', 'agl', 'grdht', 'zi', 'temp', 'pres')
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
    cat(paste(nrow(receptor), 'unique raob station + raob time\n'))
    err.info <- NULL

    for (i in 1:nrow(receptor)){
      cat(paste('working on', i, 'unique raob\n'))

      # obtain ground hgts, if no rds found, will start to generate trajec
      tmp.info <- get.ground.hgt(varsiwant = var2, met_loc = met.path,
        met_file_format = met.format, n_hours = 1, receptor = receptor[i,],
        rundir = rundir, timeout = 10 * 60, met_files = met.files, 
        run_trajec = T)
      grdhgt <- tmp.info$zsfc

      # METHOD 2 to int winds --
      # cal AGL based on raobiosonde sounding hgt (as ASL) and
      # interpolated ground HGT (from GDAS)
      sel.raob <- raob[recpstr == sort.recpstr[i],]
      sel.agl  <- sel.raob$hgt - grdhgt

      ## !!! get rid of negative agl, also max GDAS level only goes up to 25km
      # if above remove!!!
      pos.agl  <- sel.agl [sel.agl > 0 & sel.agl < 25000]
      pos.raob <- sel.raob[sel.agl > 0 & sel.agl < 25000,]

      # if no postive AGL
      if (length(pos.agl) == 0) next
      int.info <- get.ground.hgt(varsiwant = var2, met_loc = met.path,
        met_file_format = met.format, n_hours = 1, receptor = receptor[i,],
        rundir = rundir, timeout = 20 * 60, r_zagl = pos.agl, run_trajec = T) %>%
        mutate(zagl = pos.agl)

      # merge obs and sim
      merge.info <- cbind(pos.raob, int.info) %>%
        rename(u.met = ubar, v.met = vbar, w.met = wbar, temp.met = temp)

      # calculate ws and wd (FROM which dir, degree from true North)
      merge.info <- merge.info %>% mutate(
          ws.met = sqrt(u.met^2 + v.met^2),
          wd.met = atan2(u.met/ws.met, v.met/ws.met) * 180 / pi + 180)

      # when ws == 0, wd is NA, thus, replace NA with 0
      merge.info[is.na(merge.info$wd.met), 'wd.met'] <- 0

      # calculate wind errors
      tmp.err.info <- merge.info %>% mutate(temp.err = temp.met - temp.raob,
        u.err  = u.met  - u.raob,  v.err  = v.met  - v.raob,
        ws.err = ws.met - ws.raob, wd.err = wd.met - wd.raob)

      # if wd.err is closed to -360 or 360, cycle it, DW, 08/31/2018
      tmp.err.info[abs(tmp.err.info$wd.err) > 180, 'wd.err'] <-
        abs(tmp.err.info[abs(tmp.err.info$wd.err) > 180, 'wd.err']) - 360

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
