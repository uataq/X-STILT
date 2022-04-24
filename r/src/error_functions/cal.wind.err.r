# script to interpolate the winds from GDAS, using trajwind() 
# for each receptor based on raobiosonde data
# read output model-data comparisons and make plots
# written by DW, 05/15/2017

# Add temporal interpolation of GDAS winds, DW, 05/26/2017
# optimize and convert original script to subroutine, DW, 08/29/2018
#
# 'raob.wind' is a date frame that is generated from grab.raob()
# switch to Ben's calc_trajectory, via get.met.varh(), DW
#
# fix wd err, if wd.err is closed to -360 or 360, cycle it, DW, 08/31/2018
# add surface wind from ASOS, DW, see cal.met.wind.asos(), 09/19/2018 
# minor bug occur when writing table into a txt file, DW, 12/04/2018
# update the code with latest changes to STILT, DW, 02/27/2021 

#' add max height @param maxagl for calculating wind error, default is 6 km

cal.wind.err = function(err_file, met, met_path, met_file_format, xstilt_wd, 
                        err_path, site, timestr, simstep_namelist = NULL, 
                        raob_fn = NULL, raob.nhrs = -24, maxagl = 6000){

  # loop over each time period
  # grab RAOB and rename
  raob = grab.raob(raob_fn, timestr, err_path, nhrs = raob.nhrs, 
                   overwriteTF = T) %>%
         rename(temp.raob = temp, u.raob = u, v.raob = v, 
                ws.raob = ws, wd.raob = wd)
  
  # compute unique receptors (time, lat, lon)
  recpstr = paste(raob$timestr, raob$lat, raob$lon)
  sort.recpstr = sort(unique(recpstr))
  uni.recp = matrix(unlist(strsplit(sort.recpstr,' ')), ncol = 3, byrow = T)
  colnames(uni.recp) = list('time', 'lat', 'lon')

  # compute receptor
  recp.df = data.frame(lati = as.numeric(uni.recp[, 'lat']),
                       long = as.numeric(uni.recp[, 'lon'])) %>% 
            mutate(run_time = as.POSIXct(uni.recp[, 'time'], 'UTC', '%Y%m%d%H'))

  # ------------------------------------------------------------------------
  # *** temz for temp at certain height; temp for temp at the lowest height 
  # Do not change the following input variables for this simulation!!
  # only release one particle, but turn on turbulance 'nturb'
  rundir = file.path(err_path, paste0('out_wind_', timestr, '_', met))
  var2 = c('time', 'indx', 'long', 'lati', 'zagl','zsfc','mlht','temz','pres')
  simstep_namelist$delt = 1
  simstep_namelist$nturb = T
  simstep_namelist$numpar = 1
  simstep_namelist$varsiwant = var2

  # Ensure necessary files and directory structure are established in the
  # current rundir
  dir.create(rundir, showWarnings = FALSE)

  # change the path of hymodelc executable, DW, 07/31/2018
  # since we added another version of hymodelc (AER_NOAA_branch) under exe
  # reorganize HYSPLIT dependencies, DW, 02/19/2019
  link_files = dir(file.path(xstilt_wd, 'exe'))
  if (!file.exists(file.path(rundir, link_files))[1])
      file.symlink(file.path(xstilt_wd, 'exe', link_files), rundir)
  # ------------------------------------------------------------------------

  # write header to txt file that stores wind error statistics
  header = c('timestr', 'lat', 'lon', 'elev', 'pres', 'hgt', 'temp.raob', 
             'ws.raob', 'wd.raob', 'u.raob', 'v.raob', 'u.met', 'v.met', 
             'w.met', 'zsfc', 'pres.met', 'temp.met', 'ws.met', 'wd.met', 
             'temp.err', 'u.err', 'v.err', 'ws.err', 'wd.err')     
  write(header, file = err_file, append = F,sep = ',',ncolumns = length(header))

  # loop over each unique location + time
  cat(paste('\nTotal', nrow(recp.df), 'unique raob station & raob time\n'))
  err.info = NULL
  for (i in 1 : nrow(recp.df)) {

    cat(paste('# ----- working on #', i, 'unique raob loc & time ---- #\n'))
    # no need to calculate footprint per trajec, thus turn off hnf_plume 
    tmp.info = get.met.varh(receptor = recp.df[i,], agl = 1, run_trajec = T, 
                            namelist = simstep_namelist, rundir = rundir, 
                            hnf_plume = F, met_file_format, met_path, 
                            n_hours = 1, timeout = 5 * 60)

    # add error message, DW, 11/16/2018
    if (is.null(tmp.info)) { 
      cat('NO met vars extracted (e.g., recp outside met fields)...\n'); next}
    grdhgt = tmp.info$zsfc

    # METHOD 2 to int winds --
    # cal AGL based on raobiosonde sounding hgt (as ASL) and
    # interpolated ground HGT (from GDAS)
    sel.raob = raob[recpstr == sort.recpstr[i],]
    #sel.agl = sel.raob$hgt - grdhgt

    ## get rid of negative agl; max GDAS level only goes up to maxagl 
    # default value for maxagl is 6 km (since we care most about PBL)
    pos.raob = sel.raob %>% mutate(zagl = hgt - grdhgt) %>% 
               filter(zagl > 0 & zagl <= maxagl)
    pos.zagl = pos.raob$zagl 
    
    # if no postive AGL, next station
    if (length(pos.zagl) == 0) next

    # otherwise extract wind info from met fields using STILT
    int.info = NULL 
    for (agl in pos.zagl) {
      
      tmp.int = get.met.varh(receptor = recp.df[i,], agl = agl, 
                             run_trajec = T, namelist = simstep_namelist, 
                             rundir = rundir, hnf_plume = F, met_file_format, 
                             met_path, n_hours = 1, timeout = 5 * 60) 
      #print(tmp.int)
      # add error message, DW, 11/16/2018
      if (is.null(tmp.int)) { 
        cat('NO met vars extracted (e.g., recp outside met fields)...\n')
        next 
      }

      int.info = rbind(int.info, tmp.int %>% mutate(zagl = agl)) 
    } # end for

    # merge obs and sim
    merge.info = left_join(pos.raob, int.info %>% rename(pres.sim = pres), 
                           by = 'zagl') %>% 
                 rename(u.met = ubar, v.met = vbar, w.met = wbar, 
                        temp.met = temz) %>%

                 # calculate ws and wd (FROM which dir, degree from true North)
                 mutate(ws.met = sqrt(u.met^2 + v.met^2),
                        wd.met = atan2(u.met/ws.met, v.met/ws.met)*180/pi + 180)

    # when ws == 0, wd is NA, thus, replace NA with 0
    merge.info[is.na(merge.info$wd.met), 'wd.met'] = 0

    # calculate wind errors
    tmp.err.info = merge.info %>% mutate(temp.err = temp.met - temp.raob, 
                                         u.err = u.met - u.raob,  
                                         v.err = v.met - v.raob, 
                                         ws.err = ws.met - ws.raob, 
                                         wd.err = wd.met - wd.raob)

    # if wd.err is closed to -360 or 360, cycle it, DW, 08/31/2018
    tmp.err.info[abs(tmp.err.info$wd.err) > 180, 'wd.err'] =
      abs(tmp.err.info[abs(tmp.err.info$wd.err) > 180, 'wd.err']) - 360

    # minor bug occur when writing table, DW, 12/04/2018
    write.table(tmp.err.info, file = err_file, append = T, sep = ',',
                quote = F, row.names = F, col.names = F)
    
    err.info = rbind(err.info, tmp.err.info)
  } # end for i

  
  return(err.info)
}

# end of function
