### subroutine to define background region with forward-time run from a box
# see 'ggplot.forward.trajec.r' for more info on input variables
# by DW, 11/02/2017

# try box receptors/sources,'dxyp' in trajecmulti(), DW, 11/08/2017
# add time windows for releasing particles, DW, 11/09/2017
# clear things up and add more sites, DW, 04/16/2018

# *** default output directory will be automatically created and under 'xstilt_wd',
# in different Copies; make sure different nummodels are used for different
# overpass events, DW, 07/31/2018
# add customized data filtering, DW, 08/20/2018
# remove data filtering, always use QF = 0 for background estimates, DW, 10/31/2018
# rename output folder names, DW, 01/25/2019 


# ---------------------- function inside run.forward.trajec ------------------ #
# can work on multiple orbits, DW, 04/30/2021
convert.timestr2ident = function(timestr, dtime, site_lon, site_lat, agl, numpar, dxyp) {
  
  # get lat/lon for city center as well as time string
  clon = signif(site_lon, 6)
  clat = signif(site_lat, 6)
  lonstr = ifelse(clon > 0, paste0(clon, 'E'), paste0(abs(clon), 'W'))
  latstr = ifelse(clat > 0, paste0(clat, 'N'), paste0(abs(clat), 'S'))
  recp.date = as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')

  # add an additional half an hour to make sure we capture the overpass hour
  # DW, 03/04/2021
  # modify when multiple orbits are found in one day, DW, 04/30/2021
  date.list = do.call(list, lapply(recp.date, FUN = function(x) 
  { 
    rel.date = x + c(dtime, 0.5) * 60 * 60   # in second

    # calculate each release times when continuously releasing particles
    format.date = strsplit.to.df(format(rel.date, '%Y-%m-%d-%H-%M-%S'), sep = '-') 
    format.date = format.date %>% mutate_all(funs(as.numeric), colnames(format.date))
    colnames(format.date) = c('yr', 'mon', 'day', 'hr', 'min', 'sec')
    format.date$lat = clat; format.date$lon = clon 

    # create traj names for each release time, in form of
    # YYYYxMMxDDxHHxmmxNxWxmAGLxPxdeg that works with Trajecmulti()
    format.date$ident = paste0(paste(format(rel.date, '%Yx%mx%dx%Hx%M'), 
                               latstr, lonstr, formatC(agl, width = 5, flag = 0), 
                               paste0(numpar, 'P'), dxyp, sep = 'x'), '.rds')
    return(format.date)
  }))
 
  return(date.list)
} # end of convert.timestr2ident



# ---------------------- function inside run.forward.trajec ------------------ #
find.all.metfiles = function(timestr, dtime, met_file_format, met_path, nhrs) {

  # get met files, run_time in form of 'yyyy-mm-dd HH:MM:SS (UTC)'
  recp.date = as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
  rel.date = recp.date + dtime * 60 * 60  # in second

  met.files = NULL
  for (r in 1:length(rel.date))
    met.files = c(met.files, find_met_files(t_start = rel.date[r],
                                             met_file_format = met_file_format,
                                             n_hours = nhrs, met_path = met_path))
  met.files = basename(unique(met.files))
  if (length(met.files) == 0) {
    cat('find.all.metfiles(): No meteo fields found...please check\n'); return()}

  return(met.files)
} # end of find.all.met.files



# ---------------------- run.forward.trajec ()  ----------------------------- #
# timestr in YYYYMMDD if provided with path for satellite data 
# otherwise timestr should have format of YYYYMMDDHH
run.forward.trajec = function(site, site_lon, site_lat, timestr, 
                              run_trajec = T, run_hor_err = F, 
                              run_wind_err = F, run_ver_err = F, xstilt_wd, 

                              # params for generating trajec
                              store_path, box.len, dtime.from, 
                              dtime.to, dtime.sep, nhrs, delt, agl, numpar, 
                               
                              # met fields
                              met, met_res, met_file_format, met_path, 

                              # horizontal and vertial trans error input
                              # zi scaling factor
                              raob_path, siguverr = NA, nhrs.zisf = 24, zisf, 
                               
                              # add a few params to determine the overpass hour
                              # if satellite paths are provided, default = NA
                              # sensor_var is needed if grabbing data from OCO
                              sensor = c(NA, 'OCO', 'TROPOMI')[1], 
                              sensor_path = NA, sensor_ver = NA, sensor_gas = NA){
  

  setwd(xstilt_wd); source('r/dependencies.r', local = T) # source all functions
  cat(paste('\n\n## ----- Working on', site, 'on', timestr, '----- ##\n'))

  # add a portion to get receptor hour based on satellite overpass hour,
  # if satellite files are given, DW, 03/04/2021
  err_message = 'run.forward.trajec(): NO observation found...\n'
  traj_path = file.path(store_path, 'out_forward')

  if ( nchar(timestr) != 10 ) {
    
    if ( !is.na(sensor_path) ) {
      cat(paste('run.forward.trajec(): since no hour is provided,',
                'figuring out the overpass hour based on satellite data on', 
                timestr, '; it takes a while\n'))
      
      # define spatial domain for grabbing satellite observations
      lon_lat = data.frame(citylon = site_lon, citylat = site_lat, 
                           minlon = site_lon - 1, maxlon = site_lon + 1, 
                           minlat = site_lat - 1, maxlat = site_lat + 1)

      obs_df = get.sensor.obs(site, timestr, sensor, sensor_gas, sensor_fn = NULL,
                              sensor_path, sensor_ver, qfTF = F, tropomi_qa = 0,
                              lon_lat)
      
      if ( is.null(obs_df) ) return()
      if ( grepl('OCO', sensor) ) obs_df = obs_df %>% rename(time_utc = timestr)

      # now in YYYYMMDDHH
      timestr = unique(substr(obs_df$time_utc, 1, 10)); print(timestr)  
      if (length(timestr) > 1) cat('run.forward.trajec(): Found multiple orbits per day\n')

    } else stop('run.forward.trajec(): Invalid satellite sensor...\n')
  }   # end if 
  # stop(paste('Incorrect length of @param timestr (YYYYMMDDHH)',
  #            'and no satellite paths available. Please check...\n'))


  # ------------------------------ STEP 3 --------------------------------- #
  ## get horizontal transport error component if run_hor_err = T
  # if run_wind_err = T, prepare RAOB data and compute model-data wind comparisons
  # *** if you do not want to run wind error analysis, set run_wind_err to FALSE 
  # and prescribe a horizontal wind error, e.g., siguverr = 3 (with unit of m/s)
  # for more info, please see get.uverr.r
  
  # define spatial domain for calculating wind error or grabbing overpass hours
  err.lon.lat = data.frame(citylon = site_lon, citylat = site_lat, 
                           minlon = site_lon - 5, maxlon = site_lon + 5, 
                           minlat = site_lat - 5, maxlat = site_lat + 5)

  # determine dxpy using the desired box size and grid length of met fields, 
  # since final box length = 2 * dxyp * met_res, DW, 06/30/2020
  dxyp = box.len / met_res / 2   

  # reformat trajec info that fits STILT version 1, this script needs to updated later
  dtime   = seq(dtime.from, dtime.to, dtime.sep)     # vector in hours
  results = convert.timestr2ident(timestr, dtime, site_lon, site_lat, agl, numpar, dxyp)

  # parasm for getting error stats ----------------------------------
  namelist = list(run_hor_err = run_hor_err, run_ver_err = run_ver_err, 
                  run_wind_err = run_wind_err, zisf = zisf, odiac_path = NULL, 
                  raob_path = raob_path, nhrs = nhrs, met = met, 
                  met_path = met_path, met_file_format = met_file_format, 
                  agl = c(0, 100), store_path = store_path)

  # if running trajec
  if (run_trajec) {

    # in case there are multiple orbits in a single day, DW, 05/01/2021
    for (tt in 1 : length(timestr)) {
      
      # path that will store the generated wind error statistics if run_hor_err = T
      errlist = config_trans_err(namelist, site, lon_lat = err.lon.lat, 
                                 timestr = timestr[[tt]], xstilt_wd)
      hor_err = errlist$hor_err 
      pbl_err = errlist$pbl_err 
      date_df = results[[tt]]

      # create out_forward dir for generating forward trajec
      traj_dir = file.path(traj_path, paste0('out_forward_', timestr[[tt]], '_', site))
      if (!is.na(sensor)) traj_dir = paste0(traj_dir, '_', sensor); print(traj_dir)

      # delete previous directories and then create new one
      system(paste0('rm -rf ', traj_dir), ignore.stderr = T)
      dir.create(traj_dir, showWarnings = FALSE, recursive = T)

      # linking AER_NOAA_branch's hymodelc and other executables to outpath
      exes = list.files(file.path(xstilt_wd, 'exe'))
      file.symlink(file.path(xstilt_wd, 'exe', exes), traj_dir)
      
      # if using multiple receptors, or box of receptors or sources, turn it on,
      # then call updated Trajecmulti() instead of Trajec()
      varstrajec = c('time', 'index', 'lat', 'lon', 'agl', 'grdht', 'foot',
                     'sampt', 'dmass', 'zi', 'pres')  # trajec output variables

      #### try box receptors/sources, DW, 11/08/2017
      # the updated Trajecmulti() and fortran codes will randomly place receptors
      # according to dxyp
      cat('run.forward.trajec(): Generating forward trajec...\n')
      met_fns = find.all.metfiles(timestr[[tt]], dtime, met_file_format, met_path, nhrs)
      
      Trajecmulti(yr = date_df$yr - 2000, mon = date_df$mon, day = date_df$day, 
                  hr = date_df$hr, mn  = date_df$min, outname = date_df$ident, 
                  numpar = numpar, lat = date_df$lat, lon = date_df$lon,
                  dxyp = rep(dxyp, nrow(date_df)), 
                  dzp  = rep(0, nrow(date_df)),
                  agl  = rep(agl, nrow(date_df)), 
                  nhrs = rep(nhrs, nrow(date_df)),
                  nummodel = timestr[[tt]], 
                  metd = c('fnl', 'awrf'), 
                  outpath = traj_dir, 
                  overwrite = run_trajec,
                  metfile = met_fns, metlib = paste0(met_path, '/'),
                  doublefiles = T, 
                  rundir = dirname(traj_dir), 
                  rundirname = basename(traj_dir), 
                  varsout = varstrajec, 
                  siguverr = hor_err$siguverr, 
                  TLuverr = hor_err$TLuverr, 
                  zcoruverr = hor_err$zcoruverr, 
                  horcoruverr = hor_err$horcoruverr,
                  hymodelc.exe = './hymodelc.aer', # use the AER version of hymodelc
                  setup.list = list(DELT = delt, VEGHT = 0.5)) %>% invisible()
    }  # end for tt

  } else cat('run.forward.trajec(): do nothing, check run_* flags...\n')

} # end of script


