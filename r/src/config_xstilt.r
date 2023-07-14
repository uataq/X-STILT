# STILT R Executable
# For documentation, see https://uataq.github.io/stilt/
# Ben Fasoli
# edited by Dien Wu for X-STILT


# preserve original 'run_stilt.r' and modify 'run_stilt.r' as a subroutine
# add ziscale, DW, 07/25/2018
# add two before_* functions for getting ground height and AK PW weighting, DW, 01/24/2019
# allows for generate footprints with different resolutions, DW, 02/11/2019 
# add additional variables due to changes in simulation_step()'s argument, DW, 07/02/2020
# add/remove variables as migrating to STILT-HYSPLIT, DW, 07/03/2020 
# refactoring: move some of the calculation in the main script to this function
# DW, 10/26/2020

config_xstilt = function(namelist){

  # User inputs ----------------------------------------------------------------
  xstilt_wd = namelist$xstilt_wd
  lib.loc   = .libPaths()[1]
  site      = namelist$site
  timestr   = namelist$timestr[[1]]
  lon_lat   = namelist$lon_lat[[1]]
  obs_path  = namelist$obs_path
  obs_fn    = namelist$obs_fn
  obs_sensor = namelist$obs_sensor
  obs_species = unlist(namelist$obs_species)
  store_path = namelist$store_path

  # make sure obs_* are set to NA and no AK weighting for ideal runs
  if (is.na(obs_sensor)) { 
    obs_path = obs_species = NA
    namelist$ak_wgt = F 

    # *** For ideal simulations without satellite data, set to NA or FALSE 
    # receptors will be placed based on lati/long info from receptor_demo.csv
    num_bg_lat = num_bg_lon = num_nf_lat = num_nf_lon = num_jitter = NA
    jitterTF = FALSE
  }

  if ( obs_sensor == 'TROPOMI' & length(obs_species) > 1) 
    stop('Please only choose one atmospheric species per time for TROPOMI, denoted by @param obs_species...\n')


  # Model control -------------------------------------------------------------
  rm_dat      = T
  timeout     = namelist$timeout  # in sec
  run_slant   = namelist$run_slant 
  run_trajec  = namelist$run_trajec
  run_foot    = namelist$run_foot
  run_sim     = namelist$run_sim
  run_hor_err = namelist$run_hor_err
  run_ver_err = namelist$run_ver_err
  run_wind_err = namelist$run_wind_err
  
  cat('\n --- Configuring --- \n')
  if (run_trajec) cat('Need to generate trajec...\n')
  if (run_foot)   cat('Need to generate footprint...\n')
  if (!run_trajec & !run_foot & !run_sim) 
    cat('NO calculations will be performed, please check run_* flags..\n')
  if (run_sim) { run_trajec = F; run_foot = F }   

  # if running trajec/foot, no sim allowed
  if (run_trajec | run_foot) run_sim = F   

  # line source for agl
  agl     = c(namelist$minagl, namelist$maxagl) 
  n_hours = namelist$nhrs
  numpar  = namelist$numpar

  varsiwant = namelist$varstrajec
  if (length(varsiwant) == 0)
    varsiwant = c('time', 'indx', 'long', 'lati', 'zagl', 'zsfc', 'foot', 
                  'mlht', 'dens', 'samt', 'sigw', 'tlgr', 'temp', 'pres')


  # met fields ----------------------------------------------------------------
  met       = namelist$met
  met_path  = namelist$met_path
  n_met_min = namelist$n_met_min
  met_file_format = namelist$met_file_format
  met_subgrid_buffer = namelist$met_subgrid_buffer
  met_subgrid_enable = namelist$met_subgrid_enable
  met_subgrid_levels = namelist$met_subgrid_levels

  # Footprint params ----------------------------------------------------------
  xmn  = round(lon_lat$site_lon) - namelist$foot_dlon
  xmx  = round(lon_lat$site_lon) + namelist$foot_dlon
  ymn  = round(lon_lat$site_lat) - namelist$foot_dlat
  ymx  = round(lon_lat$site_lat) + namelist$foot_dlat
  xres = namelist$foot_res
  yres = namelist$foot_res
  xres2 = namelist$foot_res2
  yres2 = namelist$foot_res2
  time_integrate  = namelist$time_integrate
  time_integrate2 = namelist$time_integrate2
  hnf_plume       = namelist$hnf_plume
  smooth_factor   = namelist$smooth_factor
  projection      = namelist$projection 
  

  # Create output directory ------------------------------------------------
  outlist = create.outwd(timestr, obs_species, obs_sensor, obs_path, obs_fn, 
                         lon_lat, store_path, met, run_hor_err)
  obs_info = outlist$obs_info         # get obs file name if satellite is used
  output_wd = outlist$output_wd

  if ( length(unique(obs_info$fn)) > 1 ) {  # if multiple tracks per day
    track_indx = which(obs_info$tot.count == max(obs_info$tot.count))
    obs_fn = obs_info$fn[track_indx]
    if (length(output_wd) > 1) output_wd = output_wd[track_indx]

  } else obs_fn = unique(obs_info$fn)
  cat(paste('Obs file -', obs_fn, '\n'))


  # Compute XCO2 using existing traj & foot ------------------------------------
  # calculate XCO2 concentration and its error (need trajec and footprint ready)
  if ( !run_trajec & !run_foot & run_sim & 'CO2' %in% obs_species ) {

    # ------------------------  Horizontal trans error ---------------------- #
    ### simulate transport error in XCO2 due to met errors, DW, 07/25/2018
    # requires two sets of trajectories before running the following section:
    if (run_hor_err) { # this does not need footprint

      ## call cal.trans.err() to estimate receptor-level trans err [ppm]
      # get actual ppm error, need to have error statistics ready
      # see cal.trajfoot.stat() in called before_footprint_xstilt.r for err stat
      cat('Start simulations of XCO2 error due to horizontal trans err...\n')
      result = cal.trans.err(site, timestr, workdir = xstilt_wd, 
                             outdir = output_wd, store_path, met)
      if (is.null(result)) 
        stop('No results calculated, check cal.trans.err()\n')

    } else {

      # ------------------- XCO2 or emiss error simulations ----------------- #
      # requires trajec/footprints ready for running sim, DW, 06/04/2018  
      # call func to match ODIAC emissions with xfoot & sum up to get 'dxco2.ff'
      cat('Start simulations of XCO2.ff or its error due to emiss err...\n')
      result = run.xco2ff.sim(site, timestr, vname = namelist$odiac_ver, 
                              emiss.path = namelist$odiac_path, 
                              output_wd, foot.res = xres, xstilt_wd, store_path,
                              nhrs = n_hours, obs_sensor, met, 
                              run_emiss_err = namelist$run_emiss_err, 
                              edgar.file = namelist$edgar_file, 
                              ffdas.file = namelist$ffdas_file)
      
      if (is.null(result)) stop('No results calculated...\n')
      return(result)
    } # end if run_hor_err

  } else if ( run_sim & !'CO2' %in% obs_species ) 
    stop('You will need to write your own scripts in coupling non-CO2 emissions with footprint...\n\n')
  # end if run_sim


  # Receptor locations & time -------------------------------------------------
  # IF for ideal simulation, read from receptor file 
  #if ( is.na(obs_sensor) ) {
  if ( !is.na(namelist$recp_fn) ) {
    receptors = read.table(namelist$recp_fn, header = T, sep = ',')
    if ( 'time' %in% colnames(receptors) ) {    # if there is a `time` column
      time_string = receptors$time 
      receptors$time = NULL
    } else if ( !is.na(namelist$timestr) ) {  # no time column but timestr ava
      time_string = namelist$timestr
    } else stop('Missing receptor time strings...please check @param recp_fn or @param timestr\n')

    time_string = as.character(time_string)
    nchart = nchar(time_string[1])
    if (nchart == 8 ) formatt = '%Y%m%d'
    if (nchart == 10) formatt = '%Y%m%d%H'
    if (nchart == 12) formatt = '%Y%m%d%H%M'
    if (nchart == 14) formatt = '%Y%m%d%H%M%S'
    if (!nchart %in% c(8, 10, 12, 14)) 
      stop('Incorrect form of time string...please check @param recp_fn or @param timestr\n')

    receptors$run_time = as.POSIXct(time_string, 'UTC', format = formatt) 
    receptors$zagl = list(agl)
    
  } else {    
    
    # IF for simulations using satellite or ground-based X data
    # obtain receptors' locations based on obs availability
    receptors = get.recp.sensorv2(timestr, 
                                  obs_filter = unlist(namelist$obs_filter), 
                                  obs_fn, obs_sensor, obs_path, obs_species, 
                                  lon_lat, jitterTF = namelist$jitterTF, 
                                  num_jitter = namelist$num_jitter, numpar, 
                                  nf_dlat = namelist$nf_dlat, 
                                  nf_dlon = namelist$nf_dlon, 
                                  num_bg_lat = namelist$num_bg_lat, 
                                  num_bg_lon = namelist$num_bg_lon, 
                                  num_nf_lat = namelist$num_nf_lat, 
                                  num_nf_lon = namelist$num_nf_lon, 
                                  agl, run_trajec, run_slant, output_wd)

  } # end if
  
  cat(paste('Done with receptor setup...total', nrow(receptors), 
            'receptor(s)..\n\n'))

  # if for slant column release, be explicit about the output dir names
  if ( is.list(receptors$lati) & is.list(receptors$long) ) {
    sim_id_format = paste0('%Y%m%d%H%M_', receptors$long0, '_', receptors$lati0,
                           '_', ifelse(length(receptors$zagl[[1]]) > 1, 'X', r_zagl))
    simulation_id = strftime(receptors$run_time, sim_id_format, 'UTC')
  } else simulation_id = NA


  # Transport and dispersion settings, use default setting in STILT-R v2
  capemin     = -1
  cmass       = 0
  conage      = 48
  cpack       = 1
  delt        = 2
  dxf         = 1
  dyf         = 1
  dzf         = 0.01
  efile       = ''
  emisshrs    = 0.01
  frhmax      = 3
  frhs        = 1
  frme        = 0.1
  frmr        = 0
  frts        = 0.1
  frvs        = 0.01
  hscale      = 10800
  ichem       = 8
  idsp        = 2
  initd       = 0
  k10m        = 1
  kagl        = 1
  kbls        = 1; if (toupper(met) == 'NARR') kbls = 2
  kblt        = 5
  kdef        = 0
  khinp       = 0
  khmax       = 9999
  kmix0       = 150
  kmixd       = 3
  kmsl        = 0
  kpuff       = 0
  krand       = 4
  krnd        = 6
  kspl        = 1
  kwet        = 1
  kzmix       = 0
  maxdim      = 1
  maxpar      = numpar
  mgmin       = 10
  mhrs        = 9999
  nbptyp      = 1
  ncycl       = 0
  ndump       = 0
  ninit       = 1
  nstr        = 0
  nturb       = 0
  nver        = 0
  outdt       = 0
  p10f        = 1
  pinbc       = ''
  pinpf       = ''
  poutf       = ''
  qcycle      = 0
  rhb         = 80
  rht         = 60
  splitf      = 1
  tkerd       = 0.18
  tkern       = 0.18
  tlfrac      = 0.1
  tout        = 0
  tratio      = 0.75
  tvmix       = 1
  veght       = 0.5
  vscale      = 200
  vscaleu     = 200
  vscales     = -1
  wbbh        = 0
  wbwf        = 0
  wbwr        = 0
  wvert       = FALSE
  w_option    = 0
  zicontroltf = 0
  z_top       = 25000
  if (toupper(met) == 'NARR') z_top = 15000 
  if (toupper(met) == 'HRRR') z_top = 20000

  # Aggregate STILT/HYSPLIT namelist
  simstep_namelist = list(capemin = capemin, cmass = cmass, conage = conage,
                          cpack = cpack, delt = delt, dxf = dxf, dyf = dyf,
                          dzf = dzf, efile = efile, frhmax = frhmax, 
                          frhs = frhs, frme = frme, frmr = frmr, frts = frts, 
                          frvs = frvs, hnf_plume = hnf_plume, hscale = hscale, 
                          ichem = ichem, idsp = idsp, initd = initd, 
                          k10m = k10m, kagl = kagl, kbls = kbls, kblt = kblt, 
                          kdef = kdef, khinp = khinp, khmax = khmax, 
                          kmix0 = kmix0, kmixd = kmixd, kmsl = kmsl,
                          kpuff = kpuff, krand = krand, krnd = krnd, 
                          kspl = kspl, kwet = kwet, kzmix = kzmix, 
                          maxdim = maxdim, maxpar = maxpar, mgmin = mgmin, 
                          ncycl = ncycl, ndump = ndump, ninit = ninit,
                          nstr = nstr, nturb = nturb, numpar = numpar, 
                          nver = nver, outdt = outdt, p10f = p10f, 
                          pinbc = pinbc, pinpf = pinpf, poutf = poutf, 
                          qcycle = qcycle, rhb = rhb, rht = rht,
                          splitf = splitf, tkerd = tkerd, tkern = tkern, 
                          tlfrac = tlfrac, tout = tout, tratio = tratio, 
                          tvmix = tvmix, varsiwant = varsiwant, veght = veght,
                          vscale = vscale, vscaleu = vscaleu, 
                          vscales = vscales, wbbh = wbbh, wbwf = wbwf, 
                          wbwr = wbwr, winderrtf = 0, wvert = wvert, 
                          zicontroltf = zicontroltf)
  
  # Get error stats if needed -------------------------------------------------
  # Calculating error stats for horizontal and vertical transport errors
  errlist = config_trans_err(namelist, site, lon_lat, timestr, xstilt_wd, 
                             simstep_namelist)
  
  hor_err  = errlist$hor_err 
  pbl_err  = errlist$pbl_err 
  emiss_fn = errlist$emiss_fn 
  horcoruverr = hor_err$horcoruverr
  siguverr    = hor_err$siguverr
  tluverr     = hor_err$TLuverr
  zcoruverr   = hor_err$zcoruverr
  horcorzierr = pbl_err$horcorzierr
  sigzierr    = pbl_err$sigzierr
  tlzierr     = pbl_err$TLzierr
  zicontroltf = pbl_err$zicontroltf
  ziscale     = pbl_err$ziscale         # prescribe PBL scaling, vector form
  #ziscale    = rep(list(rep(1, 24)), nrow(receptors))

  # customized functions
  before_footprint = before_footprint_xstilt
  before_trajec    = before_trajec_xstilt

  # Parallel simulation settings
  slurm   = namelist$slurm
  n_cores = namelist$n_cores
  n_nodes = namelist$n_nodes

  # OR specify the ammount of memory per cpu your job requires  
  slurm_options = list(time = namelist$job_time, 
                       account = namelist$slurm_account, 
                       partition = namelist$slurm_partition)
  if (!is.null(namelist$mem_per_node)) slurm_options$mem = namelist$mem_per_node
  
  # specify job names
  jobname = paste('XSTILT', site, timestr, obs_sensor, obs_species, sep = '_')
  
  if ( is.na(obs_sensor) )      # no obs involved - meaning no AK 
    jobname = paste('XSTILT', site, timestr, 'ideal', sep = '_')
  
  if ( length(timestr) > 1 )    # multiple receptor times
    jobname = paste('XSTILT', site, min(timestr), max(timestr), 
                    obs_sensor, obs_species, sep = '_')

  if ( length(obs_species) > 1 )    # multiple receptor times
    jobname = paste('XSTILT', site, timestr, obs_sensor, 'multi', sep = '_')
  
  if ( length(timestr) > 1 & length(obs_species) > 1 )
    jobname = paste('XSTILT', site, min(timestr), max(timestr), 
                    obs_sensor, 'multi', sep = '_')

  if (run_hor_err | run_ver_err) jobname = paste0(jobname, '_error')


  # Startup messages -----------------------------------------------------------
  message('\n\nInitializing X-STILT')
  message('Number of receptors: ', nrow(receptors))
  message('Number of parallel threads: ', n_nodes * n_cores)


  # Source dependencies --------------------------------------------------------
  setwd(xstilt_wd)
  source('r/dependencies.r')


  # Structure out directory ----------------------------------------------------
  # Outputs are organized in three formats. by-id contains simulation files by
  # unique simulation identifier. particles and footprints contain symbolic links to the particle trajectory and footprint files in by-id
  system(paste0('rm -r ', output_wd, '/footprints'), ignore.stderr = T)
  if (run_trajec) {
    system(paste0('rm -r ', output_wd, '/by-id'), ignore.stderr = T)
    system(paste0('rm -r ', output_wd, '/met'), ignore.stderr = T)
    system(paste0('rm -r ', output_wd, '/particles'), ignore.stderr = T)
  }
  for (d in c('by-id', 'particles', 'footprints')) {
    d = file.path(output_wd, d)
    if (!file.exists(d))
      dir.create(d, recursive = T)
  }

  if (!run_trajec & !run_foot) stop('no need for parallel computing since run_trajec = F and run_foot = F ...\n')

  # removed 'isot', 'iconvert', 'outfrac', 'random' -- 
  # not used by the latest STILT
  output = xstilt_apply(FUN = simulation_step,
                        slurm = slurm, 
                        slurm_options = slurm_options,
                        n_cores = n_cores,
                        n_nodes = n_nodes,
                        jobname = jobname, 
                        before_footprint = list(before_footprint),
                        before_trajec = list(before_trajec),
                        simulation_id = simulation_id,
                        capemin = capemin,
                        cmass = cmass,
                        conage = conage,
                        cpack = cpack,
                        delt = delt,
                        dxf = dxf,
                        dyf = dyf,
                        dzf = dzf,
                        efile = efile,
                        emisshrs = emisshrs,
                        frhmax = frhmax,
                        frhs = frhs,
                        frme = frme,
                        frmr = frmr,
                        frts = frts,
                        frvs = frvs,
                        hnf_plume = hnf_plume,
                        horcoruverr = horcoruverr,
                        horcorzierr = horcorzierr,
                        hscale = hscale,
                        ichem = ichem,
                        idsp = idsp,
                        initd = initd,
                        k10m = k10m,
                        kagl = kagl,
                        kbls = kbls,
                        kblt = kblt,
                        kdef = kdef,
                        khinp = khinp,
                        khmax = khmax,
                        kmix0 = kmix0,
                        kmixd = kmixd,
                        kmsl = kmsl,
                        kpuff = kpuff,
                        krand = krand,
                        krnd = krnd,
                        kspl = kspl,
                        kwet = kwet,
                        kzmix = kzmix,
                        lib.loc = lib.loc,
                        maxdim = maxdim,
                        maxpar = maxpar,
                        met_file_format = met_file_format,
                        met_path = met_path,
                        met_subgrid_buffer = met_subgrid_buffer,
                        met_subgrid_enable = met_subgrid_enable,
                        met_subgrid_levels = met_subgrid_levels,
                        mgmin = mgmin,
                        n_hours = n_hours,
                        n_met_min = n_met_min,
                        ncycl = ncycl,
                        ndump = ndump,
                        ninit = ninit,
                        nstr = nstr,
                        nturb = nturb,
                        numpar = numpar,
                        nver = nver,
                        outdt = outdt,
                        output_wd = output_wd,
                        p10f = p10f,
                        pinbc = pinbc,
                        pinpf = pinpf,
                        poutf = poutf,
                        wion = projection,
                        qcycle = qcycle,
                        r_run_time = receptors$run_time,
                        r_lati = receptors$lati,
                        r_long = receptors$long,
                        r_zagl = receptors$zagl,
                        rhb = rhb,
                        rht = rht,
                        rm_dat = rm_dat,
                        run_foot = run_foot,
                        run_trajec = run_trajec,
                        siguverr = siguverr,
                        sigzierr = sigzierr,
                        smooth_factor = smooth_factor,
                        splitf = splitf,
                        stilt_wd = xstilt_wd,
                        time_integrate = time_integrate,
                        timeout = timeout,
                        tkerd = tkerd,
                        tkern = tkern,
                        tlfrac = tlfrac,
                        tluverr = tluverr,
                        tlzierr = tlzierr,
                        tout = tout,
                        tratio = tratio,
                        tvmix = tvmix,
                        varsiwant = list(varsiwant),
                        veght = veght,
                        vscale = vscale,
                        vscaleu = vscaleu,
                        vscales = vscales,
                        w_option = w_option,
                        wbbh = wbbh,
                        wbwf = wbwf,
                        wbwr = wbwr,
                        wvert = wvert,
                        xmn = xmn,
                        xmx = xmx,
                        xres = xres,
                        ymn = ymn,
                        ymx = ymx,
                        yres = yres,
                        zicontroltf = zicontroltf,
                        ziscale = ziscale,
                        z_top = z_top,
                        zcoruverr = zcoruverr, 

  # pass additional variables to stilt_apply and then to simulation_step() 
  # needed for before_*_xstilt() for X-STILT, DW, 07/01/2021
                        met = met, 
                        obs_sensor = obs_sensor, 
                        obs_species = list(obs_species),
                        obs_fn = obs_fn, 
                        ak_wgt = namelist$ak_wgt, 
                        pwf_wgt = namelist$pwf_wgt, 
                        xres2 = list(namelist$foot_res2), 
                        yres2 = list(namelist$foot_res2), 
                        time_integrate2 = list(time_integrate2), 
                        foot_nhrs = namelist$foot_nhrs,
                        run_hor_err = run_hor_err,
                        emiss_fn = emiss_fn, 
                        ct_ver = namelist$ct_ver, 
                        ctflux_path = namelist$ctflux_path, 
                        ctmole_path = namelist$ctmole_path)

}


# for debugging
if (F) {
  
  X = 1
  r_run_time = receptors$run_time[X]
  r_lati = receptors$lati[X]
  r_long = receptors$long[X]
  r_zagl = receptors$zagl[X]
  simulation_id = simulation_id[X]

  args = list(obs_sensor = obs_sensor, 
              obs_species = list(obs_species),
              obs_fn = obs_fn, 
              ak_wgt = namelist$ak_wgt, 
              pwf_wgt = namelist$pwf_wgt, 
              xres2 = list(namelist$foot_res2), 
              yres2 = list(namelist$foot_res2), 
              time_integrate2 = list(time_integrate2), 
              foot_nhrs = namelist$foot_nhrs,

              run_hor_err = run_hor_err,
              emiss_fn = NA, 
              ct_ver = namelist$ct_ver, 
              ctflux_path = namelist$ctflux_path, 
              ctmole_path = namelist$ctmole_path)

  rundir = file.path(output_wd, 'by-id', simulation_id)
  output = list()
  output$file = file.path(rundir, paste0(simulation_id, '_traj.rds'))
  output = readRDS(output$file)
  


}

