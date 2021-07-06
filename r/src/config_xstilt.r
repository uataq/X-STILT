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
  timestr   = namelist$timestr
  lon_lat   = namelist$lon_lat[[1]]
  obs_ver   = namelist$obs_ver
  obs_path  = namelist$obs_path
  obs_sensor = namelist$obs_sensor
  obs_species = namelist$obs_species
  store_path = namelist$store_path

  # make sure obs_* are set to NA and no AK weighting for ideal runs
  if (is.na(obs_sensor)) { obs_ver = obs_path = obs_species = NA; namelist$ak_wgt = FALSE }

  # Model control -------------------------------------------------------------
  rm_dat      = T
  timeout     = namelist$timeout  # in sec
  run_trajec  = namelist$run_trajec
  run_foot    = namelist$run_foot
  run_sim     = namelist$run_sim
  run_hor_err = namelist$run_hor_err
  run_ver_err = namelist$run_ver_err
  run_wind_err = namelist$run_wind_err
  
  if (run_trajec) cat('Need to generate trajec...\n')
  if (run_foot)   cat('Need to generate footprint...\n')
  if (!run_trajec & !run_foot & !run_sim) 
    cat('NO calculations will be performed, please check run_* flags..\n')
  if (run_sim) { run_trajec = F; run_foot = F }   
  if (run_trajec | run_foot) run_sim = F   # when running trajec/foot, no sim allowed  

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
  xmn  = round(lon_lat$citylon) - namelist$foot_dlon
  xmx  = round(lon_lat$citylon) + namelist$foot_dlon
  ymn  = round(lon_lat$citylat) - namelist$foot_dlat
  ymx  = round(lon_lat$citylat) + namelist$foot_dlat
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
  outlist = create.outwd(timestr, obs_species, obs_sensor, obs_path, lon_lat, 
                         store_path, met, run_hor_err)
  obs_fn  = outlist$obs_info$fn    # get observation file name if satellite is used
  output_wd = outlist$output_wd

  # Compute XCO2 using existing traj & foot ------------------------------------
  # calculate XCO2 concentration and its error (need trajec and footprint ready)
  if ( !run_trajec & !run_foot & run_sim & obs_species == 'CO2') {

    # ------------------------  Horizontal trans error ---------------------- #
    ### simulate transport error in XCO2 due to met errors, DW, 07/25/2018
    # requires two sets of trajectories before running the following section:
    if (run_hor_err) { # this does not need footprint

      ## call function cal.trans.err() to estimate receptor-level trans err [ppm]
      # get actual ppm error, need to have error statistics ready
      # see cal.trajfoot.stat() in called before_footprint_xstilt.r for err stat
      cat('Start simulations of XCO2 error due to horizontal trans err...\n')
      result = cal.trans.err(site, timestr, workdir = xstilt_wd, 
                             outdir = output_wd, store_path, met)
      if (is.null(result)) stop('No results calculated, check cal.trans.err()\n')

    } else {

      # ------------------- XCO2 or emiss error simulations ----------------- #
      # requires trajec and footprints ready for running things below, DW, 06/04/2018  
      # call func to match ODIAC emissions with xfoot & sum up to get 'dxco2.ff'
      cat('Start simulations of XCO2.ff or its error due to emiss err...\n')
      result = run.xco2ff.sim(site, timestr, vname = namelist$odiac_ver, 
                              tiff.path = namelist$odiac_path, 
                              output_wd, foot.res = xres, xstilt_wd, store_path,
                              nhrs = n_hours, obs_sensor, obs_ver, met, 
                              run_emiss_err = namelist$run_emiss_err, 
                              edgar.file = namelist$edgar_file, 
                              ffdas.file = namelist$ffdas_file)
      if (is.null(result)) stop('No results calculated, check run.xco2ff.sim()\n')
      return(result)
    } # end if run_hor_err

  } else if (run_sim & obs_species != 'CO2') 
    stop('Currently, the model only calculates the FF signals of XCO2...\n\n')
  # end if run_sim


  # Receptor locations & time -------------------------------------------------
  # IF for ideal simulation, read from receptor file 
  if (is.na(obs_sensor)) {
    receptors = read.table(namelist$recp_fn, header = T, sep = ',')
    
    if (!is.na(namelist$timestr)) { # if no time column found, use @param timestr
      receptors$run_times = timestr 
    } else receptors = receptors %>% rename(run_times = time)
    
    nchart = nchar(receptors$run_times[1])
    if (nchart == 8 ) formatt = '%Y%m%d'
    if (nchart == 10) formatt = '%Y%m%d%H'
    if (nchart == 12) formatt = '%Y%m%d%H%M'
    if (nchart == 14) formatt = '%Y%m%d%H%M%S'
    if (!nchart %in% c(8, 10, 12, 14)) 
      stop('Incorrect form of time string...please check @param recp_fn or @param timestr\n')

    receptors$run_time = as.POSIXct(receptors$run_times, 'UTC', format = formatt) 
    receptors$zagl = list(agl)

  } else {    # IF for simulations using satellite data

    peak_lat = c(lon_lat$citylat - namelist$urban_dlat, 
                 lon_lat$citylat + namelist$urban_dlat)
    receptors = get.recp.sensor(timestr, obs_filter = unlist(namelist$obs_filter), 
                                obs_fn, obs_sensor, obs_ver, obs_path, lon_lat, 
                                selTF = namelist$selTF, jitterTF = namelist$jitterTF, 
                                num_jitter = namelist$num_jitter, peak_lat, 
                                num_bg = namelist$num_bg, 
                                num_peak = namelist$num_peak, 
                                agl, run_trajec, output_wd)
  } # end if
  cat(paste('Done with receptor setup...total', nrow(receptors), 'receptor(s)..\n\n'))


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
  kmix0       = 250
  kmixd       = 3
  kmsl        = 0
  kpuff       = 0
  krand       = 4
  krnd        = 6
  kspl        = 1
  kwet        = 1
  kzmix       = 0
  maxdim      = 1
  maxpar      = namelist$numpar
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
  z_top       = 25000; if (toupper(met) == 'NARR') z_top = 15000 

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
                          vscale = vscale, vscaleu = vscaleu, vscales = vscales, 
                          wbbh = wbbh, wbwf = wbwf, wbwr = wbwr, winderrtf = 0, 
                          wvert = wvert, zicontroltf = zicontroltf)
  
  # Get error stats if needed -------------------------------------------------
  # Calculating error stats for horizontal and vertical transport errors
  errlist  = config_trans_err(namelist, site, lon_lat, timestr, xstilt_wd, simstep_namelist)
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
  #mem_per_cpu = 10 * 1024    
  slurm_options = list(time = namelist$job_time, 
                       account = namelist$slurm_account, 
                       partition = namelist$slurm_partition)
  if (!is.null(namelist$mem_per_node)) slurm_options$mem = namelist$mem_per_node
  jobname = paste0('XSTILT_', site, '_', timestr, '_', obs_sensor, '_', obs_species)
  if (is.na(obs_sensor)) jobname = paste0('XSTILT_', site, '_', timestr, '_ideal')


  # Startup messages -----------------------------------------------------------
  message('\n\nInitializing X-STILT')
  message('Number of receptors: ', nrow(receptors))
  message('Number of parallel threads: ', n_nodes * n_cores)


  # Source dependencies --------------------------------------------------------
  setwd(xstilt_wd)
  source('r/dependencies.r')


  # Structure out directory ------------------------------------------------------
  # Outputs are organized in three formats. by-id contains simulation files by
  # unique simulation identifier. particles and footprints contain symbolic links
  # to the particle trajectory and footprint files in by-id
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
                        lib.loc = lib.loc,
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
                        obs_species = obs_species,
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

  args = list(obs_sensor = obs_sensor, 
              obs_species = obs_species,
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

