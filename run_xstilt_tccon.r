#' create a namelist for running X-STILT trajectories
#' @author Dien Wu
#'
# ---------------------------------------------------------------------------
#' A. run_trajec == T or run_foot == T: run trajectory or footprint, simulations
#' A1. if run_hor_err == T or run_ver_err == T 
#'     X-STILT generates trajec with perturbations from wind or PBL errors.
#'     along with traj-level CO2 and its error statistics (see AA1 for stats)
# -----------------------------------------------------------------------------
#' AA1. if run_hor_err == T & run_wind_err == T
#'      X-STILT calculates wind uncertainties and biases using NOAA radiosonde
#'      data, one needs to download FSL data from https://ruc.noaa.gov/raobs/
#'      and put data under @param raob_path as @param raob_fn
# -----------------------------------------------------------------------------
#' A2. if run_emiss_err == T 
#'     X-STILT generates trajec and footprint with a resolution of 0.1deg, 
#'     since the absolute emission error has res of 0.1deg. 
#' -----------------------------------------------------------------------------
#' B. run_trajec == F & run_foot == F & run_sim = T: requires trajec and foot
#'    X-STILT models XCO2 enhancements based on FFCO2 emission from 1km ODIAC;
#'                or XCO2 errors due to emission errors (run_emiss_err == T);
#'                or XCO2 errors due to transport errors (run_hor_err == T).
#' ---------------------------------------------------------------------------
#' Instructions for ideal column simulations 
#' *** AK is simply treated as 1 (@param ak_wgt == FALSE)
#' *** PW is calculated using modeled met variables (@param pwf_wgt == TRUE)
#' *** The user NEEDS to provide a txt or csv file for receptor locations
#' *** correct file should include long, lat, and time (optional, YYYYMMDDHH). 
#' *** You can also assign your receptor time to @param timestr. HOWEVER it'll 
#' *** overwrite the existing time column in the receptor file if there is. 
# ---------------------------------------------------------------------------
# 
# Latest major changes, DW
# 1. upgrade to STILT-HYSPLITv5
# 2. allow for ideal sim w/o obs using a csv file provided by users, 12/06/2020
# 3. allow for working with TROPOMI data, 07/03/2021
# 4. allow for sounding selection (OCO and TROPOMI) within NEAR- or FAR-FIELD, 
#    see comments around L53-60, L154-160), 06/17/2022
# ---------------------------------------------------------------------------

# !!!! to run this code in R session: source('run_xstilt.r')
options(timeout = max(1500, getOption('timeout')))   # for download.file

# ----------- Dependencies and API key for geolocation (must-have) ---------- #
## go to X-STILT dir and source functions and libraries
homedir = '/central/home/dienwu'
xstilt_wd = file.path(homedir, 'X-STILT') 
setwd(xstilt_wd); source('r/dependencies.r')

# --------------------- Site location params (must-have) --------------------- #
site = 'Lamont'

# define entire- and near-field domain (usually for high observations)
# for selecting observations, if num_* is not NA
# ENTIRE domain - site_lat +/- dlat, site_lon +/- dlon 
# NEAR-FIELD domain - site_lat +/- nf_dlat, site_lon +/- nf_dlon
dlon = 1            # e.g., dlon = 0.5 means 1 x 1 degree box around the site
dlat = 1            # dlat/dlon in degrees 
nf_dlon = 0.3
nf_dlat = 0.3

# - automatically obtain lat/lon (requires Google API key in insert_ggAPI.csv)
# - OR mannually insert coordinates, e.g., site.loc = data.frame(lon = , lat = )
site_loc = data.frame(lon = -97.486, lat = 36.604)
lon_lat = get.lon.lat(site, dlon, dlat, site.loc = site_loc)


# ------------------------ I/O params (must-have) --------------------------- #
# as a default, all input data are stored under input_path, change if you like
input_path  = '/central/groups/POW'
store_path  = file.path(input_path, 'XSTILT_output', site)

# *** modify the info path/directory that stores column data -------
obs_sensor  = c('OCO-2', 'OCO-3', 'TROPOMI', 'TCCON', 'EM27', NA)[4]

# retrieval ver, e.g., V11r for OCO-2, V10p4r for OCO-3, GGG2020 for TCCON
# if not applicable, set it to NA ---------
obs_ver = c('V11r', 'V10p4r', 'GGG2020', NA)[3]     

# currently only allow one species per time for OCO and TROPOMI
# fine with multiple species for TCCON or EM27/SUN --------
obs_species = c('CO2', 'CO', 'NO2', 'CH4', 'N2O')[1]   

# obs paths and files (including OCO-2/3, OCO SIF, TROPOMI, TCCON)
oco_path = file.path(input_path, obs_sensor, paste0('L2_Lite_FP_', obs_ver))
sif_path = file.path(input_path, obs_sensor, paste0('L2_Lite_SIF_', obs_ver))
trp_path = file.path(input_path, obs_sensor, obs_species)
tcn_path = file.path(input_path, obs_sensor, obs_ver)

# *** one need to specific the obs file if for TCCON and EM27 ------------
#obs_fn = NA   # obs_fn = NA is the default
obs_fn = file.path(tcn_path, 'oc20110416_20221022.public.qc.nc')

# -----------------------------------------------------------------------------
# *** if obs_sensor is NA, perform ideal runs without satellite dependence, 
# *** The user needs to provide a .csv/.txt file for receptor info
# Such file should contain desired longitude ('long') and latitude ('lati'), 
# A third column for receptor time ('time') is optional OR enter the 
# receptor time as `timestr` below. Each column needs to be separated by comma. 
# -----------------------------------------------------------------------------
if ( is.na(obs_sensor) ) {   # X-simulations WOUT satellite data
  recp_fn = 'receptor_demo.csv'                            # <path/filename>
  obs_ver = obs_species = obs_path = sif_path = NA
} else {                   # X-simulations WITH satellite data
  recp_fn = NA
  obs_path = ifelse(obs_sensor == 'TROPOMI', trp_path, 
             ifelse(obs_sensor == 'TCCON', tcn_path, oco_path))
}

# (optional) paths for radiosonde (for transport error analysis) and ODIAC ----
raob_path = file.path(input_path, 'XSTILT_dep/RAOB')          # NOAA radiosonde
raob_fn   = list.files(raob_path, '.tmp', full.names = T)
odiac_ver = c('2019', '2020')[2]                             # ODIAC version
odiac_path = file.path(input_path, 'ODIAC', paste0('ODIAC', odiac_ver))  


# ---------------------- obtaining overpass time string --------------------- #
#' @param obs_sensor == TROPOMI - daily obs, time string provided by user
#' @param obs_sensor == NA - ideal run w/o obs, time string provided by user
#' @param timestr in format of YYYYMMDD or YYYYMMDDHH (either works)
timestr = '20200701'

#' @param obs_sensor == OCO-2/3 - can help search for overpasses with #
#' sufficient data over the entire area and/or near-field area around site
if ( grepl('OCO', obs_sensor) ) 
  timestr = get_timestr(site, lon_lat, obs_sensor, obs_ver, obs_path, 
                        store_path, recp_fn, plotTF = F, qfTF = T, nfTF = T, 
                        nf_dlon, nf_dlat, sif_path, searchTF = F, 
                        date_range = c('20140101','20231231'))

#' @param obs_sensor == TCCON - distribute obs times per thred
#' insert a valid time range, e.g., work on TCCON from 20 UTC to 21 UTC in 2020
if ( grepl('TCCON', obs_sensor) ) {
  sep = c('week', 'day', 'hour', 'min')[4]
  dates = seq(as.POSIXct('2020-07-01 20:00:00', tz = 'UTC'), 
              as.POSIXct('2020-07-01 20:59:59', tz = 'UTC'), by = sep)
  timestr = unique(format(dates, tz = 'UTC', format = '%Y%m%d%H'))
}

cat('Done with choosing cities & overpasses...\n')
# --------------------------------------------------------------------------- #


# --------------------- Basis X-STILT flags (must-have) --------------------- #
run_trajec    = T    # if to generate trajec; T: may overwrite existing trajec
run_slant     = T    # recalc recp loc based on SZA, VZA, AZA by obs geometry
run_foot      = F    # if to generate footprint
run_hor_err   = F    # T: set error parameters
run_ver_err   = F    # T: set error parameters
run_emiss_err = F    # T: get XCO2 error due to prior emiss err
run_wind_err  = F    # T: calc wind error based on RAOB 
run_sim       = F    # T: calc XFF or error with existing foot (only for CO2)
nhrs          = -1  # number of hours backward (-) or forward (+)

# output variable names required in trajec.rds
varstrajec = c('time', 'indx', 'lati', 'long', 'zagl', 'zsfc', 'foot', 'samt',
               'dmas', 'mlht', 'temz', 'pres', 'sigw', 'tlgr', 'dens', 'sphu')


# ---------------------- Receptor params (must-have) ------------------------ #
# line source for column release according to HYSPLITv5, DW, 09/17/2020
# all particles are roughly evenly distributed between minagl and maxagl
minagl = 0                # min release height in meter AGL
maxagl = 3000             # max release height in meter AGL
numpar = 3000             # total number of particles between minagl and maxagl

# Receptor selection (params are set as NA for ideal sim) ---------------------
# T: create receptors within each satellite sounding besides centered lats/lons
# F: place receptors ONLY at the centered lat/lon of a sounding, DW, 07/02/2021
jitterTF   = FALSE            # T: works better for TROPOMI with larger polygons
num_jitter = 5                # number of additional receptors per sounding

#' Only place receptors for soundings that qualify @param obs_filter
#' here are some choices for OCO-2/3 and TROPOMI (uncomment the one you need)
#obs_filter = c('QF', 0)             # select OCO soundings with QF = 0
#obs_filter = c('QA', 0.5)            # select TROPOMI soundings with QA >= 0.5 
obs_filter = NULL                   # use all soundings regardless of QF or QA

#' evenly select soundings within near-field or far-field (background) regions
#' near-field: defined by @param nf_dlat & @param nf_dlon 
#' far-field: defined by @param dlat & @param dlon except for near-field region
# num of soundings/receptors along lat or lon within farfield or nearfield
# if ALL are NA - NO need to select soundings (see github README for more clues)
num_bg_lat = num_bg_lon = num_nf_lat = num_nf_lon = NA                
#num_bg_lat = 5; num_bg_lon = 5; num_nf_lat = 10; num_nf_lon = 10


# ------------------- ARL format meteo params (must-have) -------------------- #
# see STILTv2 https://uataq.github.io/stilt/#/configuration?id=meteorological-data-input
met = c('gfs0p25', 'hrrr', 'wrf27km')[2]         # choose met fields
met_path = file.path(homedir, met)     
met_file_format = '%Y%m%d'                       # met file name convention

#' if using your own metfields with no need to download files from ARL, 
#' set @param n_met_min to the min # of files needed and @param selfTF to TRUE
#' otherwise set @param n_met_min as NA and @param selfTF to FALSE
n_met_min = download.met.arl(timestr, met_file_format, nhrs, met_path, met, 
                             run_trajec, n_met_min = NA, selfTF = FALSE)

# OPTION for subseting met fields if met_subgrid_enable is on, 
# useful for large met fields like GFS or HRRR
met_subgrid_buffer = 0.1   # Percent to extend footprint area for met subdomain
met_subgrid_enable = F   

# if set, extracts the defined number of vertical levels from the meteorological 
# data files to further accelerate simulations, default is NA
met_subgrid_levels = NA    
cat('Done with params for receptor and meteo- setup...\n')


# -------------------- Footprint params (must-have) ------------------------- #
# whether weight footprint by AK and PW for column simulations 
# if F, AK = 1; if T, look for AK at the closest sounding from real data
ak_wgt  = TRUE         # *** if obs_sensor is NA, ak_wgt is forced as FALSE 
pwf_wgt = TRUE
if (is.na(obs_sensor)) ak_wgt = FALSE 

# footprint spatial domain defined as site.lat +/- foot_dlat and 
# site.lon +/- foot_dlon in degrees, 10 here meaning 20 x 20deg box around site
foot_dlat = 10     
foot_dlon = 10 

# ---------------------------------------------------------------------------- #
# *** You can run multiple sets of footprints with different spatiotemporal 
#     resolution from the same trajec at one time (see below), DW, 07/01/2021
#' @param MAIN_RUN (e.g., 1km time-integrated footprint) -----------------------
foot_res  = 1/120       # spatial resolution in degree
foot_nhrs = nhrs        # if foot_nhrs < nhrs, subset trajec before footprint
time_integrate = TRUE   # F: hourly foot; T: time-integrated foot

#' @param OPTIONAL_RUN with diff configurations (e.g., hourly 0.05 deg foot) ---
#' both params can be a vector, footprint filename contains res info
#' if no need for alternative runs, set @param foot_res to NA, DW, 02/11/2019
foot_res2 = c(NA, 1/10, 1/20, 1)[1]    # spatial res in degree
time_integrate2 = FALSE                # ndim(time_integrate2) = ndim(foot_res2)

# ---------------------------------------------------------------------------- #
# other neccesary footprint params using STILTv2 (Fasoli et al., 2018)
hnf_plume     = TRUE                # T: hyper near-field (HNP) for mixing hgts
smooth_factor = 1                   # Gaussian smooth factor, 0 to disable
projection    = '+proj=longlat'
cat('Done with params for error analysis and footprint setup...\n')
# ---------------------------------------------------------------------------- #


# ------- Transport error params (only needed for CO2 error, OPTIONAL) ------- #
# if run_hor_err = T, require ODIAC and CT fluxes and mole fractions
# to calculate horizontal transport error of total CO2, DW, 07/28/2018
if (run_hor_err) {
  ct_ver      = ifelse(substr(timestr, 1, 4) >= '2016', 'v2017', 'v2016-1')
  ct_path     = file.path(input_path, 'CT-NRT', ct_ver)
  ctflux_path = file.path(ct_path, 'fluxes/optimized')
  ctmole_path = file.path(ct_path, 'molefractions/co2_total')
} else ct_ver = ctflux_path = ctmole_path = NA       # end if run_hor_err

# if run_ver_err = T, calculate vertical transport error
# if run_ver_err = F, set zisf = 1 (default)
zisf = c(0.6, 0.8, 1.0, 1.2, 1.4)[3]; if (!run_ver_err) zisf = 1.0

# if run_emiss_err = T, calculating XCO2 error due to emiss error, 
# need EDGAR and FFDAS files to compute emission spread, DW, 10/21/2018
if (run_emiss_err) { 
  edgar_file = file.path(input_path, 'EDGAR/v42_CO2_2008_TOT.0.1x0.1.nc')
  ffdas_path = file.path(input_path, 'FFDAS')
  ffdas_file = list.files(ffdas_path, 'totals')
  ffdas_file = file.path(ffdas_path, ffdas_file[grep('2008', ffdas_file)])
} else edgar_file = ffdas_file = NA                   # end if run_emiss_err


# ----------------------------- SLURM params -------------------------------- #
# avoid using too many e.g., > 10 cores per node -> memory limits
# set slurm to FALSE, if system does not support slurm
# *** there is a current bug with the rslurm package, you may need to 
# mannually pull the development version from github by doing:
#devtools::install_github('SESYNC-ci/rslurm')               # DW, 11/6/2020
slurm = T                           # T: SLURM parallel computing
n_nodes = 12
n_cores = 10
if (!slurm) n_nodes = n_cores = 1
timeout = 12 * 60 * 60              # time allowed before terminations in sec
job_time = '12:00:00'               # total job time
slurm_account = 'pow'
slurm_partition = 'expansion'

# *** IF YOU EVER RAN INTO OOM-KILL ERROR, YOU CAN ENLARGE THE MAX MEM HERE *** 
# The ammount of memory per node you need in MB, extending to 10 GB per core
# given 180+ GB max memory on each Caltech cluster with 32 cores
# **** PLEASE adjust mem_per_* based on your machine
mem_per_core = 10                                 # max memory per core in GB
mem_per_node = n_cores * mem_per_core * 1024      # max mem per node now in MB 


# namelist required for generating trajec
namelist = list(ak_wgt = ak_wgt, ct_ver = ct_ver, ctflux_path = ctflux_path, 
                ctmole_path = ctmole_path, edgar_file = edgar_file, 
                ffdas_file = ffdas_file, foot_res = foot_res, 
                foot_res2 = list(foot_res2), foot_nhrs = foot_nhrs, 
                foot_dlat = foot_dlat, foot_dlon = foot_dlon, 
                hnf_plume = hnf_plume, jitterTF = jitterTF, 
                job_time = job_time, lon_lat = list(lon_lat), 
                mem_per_node = mem_per_node, met = met,                 
                met_file_format = met_file_format, met_path = met_path, 
                met_subgrid_buffer = met_subgrid_buffer, 
                met_subgrid_enable = met_subgrid_enable, 
                met_subgrid_levels = met_subgrid_levels, minagl = minagl, 
                maxagl = maxagl, nhrs = nhrs, n_cores = n_cores, 
                n_met_min = n_met_min, n_nodes = n_nodes, nf_dlat = nf_dlat, 
                nf_dlon = nf_dlon, num_jitter = num_jitter, 
                num_bg_lat = num_bg_lat, num_bg_lon = num_bg_lon, 
                num_nf_lat = num_nf_lat, num_nf_lon = num_nf_lon, 
                numpar = numpar, obs_filter = list(obs_filter), obs_fn = obs_fn,
                obs_path = obs_path, obs_sensor = obs_sensor, 
                obs_species = obs_species, odiac_path = odiac_path, 
                odiac_ver = odiac_ver, projection = projection, 
                pwf_wgt = pwf_wgt, raob_fn = raob_fn, recp_fn = recp_fn, 
                run_emiss_err = run_emiss_err, run_foot = run_foot, 
                run_hor_err = run_hor_err, run_sim = run_sim, 
                run_slant = run_slant, run_trajec = run_trajec, 
                run_ver_err = run_ver_err, run_wind_err = run_wind_err, 
                site = site, slurm = slurm, slurm_account = slurm_account,
                slurm_partition = slurm_partition,
                smooth_factor = smooth_factor, store_path = store_path, 
                time_integrate = time_integrate, 
                time_integrate2 = list(time_integrate2), timeout = timeout, 
                timestr = timestr, varstrajec = varstrajec, 
                xstilt_wd = xstilt_wd, zisf = zisf)
cat('Done with creating namelist...\n')
stop()
config_xstilt(namelist)  # start X-STILT, either calc traj, foot or simulation
q('no')
# end of main script
# ---------------------------------------------------------------------------- #