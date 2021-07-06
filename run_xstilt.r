#' create a namelist for running X-STILT trajectories
#' @author Dien Wu
#'
# ---------------------------------------------------------------------------
#' A. run_trajec == T or run_foot == T: run trajectory or footprint, simulations
#' A1. if run_hor_err == T or run_ver_err == T 
#'     X-STILT generates trajec with perturbations from wind or PBL errors.
#'         along with traj-level CO2 and its error statistics to info.rds files
#' A2. if run_emiss_err == T 
#'     X-STILT generates trajec and footprint with a resolution of 0.1deg, 
#'     since the absolute emission error has res of 0.1deg. 
#'
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
#' *** You can also assign your receptor time to @param timestr. HOWEVER it will 
#' *** overwrite the existing time column in the receptor file if there is. 
# ---------------------------------------------------------------------------
# 
# Latest major changes
# ---------------------------------------------------------------------------
# upgrade to use STILT-HYSPLIT and a brand new refactoring on tht X-STILT part 
# place receptor locations according to a txt file provided by the user 
# getting rid of the OCO dependencies (see 'run_xstilt_ideal.r'), DW, 12/06/2020

# ---------------------------------------------------------------------------
# add X-STILT simulations for TROPOMI soundings, will look for the AK profiles 
# of the nearest soundings to receptor locations, DW, 07/01/2021

# add a jittering option to sample additional receptors within 
# a satellite sounding, DW, 07/03/2021


# ----------- Dependencies and API key for geolocation (must-have) ---------- #
## go to X-STILT dir and source functions and libraries
homedir   = '/central/home/dienwu'
xstilt_wd = file.path(homedir, 'models/X-STILT') 
setwd(xstilt_wd); source('r/dependencies.r')

# *** there is a current bug with the rslurm package, you may need to 
# mannually pull the development version from github by doing:
#devtools::install_github('SESYNC-ci/rslurm')               # DW, 11/6/2020

# Please insert your API in the 'insert_ggAPI.csv' for use of ggplot and ggmap
# google API can be obtained from https://console.developers.google.com/
api.key = readLines('insert_ggAPI.csv')
if (api.key == '') stop('Missing googleAPI, insert your API in insert_ggAPI.csv\n')
register_google(key = api.key)


# ------------------------ City params (must-have) -------------------------- #
cat('Enter your city name: ')
site = readLines('stdin', n = 1); print(site)

# define entire- and inner urban- domain (affect receptor density if selTF = T)
dlon = 0.5           # e.g., dlon of 0.5 means 1 x 1 degree box around the site
dlat = 0.5       
urban_dlon = 0.3     # urban box defined as city.lat +/- dlat, city.lon +/- dlon
urban_dlat = 0.3     # dlat/dlon in degrees 

# automatically obtain the city lat/lon (requires Google API key)
# for mannually insert coordinates, set city.loc = data.frame(lon = , lat = )
lon_lat = get.lon.lat(site, dlon, dlat, city.loc = NULL)


# ------------------------ I/O params (must-have) --------------------------- #
# as a default, all input data are stored under input_path, change if you like
input_path  = '/central/groups/POW'
store_path  = file.path(input_path, 'XSTILT_output', site)

# *** modify the info path/directory that stores OCO-2/3 or TROPOMI data -------
obs_sensor  = c('OCO-2', 'OCO-3', 'TROPOMI', NA)[2]
obs_ver     = c('V10r', 'VEarlyR', NA)[2]       # retrieval algo ver if there is
obs_species = c('CO2', 'CO', 'NO2', 'CH4')[1]   # can only choose one for TROPOMI
oco_path = file.path(input_path, obs_sensor, paste0('L2_Lite_FP_', obs_ver))
sif_path = file.path(input_path, obs_sensor, paste0('L2_Lite_SIF_', obs_ver))
trp_path = file.path(input_path, obs_sensor, obs_species)

# *** if obs_sensor is NA, perform ideal simulation without satellite dependence, 
# *** instead, the user needs to provide a .csv/.txt file for receptor info
# It should contain desired longitude ('long') and latitude ('lati'), 
# one can also include a third column for receptor time ('time') OR enter the 
# receptor time to `timestr` below. Each column needs to be separated by comma. 

if (is.na(obs_sensor)) {                # X-simulations WOUT satellite data
  recp_fn = 'receptor_demo.csv'                               # <path/filename>
  obs_ver = obs_species = obs_path = sif_path = NA
} else {                                # X-simulations WITH satellite data
  recp_fn = NA
  obs_path = ifelse(obs_sensor == 'TROPOMI', trp_path, oco_path)
}

# paths for radiosonde, ODIAC, chemical transport model (optional) -------------
raob_path  = file.path(input_path, 'RAOB', site)  # NOAA radiosonde
odiac_ver  = c('2019', '2020')[2]         # ODIAC version
odiac_path = file.path(input_path, 'ODIAC', paste0('ODIAC', odiac_ver))  


# ------------------ obtaining overpass time string -------------------- #
# There are several options included here, if oobs_sensor is 
# OCO-2/3 - look for overpasses with sufficient data over area around site;
# TROPOMI - since it's daily, time string is provided by user
# NA - ideal run without satellite data, time string's provided by user
# see get.timestr() for more details, DW, 07/05/2021
timestr = get_timestr(site, lon_lat, obs_sensor, obs_ver, obs_path, store_path, 
                      recp_fn, plotTF = FALSE, urbanTF = TRUE, urban_dlon, 
                      urban_dlat, sif_path, searchTF = FALSE, 
                      date.range = c('20140101', '20211231'))

cat('Done with choosing cities & overpasses...\n')
# --------------------------------------------------------------------------- #


# --------------------- Basis X-STILT flags (must-have) --------------------- #
run_trajec    = T    # whether to generate trajec; T: may overwrite existing trajec
run_foot      = T    # whether to generate footprint
run_hor_err   = F    # T: set error parameters
run_ver_err   = F    # T: set error parameters
run_emiss_err = F    # T: get XCO2 error due to prior emiss err
run_wind_err  = F    # T: calc wind error based on RAOB 
run_sim       = F    # T: calc XFF or error with existing footprint (only for CO2)
nhrs          = -12  # number of hours backward (-) or forward (+)

# output variable names required in trajec.rds
varstrajec = c('time', 'indx', 'lati', 'long', 'zagl', 'zsfc', 'foot', 'samt',
               'dmas', 'mlht', 'temz', 'pres', 'sigw', 'tlgr', 'dens', 'sphu')


# ---------------------- Receptor params (must-have) ------------------------ #
# line source for column release according to HYSPLITv5, DW, 09/17/2020
# all particles are roughly evenly distributed between minagl and maxagl
minagl = 0                # min release height in meter AGL
maxagl = 3000             # max release height in meter AGL
numpar = 3000             # total number of particles between minagl and maxagl

# Receptor selection ---------------------------------------------------------
selTF = TRUE          # T: selected soundings as receptors; F: use all soundings

# T: create receptors within each satellite sounding besides the centered lat/lon
# F: place receptors ONLY at the centered lat/lon of a sounding, DW, 07/02/2021
jitterTF   = TRUE            # T: works better for TROPOMI with larger polygons
num_jitter = 5               # number of additional receptors per sounding

#' Only place receptors for soundings that qualify @param obs_filter
#' here are some choices for OCO-2/3 and TROPOMI (uncomment the one you need)
#obs_filter = c('QF', 0)              # select OCO soundings with QF = 0
#obs_filter = c('QA', 0.5)            # select TROPOMI soundings with QA >= 0.5 
obs_filter = NULL                     # use all soundings regardless of QF or QA

#' number of soundings desired within 1 degree lat range (background or enhanced)
#' peak range is defined by @param urban_dlat (see Line 75)
#num_bg = num_peak = NA                # if NA, no need to select soundings
#num_bg = 50 ; num_peak = 150         # num of soundings/deg selected for OCO-2
num_bg = 100; num_peak = 500         # num of soundings/deg selected for OCO-3

# *** For ideal simulations, no need to rely on satellite data, set to NA or FALSE 
# receptors will only based placed based on lati/long info from receptor_demo.csv
if (is.na(obs_sensor)) { num_bg = num_peak = num_jitter = NA; selTF = jitterTF = FALSE }


# ------------------- ARL format meteo params (must-have) -------------------- #
# see STILTv2 https://uataq.github.io/stilt/#/configuration?id=meteorological-data-input
met = c('gdas0p5', 'gfs0p25', 'hrrr')[2]         # choose met fields
met_path = file.path(homedir, met)               # path of met fields
n_met_min = 1                                    # min number of files needed
met_file_format = '%Y%m%d'                       # met file name convention

# OPTION for subseting met fields if met_subgrid_enable is on, 
# useful for large met fields like GFS or HRRR
met_subgrid_buffer = 0.1   # Percent to extend footprint area for met subdomain
met_subgrid_enable = T    

# if set, extracts the defined number of vertical levels from the meteorological 
# data files to further accelerate simulations, default is NA
met_subgrid_levels = NA    
cat('Done with params for receptor and meteo- setup...\n')


# -------------------- Footprint params (must-have) ------------------------- #
# whether weight footprint by AK and PW for column simulations 
# if F, AK = 1; if T, look for AK at the closest sounding from real data
ak_wgt  = TRUE          # *** if obs_sensor is NA, ak_wgt will be forced as FALSE         
pwf_wgt = TRUE
if (is.na(obs_sensor)) ak_wgt = FALSE 

# footprint spatial domain defined as city.lat +/- foot_dlat and 
# city.lon +/- foot_dlon in degrees, 10 here meaning 20 x 20deg box around city
foot_dlat = 10     
foot_dlon = 10 

# ---------------------------------------------------------------------------- #
# *** You can run multiple sets of footprints with different spatiotemporal 
#     resolution from the same trajec at one time (see below), DW, 07/01/2021
#' @param MAIN_RUN (e.g., 1km time-integrated footprint) -----------------------
foot_res  = 1/120       # spatial resolution in degree
foot_nhrs = nhrs        # if foot_nhrs < nhrs, subset trajec before footprint
time_integrate = TRUE   # F: hourly foot; T: time-integrated foot using all foot_nhrs

#' @param OPTIONAL_RUN with diff configurations (e.g., hourly 0.05 deg foot) ---
#' both params can be a vector, footprint filename contains res info
#' if no need for alternative runs, set @param foot_res to NA, DW, 02/11/2019
foot_res2 = c(NA, 1/10, 1/20, 1)[3]     # spatial res in degree
time_integrate2 = FALSE                 # ndim(time_integrate2) = ndim(foot_res2)

# ---------------------------------------------------------------------------- #
# other neccesary footprint params using STILTv2 (Fasoli et al., 2018)
hnf_plume     = TRUE                # T: hyper near-field (HNP) for mixing hgts
smooth_factor = 1                   # Gaussian smooth factor, 0 to disable
projection    = '+proj=longlat'
cat('Done with params for error analysis and footprint setup...\n')
# ---------------------------------------------------------------------------- #


# ---------- Transport error params (only needed for CO2 error) -------------- #
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
slurm = T                           # T: SLURM parallel computing
n_nodes = 3
n_cores = 7
timeout = 12 * 60 * 60              # time allowed before terminations in sec
job_time = '12:00:00'               # total job time
slurm_account = 'pow'
slurm_partition = 'any'

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
                hnf_plume = hnf_plume, jitterTF = jitterTF, job_time = job_time, 
                lon_lat = list(lon_lat), mem_per_node = mem_per_node, met = met, 
                met_file_format = met_file_format, met_path = met_path, 
                met_subgrid_buffer = met_subgrid_buffer, 
                met_subgrid_enable = met_subgrid_enable, 
                met_subgrid_levels = met_subgrid_levels, minagl = minagl, 
                maxagl = maxagl, nhrs = nhrs, n_cores = n_cores, 
                n_met_min = n_met_min, n_nodes = n_nodes, num_bg = num_bg, 
                num_jitter = num_jitter, num_peak = num_peak, numpar = numpar, 
                obs_filter = list(obs_filter), obs_path = obs_path, 
                obs_sensor = obs_sensor, obs_species = obs_species, 
                obs_ver = obs_ver, odiac_path = odiac_path, odiac_ver = odiac_ver, 
                projection = projection, pwf_wgt = pwf_wgt, 
                raob_path = raob_path, recp_fn = recp_fn, 
                run_emiss_err = run_emiss_err, run_foot = run_foot, 
                run_hor_err = run_hor_err, run_sim = run_sim, 
                run_trajec = run_trajec, run_ver_err = run_ver_err, 
                run_wind_err = run_wind_err, selTF = selTF, site = site, 
                slurm = slurm, slurm_account = slurm_account, 
                slurm_partition = slurm_partition, smooth_factor = smooth_factor, 
                store_path = store_path, time_integrate = time_integrate, 
                time_integrate2 = list(time_integrate2), timeout = timeout, 
                timestr = timestr, urban_dlat = urban_dlat, 
                varstrajec = varstrajec, xstilt_wd = xstilt_wd, zisf = zisf)      
cat('Done with creating namelist...Configuring\n\n')
config_xstilt(namelist)  # start X-STILT, either calc traj, foot or simulation
q('no')
# end of main script
# ---------------------------------------------------------------------------- #