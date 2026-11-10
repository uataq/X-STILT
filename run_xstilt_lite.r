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
# 5. allow for working with TEMPO NO2 L2/L3 data, TROPOMI TCWV, and AWS downloads for meteo, a lite version without running error analysis, 11/09/2026,
# ---------------------------------------------------------------------------

# !!!! to run this code in R session: source('run_xstilt.r')
options(timeout = max(1500, getOption('timeout')))   # for download.file

# ----------- Dependencies and API key for geolocation (must-have) ---------- #
## go to X-STILT dir and source functions and libraries
homedir = '/home/C830382435/model'
xstilt_wd = file.path(homedir, 'X-STILT') 
source(file.path(xstilt_wd, 'r/dependencies.r'))

# Please insert your API in the 'insert_ggAPI.csv' for use of ggplot and ggmap
# google API can be obtained from https://console.developers.google.com/
api.key = readLines('insert_ggAPI_cp.csv')

# --------------------- Site location params (must-have) --------------------- #
#site = 'Bakerfield'
site = 'TwinFalls'

# define entire- and near-field domain (usually for high observations)
# for selecting observations, if num_* is not NA
# ENTIRE domain - site_lat +/- dlat, site_lon +/- dlon 
# NEAR-FIELD domain - site_lat +/- nf_dlat, site_lon +/- nf_dlon
dlon = 0.5            # e.g., dlon = 0.5 means 1 x 1 degree box around the site
dlat = 0.5            # dlat/dlon in degrees 
nf_dlon = 0.3
nf_dlat = 0.3

# - automatically obtain lat/lon (requires Google API key in insert_ggAPI.csv)
# - OR mannually insert coordinates, e.g., site.loc = data.frame(lon = , lat = )
lon_lat = get.lon.lat(site, dlon, dlat, site.loc = NULL, api.key = api.key)


# ------------------------ I/O params (must-have) --------------------------- #
# as a default, all input data are stored under input_path, change if you like
input_path = '/home/C830382435/data'
store_path = file.path('/home/C830382435/XSTILT_output', site)

# *** modify path/directory that stores OCO-2/3/TROPOMI/TEMPO data -------
# !!! cannot work with TEMPO L3, since it does not have AK in the file
obs_sensor = c('OCO-2', 'OCO-3', 'TROPOMI', 'TEMPO', NA)[3]
obs_ver    = c('V11r', 'V10p4r', 'V04', NA)[4]   # retrieval algo ver
obs_level  = c('L2_Lite_FP', 'L2')[2]                # data level

# only 1 gas per run
# arrange your input data as <sensor>/<level>/<species>, DW, 11/05/2026
obs_species = c('CO2', 'CO', 'NO2', 'CH4', 'HCHO', 'TCWV')[6]
obs_path    = file.path(input_path, obs_sensor, obs_level, obs_species)
store_path  = file.path(store_path, obs_sensor) 
print(store_path)


# *** if obs_sensor is NA, perform ideal runs without satellite dependence, 
# *** instead, the user needs to provide a .csv/.txt file for receptor info
# It should contain desired longitude ('long') and latitude ('lati'), 
# one can also include a third column for receptor time ('time') OR enter the 
# receptor time to `timestr` below. Each column needs to be separated by comma. 
if (is.na(obs_sensor)) {   # X-simulations WOUT satellite data
  recp_fn = 'receptor_demo.csv'                               # <path/filename>
  obs_ver = obs_path = sif_path = NA
} else recp_fn = NA
obs_fn = NA 


# ---------------------- obtaining overpass time string --------------------- #
#' @param obs_sensor == TROPOMI - daily obs, time string provided by user
#' @param obs_sensor == NA - ideal run w/o obs, time string provided by user
#' @param timestr in format of YYYYMMDD or YYYYMMDDHH (either works)
#' @param obs_sensor == TEMPO - hourly obs, need to specific hours; if hours are missing, it will loop over all hourly data in the day, DW, 11/05/2026 

#timestr = '2024072521'   # central valley for NO2
timestr = '2026060219'    # idaho for TCWV
if (is.na(obs_sensor)) timestr = unique(read.csv(recp_fn)$time)

#' @param obs_sensor == OCO-2/3 - can help search for overpasses with #
#' sufficient data over the entire area and/or near-field area around site
if (grepl('OCO', obs_sensor)) 
  timestr = get_timestr(site, lon_lat, obs_sensor, obs_ver, obs_path, 
                        store_path, recp_fn, plotTF = FALSE, qfTF = TRUE, 
                        nfTF = TRUE, nf_dlon, nf_dlat, sif_path, searchTF = F,
                        date_range = c('20140101', '20221231'))
cat('Done with choosing cities & overpasses...\n')
# --------------------------------------------------------------------------- #


# --------------------- Basis X-STILT flags (must-have) --------------------- #
run_trajec    = F    # if to generate trajec; T: may overwrite existing trajec
run_slant     = F    # recalc recp loc based on SZA, VZA, AZA by obs geometry
run_foot      = T    # if to generate footprint
store_totx    = F    # T: store particles info for columns above release hgt
nhrs          = -12  # number of hours backward (-) or forward (+)

# output variable names required in trajec.rds
varstrajec = c('time', 'indx', 'lati', 'long', 'zagl', 'zsfc', 'foot', 'samt',
               'dmas', 'mlht', 'temz', 'pres', 'sigw', 'tlgr', 'dens', 'sphu')


# ---------------------- Receptor params (must-have) ------------------------ #
# line source for column release according to HYSPLITv5, DW, 09/17/2020
# all particles are roughly evenly distributed between minagl and maxagl
minagl = 0                # min release height in meter AGL
maxagl = 2000             # max release height in meter AGL
numpar = 2000             # total number of particles between minagl and maxagl

# Receptor selection (params are set as NA for ideal sim) ---------------------
# T: create receptors within each satellite sounding besides centered lats/lons
# F: place receptors ONLY at the centered lat/lon of a sounding, DW, 07/02/2021
jitterTF   = FALSE            # T: works better for TROPOMI with larger polygons
num_jitter = 5                # number of additional receptors per sounding

#' Only place receptors for soundings that qualify @param obs_filter
#' here are some choices for OCO-2/3 and TROPOMI (uncomment the one you need)
#obs_filter = c('QF', 0)            # select OCO soundings with QF = 0
#obs_filter = c('QA', 0.7)          # select TROPOMI soundings with QA >= 0.7 
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
met             = c('gfs0p25', 'HRRR', 'wrf27km')[2]  # choose met fields
met_file_tres   = c('1 day', '6 hours', '1 hour')[2]  # met temporal resolution
met_path        = file.path(homedir, met)     
met_file_format = '%Y%m%d'                            # met file name convention
download_method = c('ftp', 'aws')[2]                  # whether ftp or aws

#' if you prefer to use your own metfields without downloading files from ARL, 
#' set @param n_met_min to the min # of files needed and @param selfTF to TRUE
#' otherwise set @param n_met_min as NA and @param selfTF to FALSE
n_met_min = download.met.arl(timestr, met_file_format, nhrs, met_path, 
                             tolower(met), met_file_tres, run_trajec, 
                             n_met_min = NA, selfTF = FALSE, download_method)


#' Note that you need to have AWS CLI installed and configured on your system.
#' AWS CLI download approach in terminal, an example is provided below:
#' aws s3 cp s3://noaa-oar-arl-hysplit-pds/hrrr/2024/07/ ./ --recursive --exclude "*" --include "*20240723*" --no-sign-request

# OPTION for subseting met fields if met_subgrid_enable is on, 
# useful for large met fields like GFS or HRRR
met_subgrid_buffer = 0.1   # Percent to extend footprint area for met subdomain
met_subgrid_enable = F   

# if set, extracts the defined number of vertical levels from the
# meteorological data files to further accelerate simulations, default is NA
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
foot_res  = 1/10       # spatial resolution in degree
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


# ----------------------------- SLURM params -------------------------------- #
# avoid using too many e.g., > 10 cores per node -> memory limits
# set slurm to FALSE, if system does not support slurm
# *** there is a current bug with the rslurm package, you may need to 
# mannually pull the development version from github by doing:
#devtools::install_github('SESYNC-ci/rslurm')               # DW, 11/6/2020
slurm   = T                           # T: SLURM parallel computing
n_nodes = 6
n_cores = 10
if (!slurm) n_nodes = n_cores = 1
timeout   = 12 * 60 * 60              # time allowed before terminations in sec
job_time  = '12:00:00'                # total job time
slurm_account = 'engineering'
#slurm_partition = 'wu_all'
slurm_partition = 'coe_all'

# *** IF YOU EVER RAN INTO OOM-KILL ERROR, YOU CAN ENLARGE THE MAX MEM HERE *** 
# The ammount of memory per node you need in MB, extending to 10 GB per core
# given 180+ GB max memory on each Caltech cluster with 32 cores
# **** PLEASE adjust mem_per_* based on your machine
mem_per_core = 10                                 # max memory per core in GB
mem_per_node = n_cores * mem_per_core * 1024      # max mem per node now in MB 

# namelist required for generating trajec
namelist = list(ak_wgt = ak_wgt, ct_ver = NA, ctflux_path = NA, 
                ctmole_path = NA, edgar_file = NA, ffdas_file = NA, 
                foot_res = foot_res, foot_res2 = list(foot_res2), 
                foot_nhrs = foot_nhrs, foot_dlat = foot_dlat, 
                foot_dlon = foot_dlon, hnf_plume = hnf_plume, 
                jitterTF = jitterTF, job_time = job_time,
                lon_lat = list(lon_lat), mem_per_node = mem_per_node, met = met,
                met_file_format = met_file_format,met_file_tres = met_file_tres,
                met_path = met_path, met_subgrid_buffer = met_subgrid_buffer, 
                met_subgrid_enable = met_subgrid_enable, 
                met_subgrid_levels = met_subgrid_levels, minagl = minagl, 
                maxagl = maxagl, nhrs = nhrs, n_cores = n_cores, 
                n_met_min = n_met_min, n_nodes = n_nodes, nf_dlat = nf_dlat, 
                nf_dlon = nf_dlon, num_jitter = num_jitter, 
                num_bg_lat = num_bg_lat, num_bg_lon = num_bg_lon, 
                num_nf_lat = num_nf_lat, num_nf_lon = num_nf_lon, 
                numpar = numpar, obs_filter = list(obs_filter), obs_fn = obs_fn,
                obs_level = obs_level, obs_path = obs_path, 
                obs_sensor = obs_sensor, obs_species = obs_species, 
                odiac_path = NA, odiac_ver = NA, projection = projection, 
                pwf_wgt = pwf_wgt, raob_fn = NA, recp_fn = recp_fn, 
                run_emiss_err = F, run_foot = run_foot, run_hor_err = F, 
                run_sim = F, run_slant = run_slant, run_trajec = run_trajec, 
                run_ver_err = F, run_wind_err = F, site = site, slurm = slurm, 
                slurm_account = slurm_account,slurm_partition = slurm_partition,
                smooth_factor = smooth_factor, store_path = store_path, 
                store_totx = store_totx, time_integrate = time_integrate, 
                time_integrate2 = list(time_integrate2), timeout = timeout, 
                timestr = list(timestr), varstrajec = varstrajec, 
                xstilt_wd = xstilt_wd, zisf = 1)
cat('Done with creating namelist...\n')

config_xstilt(namelist)  # start X-STILT, either calc traj, foot or simulation
q('no')
# end of main script
# ---------------------------------------------------------------------------- #