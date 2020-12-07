#' A lite and ideal version of X-STILT without any dependencies on satellites. 
#' PLEASE use with cautions. PLEASE read the instructions after *** 
#' @author Dien Wu, 12/06/2020

#' A. run_trajec == T or run_foot == T: run trajectory or footprint, simulations
#' B. run_trajec == F & run_foot == F & run_sim = T: requires trajec and foot
#'    X-STILT models XCO2 enhancements based on FFCO2 emission from 1km ODIAC;

#' ---------------------------------------------------------------------------
#' Instructions for ideal column simulations 
#' *** AK is simply treated as 1 (@param ak_wgt == FALSE)
#' *** PW is calculated using modeled met variables (@param pwf_wgt == TRUE)

#' *** The user NEEDS to provide a txt or csv file for receptor locations
#' *** correct file should include long, lat, and time (optional, YYYYMMDDHH). 

#' *** You can also assign your receptor time to @param timestr. HOWEVER it will 
#' *** overwrite the existing time column in the receptor file if there is. 
#' ---------------------------------------------------------------------------


## go to X-STILT dir and source functions and libraries
homedir   = '/central/home/dienwu'
xstilt_wd = file.path(homedir, 'models/X-STILT') 
setwd(xstilt_wd); source('r/dependencies.r')

# *** there is a current bug with the rslurm package, you may need to 
# mannually pull the development version from github by doing:
#devtools::install_github('SESYNC-ci/rslurm')               # DW, 11/6/2020

# ------------------------------ Basis XSTILT flags -------------------------- #
site = 'Paris'

# as a default, all input data are stored under input_path
input_path = '/central/groups/POW'
store_path = file.path(input_path, 'XSTILT_output', site)
odiac_ver  = c('2016', '2017', '2018', '2019')[4]         # ODIAC version
odiac_path = file.path(input_path, 'ODIAC', paste0('ODIAC', odiac_ver))  
raob_path  = file.path(input_path, 'RAOB', site)  # NOAA radiosonde

run_trajec    = T     # whether to generate trajec; T: overwrite existing trajec
run_foot      = T     # whether to generate footprint
run_hor_err   = F     # T: need to set error parameters
run_ver_err   = F     # T: need to set error parameters
run_emiss_err = F     # T: get XCO2 error due to prior emiss err
run_sim       = F     # whether to run XCO2.ff or error using existing foot
nhrs          = -24   # number of hours backward (-) or forward (+)

# output variable names in trajec.rds
varstrajec = c('time', 'indx', 'lati', 'long', 'zagl', 'zsfc', 'foot', 'samt',
               'dmas', 'mlht', 'temz', 'pres', 'sigw', 'tlgr', 'dens', 'sphu')


# ------------------------- ARL format meteo params -------------------------- #
# line source for column release according to HYSPLITv5, DW, 09/17/2020
# all particles are roughly evenly distributed between minagl and maxagl
minagl = 0           # min release height in meter AGL
maxagl = 3000        # max release height in meter AGL
numpar = 3000        # total number of particles between minagl and maxagl

# A csv or txt file containing desired longitude ('long') and latitude ('lati'), 
# one can also include a third column for receptor time ('time') OR assign the 
# receptor time to `timestr` below. Each column needs to be separated by comma. 
recp_file = 'receptor.csv'    # <path/filename>

# if you assign `timestr` here, it'll overwrite the time column in `recp_file`; 
# if you don't want the receptor time to be overwrited, leave it as NA
timestr = '2020041211'      # YYYYMMDDHH in UTC, as character string


# -------------------------- ARL format meteo params ------------------------- #
# see https://uataq.github.io/stilt/#/configuration?id=meteorological-data-input
met       = c('gdas0p5', 'gfs0p25', 'hrrr')[2]   # choose met fields
met_path  = file.path(homedir, met)              # path of met fields
n_met_min = 1                                    # min number of files needed
met_file_format = '%Y%m%d'                       # met file name convention

# OPTION for subseting met fields if met_subgrid_enable is on, 
# useful for large met fields like GFS or HRRR
met_subgrid_buffer = 0.1   # Percent to extend footprint area for met subdomain
met_subgrid_enable = T    

# if set, extracts the defined number of vertical levels from met files
met_subgrid_levels = NA    # default as NA


# -------------------------- Transport error params -------------------------- #
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


# ----------------------------- Footprint params ---------------------------- #
# define spatial domain for final spatial footprint, with params in degrees
# e.g., long range of footprint [site.lon - foot_dlon, site.lon + foot_dlon]; 
#       lati range of footprint [site_lat - foot_dlat, site_lat + foot_dlat]
site_lon  = 2.35
site_lat  = 48.86
foot_dlat = 10      # default as 10, meaning 20 x 20deg box around city
foot_dlon = 10 

# footprint resolution in degree, 1km for ODIAC
foot_res  = 1/120  
foot_nhrs = nhrs   # if foot_nhrs < nhrs, subset trajec and then calc footprint

# (optinal) generate footprint with different resolutions other than "foot_res"
# foot_res2 can be a vector, foot filename will contain res info, DW, 02/11/2019
# if no need to generate second sets of footprints, set it to NA     
foot_res2 = c(NA, 1/10, 1/20, 1)[1]     

# whether weight footprint by AK and PW for column simulations
# for ideal simulation, set ak.wgt as FALSE; otherwise will look for OCO-2 files
ak_wgt  = FALSE               # if FALSE, AK will be set as 1 for regional sims
pwf_wgt = TRUE

# other neccesary footprint params using STILTv2 (Fasoli et al., 2018)
hnf_plume      = TRUE  # whether turn on hyper near-field (HNP) for mising hgts
smooth_factor  = 1     # Gaussian smooth factor, 0 to disable
time_integrate = TRUE  # whether to integrate footprint along time; F, hourly foot
projection     = '+proj=longlat'
cat('Done with params for error analysis and footprint setup...\n')


# ----------------------------- SLURM params -------------------------------- #
# avoid using < 10 cores per node when running trans error stat (run_hor_err) 
# along with calculating 2D foot together (run_foot), *** memory limits
n_nodes = 10
n_cores = 10
timeout = 12 * 60 * 60    # time allowed before forced terminations in second
job_time = '12:00:00'     # total job time
slurm = T                 # if TRUE, for SLURM parallel computing
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
                hnf_plume = hnf_plume, job_time = job_time, 
                mem_per_node = mem_per_node, met = met, 
                met_file_format = met_file_format, met_path = met_path, 
                met_subgrid_buffer = met_subgrid_buffer, 
                met_subgrid_enable = met_subgrid_enable, 
                met_subgrid_levels = met_subgrid_levels, 
                minagl = minagl, maxagl = maxagl, nhrs = nhrs, n_cores = n_cores, 
                n_met_min = n_met_min, n_nodes = n_nodes, numpar = numpar, 
                odiac_ver = odiac_ver, odiac_path = odiac_path, 
                projection = projection, pwf_wgt = pwf_wgt, recp_file = recp_file, 
                run_foot = run_foot, run_hor_err = run_hor_err, run_sim = run_sim, 
                run_trajec = run_trajec, run_ver_err = run_ver_err, site = site, 
                site_lon = site_lon, site_lat = site_lat, slurm = slurm, 
                slurm_account = slurm_account, slurm_partition = slurm_partition, 
                smooth_factor = smooth_factor, store_path = store_path, 
                timestr = timestr, time_integrate = time_integrate, 
                timeout = timeout, varstrajec = varstrajec, 
                xstilt_wd = xstilt_wd, zisf = zisf)      

cat('Done with creating namelist...\n\n')

config_xstilt_ideal(namelist)  # start X-STILT, either calc traj, foot or simulation
# end of main script
# ---------------------------------------------------------------------------- #