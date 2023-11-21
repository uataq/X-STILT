
# main script to simulate tropospheric NO2 and column CO and CO2
# by Dien Wu

# in order to use this script, you will need to generate STILT trajectories 

# ---------------------------
# archived updates: 
# script to estimate NOx lifetime for each STILT trajectory at each time stamp
# and calculate the chemical footprint in parallels - DW, 02/03/2021
# update emission using EDGARv6 for NOx, CO, and CO2 - DW, 11/17/2022
# speed up the solver time - DW, 11/29/2022
# add simulations using posterior emissions - DW, 01/06/2023
# add ODIAC and Vulcan CO2 emissions & refactoring - DW, 07/21/2023 
# modify mixing lenth-/time-scale (typical Dh = 1e3 to 1e4 m2 s-1) - DW, 07/25/2023

#args = commandArgs(trailingOnly = TRUE)

# CHANGE working directory ***
homedir = '/central/home/dienwu'
xstilt_wd = file.path(homedir, 'X-STILT')    # current dir
source(file.path(xstilt_wd, 'r/dependencies.r'))    # source all functions

# Please insert your API in the 'insert_ggAPI.csv' for use of ggplot and ggmap
# google API can be obtained from https://console.developers.google.com/
api.key = readLines('insert_ggAPI.csv')

# ----------------------------- choose the site ----------------------------- #
site = 'LosAngeles'
epa_name = site     # only useful if simulating over power plants using EPA emis
met  = 'hrrr'
dlon = 1.5
dlat = 1

# - automatically obtain lat/lon (requires Google API key in insert_ggAPI.csv)
# - OR mannually insert coordinates, e.g., site.loc = data.frame(lon = , lat = )
lon_lat = get.lon.lat(site, dlon, dlat, site.loc = NULL, api.key)

# ------------------------------ I/O path ----------------------------------- #
# including main working dir, satellite tNO2, xCO, xCO2 paths, and output path
input_path = '/central/groups/POW'
workdir    = '/home/dienwu/postdoc_proj/NOx'
trp_path   = file.path(input_path, 'TROPOMI')
tno2_path  = file.path(trp_path, 'NO2/L2')
tno2x_path = file.path(trp_path, 'NO2/AUX')
xco_path   = file.path(trp_path, 'CO/L2')
xch4_path  = file.path(trp_path, 'CH4/L2')
xco2_path  = file.path(input_path, 'OCO-3/L2_Lite_FP_V10p4r')
store_path = file.path(input_path, 'XSTILT_output', site)
all_timestr = get_timestr_from_outpath(store_path, met)

# emission path
edgar_path = file.path(xstilt_wd, 'data/EDGARv6')
epa_path   = file.path(workdir, 'pp/ampd')
vulc_path  = file.path(input_path, 'XSTILT_dep/Vulcan/data')
odiac_path = file.path(input_path, 'XSTILT_dep/ODIAC/v2022/2020') 

# -------------------------- choose overpass time --------------------------- #
#timestr = all_timestr[as.numeric(args[1])]
timestr = '2022012120'
print(timestr)

# ------------------ NO, CO, CO2 emission parameters ------------------------ #
aq_invent  = c('epa', 'edgar', 'odiac', 'vulcan', 'posterior')[2]  
ghg_invent = c('epa', 'edgar', 'odiac', 'vulcan', 'posterior')[2]
edgar_ym   = paste0('2018', substr(timestr, 5, 6))
eno_fn  = file.path(edgar_path, paste0('EDGARv6p1_NOx_', edgar_ym, '.tif'))
eco_fn  = file.path(edgar_path, paste0('EDGARv6p1_CO_',  edgar_ym, '.tif'))
eco2_fn = file.path(edgar_path, paste0('EDGARv6p1_CO2_', edgar_ym, '.tif'))
ech4_fn = file.path(edgar_path, paste0('EDGARv6p1_CH4_', edgar_ym, '.tif'))
tno2_path = file.path(tno2_path, substr(timestr, 1, 4))

# a second non-EDGAR inventory for CO2 emissions (e.g., vulcan for 2015, odiac for current year)
# will adopt emission ratio from EDGAR to approximate CO and NOx emissions
eco2_fn2 = NA
#eco2_fn2 = file.path(vulc_path, 'Vulcan_v3_US_annual_1km_total_mn.nc4')
#eco2_fn2 = file.path(odiac_path, 'odiac2022_1km_excl_intl_2001.tif')  

# --------------------- Chem & emiss & mixing parameters --------------------- #
overwriteTF = T
ts_chem = NA    # if NA, use NOx curve; if not, set a constant lifetime [hr]
ts_fn   = file.path(xstilt_wd, 'data/nox_curves_gapfill_v20230503.rds')
bg_nox  = NA     # if NA, grab from TROPOMI aux data; else bg NOx in ppm
eno_sf  = 1      # scaling factor of emissions for the whole field
eco_sf  = 1 
ech4_sf = 1
eco2_sf = 1     
mx_res  = 1      # mixing box in km, default of 1km
mx_hr   = 1       # mixing timescale in hr

# ------------------------- Perturbation parameters ------------------------- #
# this is useful for quantifying uncertainty or emission inversion via EnKF
perturb_emiTF = F
perturb_mixTF = F
perturb_tsTF  = F
n_perturb  = c(NA, 30)[1]   # number of perturbation
transerrTF = F 

# ----------------------------- Parallel computing -------------------------- #
slurm = T
n_nodes = 12
n_cores = 10
mem_per_node = 10
if (!slurm) n_nodes = n_cores = 1
slurm_path = '/home/dienwu/postdoc_proj/NOx/stilt_nox'
namelist = list(aq_invent = aq_invent, bg_nox = bg_nox, ech4_fn = ech4_fn,  
                eco_fn = eco_fn, eco2_fn = eco2_fn, eco2_fn2 = eco2_fn2, 
                eno_fn = eno_fn, ech4_sf = ech4_sf, eco_sf = eco_sf, 
                eco2_sf = eco2_sf, eno_sf = eno_sf, epa_name = epa_name, 
                epa_path = epa_path, ghg_invent = ghg_invent, 
                lon_lat = list(lon_lat), met = met, mem_per_node = mem_per_node,
                mx_hr = mx_hr, mx_res = mx_res, n_cores = n_cores, 
                n_nodes = n_nodes, n_perturb = n_perturb, 
                overwriteTF = overwriteTF, perturb_emiTF = perturb_emiTF, 
                perturb_mixTF = perturb_mixTF, perturb_tsTF = perturb_tsTF, 
                site = site, slurm = slurm, slurm_path = slurm_path, 
                store_path = store_path, timestr = timestr, 
                tno2_path = tno2_path, tno2x_path = tno2x_path, 
                transerrTF = transerrTF, ts_chem = ts_chem, ts_fn = ts_fn, 
                xch4_path = xch4_path, xco_path = xco_path, 
                xco2_path = xco2_path, xstilt_wd = xstilt_wd)

config_chemv2(namelist)
