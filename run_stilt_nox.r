# script to estimate NOx lifetime for each STILT trajectory at each time stamp
# and calculate the chemical footprint in parallels - DW, 02/03/2021
# update using EDGARv6 for NOx, CO, and CO2 - DW, 11/17/2022
# speed up the solver time - DW, 11/29/2022
# add simulations using posterior emissions - DW, 01/06/2023

#args = commandArgs(trailingOnly = TRUE)

# CHANGE working directory ***
homedir = '/central/home/dienwu'
xstilt_wd = file.path(homedir, 'X-STILT')    # current dir
source(file.path(xstilt_wd, 'r/dependencies.r'))    # source all functions

# Please insert your API in the 'insert_ggAPI.csv' for use of ggplot and ggmap
# google API can be obtained from https://console.developers.google.com/
api.key = readLines('insert_ggAPI.csv')

### ----------------------------------------------------------------------------
# choose the site
indx = 4
epa_name = c('LosAngeles', 'Phoenix', 'Baotou', 'New Madrid', 'Thomas Hill', 
             'Martin Lake', 'Intermountain')[indx]
met = c('hrrr', 'gfs0p25')[1]
lon_lat = get.lon.lat.epa(epa_name, 0.5, 0.5, api.key) 
site   = lon_lat$site_id
epa_tz = lon_lat$tz     # will convert to standard time

### ----------------------------------------------------------------------------
input_path = '/central/groups/POW'
workdir  = '/home/dienwu/postdoc_proj/NOx'
tau_fn   = file.path(xstilt_wd, 'data/nox_curves_gapfill_v20230503.rds')
no2_path = file.path(input_path, 'TROPOMI/NO2/L2')
aux_path = file.path(input_path, 'TROPOMI/NO2/AUX')
co_path  = file.path(input_path, 'TROPOMI/CO/L2_v2')
co2_path = file.path(input_path, 'OCO-3/L2_Lite_FP_V10p4r')
epa_path = file.path(workdir, 'pp/ampd')
edgar_path = file.path(xstilt_wd, 'data')

store_path = file.path(input_path, 'XSTILT_output', site)
out_path = list.files(store_path, 'out_20', full.names = T)
out_path = out_path[grepl('TROPOMI', out_path) & grepl(met, out_path)]
no2_dirs = list.files(out_path, 'NO2', full.names = T)
all_timestr = unique(sort(strsplit.to.df(basename(dirname(no2_dirs)))$V2))


### ----------------------------------------------------------------------------
### choose overpass time and get TROPOMI NO2 and trajec info
#timestr = all_timestr[as.numeric(args[1])]
timestr = c('2020020818', '2020061518')[2]
print(timestr)

# NO, CO, CO2 emission files ----------
ym = paste0('2018', substr(timestr, 5, 6))
eno_fn = file.path(edgar_path, paste0('EDGARv6p1_NOx_', ym, '.tif'))
eco_fn = file.path(edgar_path, paste0('EDGARv6p1_CO_', ym, '.tif'))
eco2_edgar_fn = file.path(edgar_path, paste0('EDGARv6p1_CO2_', ym, '.tif'))
byid_path = file.path(out_path[grepl(timestr, out_path)], 'NO2/by-id')

# check if trajec and satellite data are available before simulations
fns_list = prep_fns4no2(timestr, byid_path, no2_path, aux_path, co_path, 
                        co2_path, lon_lat)
traj_info = all_info = fns_list$traj_info
no2_fn = fns_list$no2_fn
co2_fn = fns_list$co2_fn
co_fn = fns_list$co_fn


### ----------------------------------------------------------------------------
run_missingTF = F
ii = 1
for (i in ii) {

    # runs         1,    2,    3,    4,    5,     
    mx_res   = c(  1,    1,    1,    1,    1)[i]   # in km
    mx_hr    = c(  3,    3,    3,    3,    3)[i]   # in hr
    emiss_sf = c(  1,    1,    1,    1,    1)[i] 
    
    # NOx lifetime treatment, input NOx curves version or a constant value
    #tau_chem = 'v20230503'
    tau_chem = NA

    # background NOx in ppm, if NA, grab from Auxillary data from TROPOMI
    bg_nox = NA 

    # if EPA, scale EDGAR PP using EPA hourly data
    aq_invent  = c('epa', 'edgar', 'epa', 'odiac', 'posterior')[i]  
    ghg_invent = c('epa', 'edgar', 'epa', 'odiac', 'posterior')[i]
    
    if (ghg_invent == 'epa' | aq_invent == 'epa') {
        epa_fn = file.path(epa_path, paste0(substr(timestr, 1, 4), 
                                            tolower(lon_lat$state_abb), 
                                            substr(timestr, 5, 6), '.csv'))
        if (length(epa_fn) == 0) stop('NO hourly EPA data found...\n')
    } else epa_fn = NA

    ### ------------------------------------------------------------------------
    # parallel computing
    slurm = T
    n_nodes = 12
    n_cores = ceiling(nrow(traj_info) / n_nodes)
    if (n_cores > 10) n_cores = 10
    if (!slurm) n_nodes = n_cores = 1       # for debugging purpose
    jobname = paste0('r', i, '_', substr(timestr, 3, 8), '_', site, '_', met)
    slurm_options = list(time = '06:00:00', account = 'pow', partition = 'any',
                         mem = 10 * n_nodes * 1024)
    
    # get output filename
    rds_patt = naming_fns(aq_invent, ghg_invent, mx_res, emiss_sf, tau_chem)
    print(rds_patt)
    if (run_missingTF) {    

        # modified time needs to be after this min time
        mtime_min = '2023-04-01 00:00:00'        
        ms_info = subset_recp(byid_path, rds_patt, all_info, mtime_min) 
        if (length(ms_info) == 0) { cat('NO missing simulations...\n'); next }
        
        traj_info = ms_info
        cat(paste('\n\n** Need to rerun', nrow(traj_info), 
                'missing/outdated receptors...\n'))
        n_cores = ceiling(nrow(traj_info) / n_nodes)
        jobname = paste0(jobname, '_ms')

    } else {    # if re-run for all receptors, delete previous versions
        rds_fns = list.files(byid_path, rds_patt, full.names = T, recursive = T)
        file.remove(rds_fns)
    }   # end if 

    setwd('/home/dienwu/postdoc_proj/NOx/stilt_nox')
    xstilt_apply(FUN = xgas_solver, slurm, slurm_options, n_nodes, n_cores, 
                 jobname, site, traj_fn = traj_info$fn, 
                 timestr = traj_info$time, recp_lon = traj_info$lon, 
                 recp_lat = traj_info$lat, bg_nox, eno_fn, eco_fn, 
                 eco2_edgar_fn, NA, eno_sf = emiss_sf, eco_sf = emiss_sf, 
                 eco2_sf = emiss_sf, perturb_emissTF = F, perturb_tsTF = F, 
                 perturb_mixTF = F, perturb_fn = NA, perturb_indx = NA, 
                 no2_fn, aux_path, co_fn, co2_fn, epa_fn, epa_name, epa_tz, 
                 aq_invent, ghg_invent, nmins = -720, xmin = NA, xmax = NA, 
                 ymin = NA, ymax = NA, mx_res, mx_hr, tau_fn, tau_chem, 
                 storeTF = T, truncTF = T, xstilt_wd)

}

