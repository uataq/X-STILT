#' Compute background XCO2, given different methods:
#' M1. Trajec-endpoint (using CarbonTracker)
#' M2H. Regional daily median (based on Hakkareinen et al., 2016)
#' M2S. Localized normal statistics (based on Silva and Arellano, 2017)
#' M3. X-STILT overpass-specific (based on Wu et al., GMDD)
#' originated from 'create_namelist_oco2-xsilt.r'
#' @author Dien Wu, 04/18/2018, latest change 07/05/2021

## source all functions and load all libraries
# CHANGE working directory ***
homedir = '/central/home/dienwu'
xstilt_wd = file.path(homedir, 'models/X-STILT') 
setwd(xstilt_wd); source('r/dependencies.r') # source all functions

# Please insert your API in the 'insert_ggAPI.csv' for use of ggplot and ggmap
# google API can be obtained from https://console.developers.google.com/
api.key = readLines('insert_ggAPI.csv'); register_google(key = api.key)


#------------------------------ STEP 1 --------------------------------------- #
site      = 'SanFrancisco'  # insert target city
lon_lat   = get.lon.lat(site = site, dlon = 1, dlat = 2)
site_lon  = lon_lat$citylon
site_lat  = lon_lat$citylat

# time string in format of YYYYMMDDHH or YYYYMMDD, can be a vector
# If HH is NOT provided, one can provide the satellite path/file
# inner function will automatically look for overpass HH
all_timestr = c(2020101317, 2021020420, 2021020719, 2021022718)
#all_timestr = c(20201013, 20210204, 20210207, 20210227)

# paths
input_path  = '/central/groups/POW'
sensor      = c('OCO-2', 'OCO-3', 'TROPOMI')[3]
sensor_gas  = c('CO2', 'CO', 'NO2', 'CH4')[4]
sensor_ver  = c('V10r', 'VEarlyR', NA)[3]
oco_path    = file.path(input_path, sensor, paste0('L2_Lite_FP_', sensor_ver))
trp_path    = file.path(input_path, sensor, sensor_gas)
sensor_path = ifelse(sensor == 'TROPOMI', trp_path, oco_path)
raob_path   = file.path(input_path, 'RAOB', site)    # path for radiosonde data
store_path  = file.path(input_path, 'XSTILT_output', site)


# --------------------------------- STEP 2 ----------------------------------- #
run_trajec   = F            # T: run forward traj, will overwrite existing
run_bg       = T            # T: calculate background and plot 2D density 
                            #    requires forward trajec to be ready
run_hor_err  = T            # T: run trajec with horizontal wind errors 
run_ver_err  = F            # T: run trajec with mixed layer height scaling
run_wind_err = F            # T: run the wind error estimates using RAOB
                            # F: use the prescribed wind RMSE in m s-1, see siguverr
siguverr     = 3            # if run_wind_err == F, specify wind RMSE value, m s-1

# set zisf = 1 if run_ver_err = F
zisf = c(0.6, 0.8, 1.0, 1.2, 1.4)[3]; if (!run_ver_err) zisf = 1.0       

# release particles from a box around the city center -------------------------
# dxyp: if > 0, randomized horizontal release locations for this receptor 
#       (xp +/- dxyp, yp +/- dxyp instead of xp, yp) in units of x, y-grid lengths
# Final box length = 2 * dxyp * met.res, DW, 06/30/2020
# default box.len = 0.5 meaning 0.5 deg x 0.5 deg box around city center
box.len = 0.3                   # specify the desired box size for recp in deg

# time window for continuously release particles ------------------------------  
dtime.from = -10                # FROM 10 hours before the overpass time (-10)
dtime.to   = 0                  # TILL the overpass time (dtime.to = 0) 
dtime.sep  = 0.5                # WITH 30-min interval (dtime.sep = 0.5 hr)

# numbers of hours before terminating the generation of forward-time particles, 
# forward run, nhrs should always be positive, default is 12 hours 
nhrs = 12                

# release particles from fixed levels, recort particles for every 2mins
delt   = 2          # in mins
agl    = 10         # in mAGL, if for power plants, use stack height
numpar = 1000       # particle # per each time window

# meteo info ------------------------------------------------------------------
met_indx = 1
met      = c('gfs0p25', 'hrrr')[met_indx]
met_res  = c(0.25, 0.027)[met_indx]         # horizontal grid spacing ~in degree
met_path = file.path(homedir, met) 
met_file_format = c('%Y%m%d', '%Y%m%d')[met_indx]
cat(paste('\n\nDone with settings with receptor and meteo fields...\n'))


# -------------------------------- STEP 3 ------------------------------------ #
if (run_trajec) {       # parallel running, DW, 11/06/2019

    n_cores = 8
    n_nodes = ceiling(length(all_timestr) / n_cores)
    slurm_options = list(time = '04:00:00', account = 'pow', partition = 'any')
    jobname = paste0('XSTILT_forward_', site, '_', sensor)
    xstilt_apply(FUN = run.forward.trajec, slurm = T, slurm_options, 
                 n_nodes, n_cores, jobname, site, site_lon, site_lat, 
                 timestr = all_timestr, run_trajec, run_hor_err, run_wind_err, 
                 run_ver_err, xstilt_wd, store_path, box.len, dtime.from, 
                 dtime.to, dtime.sep, nhrs, delt, agl, numpar, met, met_res, 
                 met_file_format, met_path, raob_path, siguverr, nhrs.zisf = 24, 
                 zisf, sensor, sensor_path, sensor_ver, sensor_gas)
    q('no')
}   # end of running forward trajec


# -------------------------------- STEP 4 ------------------------------------ #
#' @param for calculating background based on 2D kernel density 
# threshold for outmost boundary of modeled urban plume (smaller td, more outwards)
td      = 0.1                            # range from 0.1 to 1 
bg_dlat = 0.5                            # buffer for bg region, in deg-lat
bg_dlon = 0.5                            # buffer for bg region, in deg-lon
zoom    = 8                              # zoom for plotting ggmap, see ggmap()
writeTF = T                              # T: output bg info to txt file
qfTF    = T                              # T: use screened obs; F: use all obs
sensor_qa = 0.5                          # quality assurance for TROPOMI, QA > 0.5

if (run_trajec == F & run_bg) {         # need forward trajec ready

    bg_df = NULL
    for ( tt in 1 : length(all_timestr) ) {            # loop over each overpass
        timestr = all_timestr[tt]
        tmp_df = calc.bg.forward.trajec(site, timestr, sensor, sensor_path, 
                                        sensor_gas, sensor_ver, sensor_qa, qfTF,
                                        store_path, met, td, bg_dlat, bg_dlon, 
                                        zoom, api.key)
        if ( is.null(tmp_df) ) next
        bg_df = rbind(bg_df, tmp_df)
    }   # end for tt

    if (writeTF) {
        if ( grepl('OCO', sensor) ) {           # label qfTF for OCO
            fn = paste0('bg_', site, '_', sensor, '_', sensor_gas, '_qf', qfTF, '.txt')
        } else if (sensor == 'TROPOMI') {       # label qa for TROPOMI
            fn = paste0('bg_', site, '_', sensor, '_', sensor_gas, '_qa', sensor_qa, '.txt')
        } else stop('Invalid sensor names...\n'); print(fn)

        write.table(bg_df, file.path(store_path, fn), quote = F, row.names = F, sep = ',')
    }   # end if writeTF
}   # end if run_bg

# end of script