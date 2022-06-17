# script to load OCO-2 or OCO-3 or TROPOMI CO, NO2 data, DW, 03/04/2021
# will add more sensors or species soon 

#' @param timestr needs to be in form of YYYYMMDDHH
#' @param qa quality assurance of TROPOMI data
#'           e.g., if qa = 0.4, chose data with QA >= 0.4
if (F) {
    sensor = obs_sensor
    sensor_gas = obs_species
    tropomi.fn = sensor_fn = obs_fn
    tropomi.path = sensor_path = obs_path
}

# for TROPOMI vertical column density (VCD), convert to mixing ratio in ppb 
#' get rid of @param sensor.ver for grabbing warn levels, DW, 07/23/2021
get.sensor.obs = function(site, timestr, 
                          sensor = c('OCO-2', 'OCO-3', 'TROPOMI')[1], 
                          sensor_gas = c('CO2', 'CO', 'NO2', 'CH4', 'HCHO')[1], 
                          sensor_fn = NULL, sensor_path, qfTF = T, 
                          tropomi_qa = c(0, 0.4, 0.5, 0.7, 1)[1], 
                          lon_lat = NULL) {

    # --------------------------- STEP 1 ------------------------------------- #
    # grab satellite observations (OCO-2/3, TROPOMI CO and NO2)
    if ( is.null(lon_lat) ) lon_lat = get.lon.lat(site, dlat = 1.5, dlon = 1.5)
    err_message = 'get.sensor.obs(): NO observation found...\n'
    err_message2 = 'get.sensor.obs(): NO QUANLIFIED observation found...\n'

    if ( grepl('OCO', sensor) ) {    # OCO-2 or 3 CO2
        obs = grab.oco(oco.path = sensor_path, timestr, lon_lat, 
                       oco.fn = sensor_fn) 
        if ( is.null(obs) ) { cat(err_message); return() }
        if ( nrow(obs) == 0 ) { cat(err_message); return() }
        obs$datestr = as.POSIXlt(as.character(obs$time), 'UTC', 
                                 format = '%Y-%m-%d %H:%M:%S')
        if (qfTF) obs = obs %>% filter(qf == 0)
    }
    
    # TROPOMI CO VCD and mixing ratio XCO in ppb
    if ( grepl('TROPOMI', sensor) & sensor_gas == 'CO' ) {      
        obs = grab.tropomi.co(tropomi.path = sensor_path, timestr, lon_lat, 
                              tropomi.fn = sensor_fn) 
        if ( is.null(obs) ) { cat(err_message); return() }
        if ( nrow(obs) == 0) { cat(err_message); return() }
        obs$datestr = as.POSIXlt(as.character(obs$time_utc), 'UTC', 
                                 format = '%Y%m%d%H%M%S')
        if (qfTF) obs = obs %>% filter(qa >= tropomi_qa)
    }
    
    # TROPOMI NO2, final mixing ratio of tropo NO2 in ppb
    if ( grepl('TROPOMI', sensor) & sensor_gas == 'NO2' ) {      
        obs = grab.tropomi.no2(tropomi.path = sensor_path, timestr, lon_lat, 
                               tropomi.fn = sensor_fn) 
        if ( is.null(obs) ) { cat(err_message); return() }
        if ( nrow(obs) == 0) { cat(err_message); return() }
        obs$datestr = as.POSIXlt(as.character(obs$time_utc), 'UTC', 
                                 format = '%Y%m%d%H%M%S')
        if (qfTF) obs = obs %>% filter(qa >= tropomi_qa)
    }

    # TROPOMI CH4
    if ( grepl('TROPOMI', sensor) & sensor_gas == 'CH4' ) {      
        obs = grab.tropomi.ch4(tropomi.path = sensor_path, timestr, lon_lat, 
                               tropomi.fn = sensor_fn) 
        if ( is.null(obs) ) { cat(err_message); return() }
        if ( nrow(obs) == 0) { cat(err_message); return() }
        obs$datestr = as.POSIXlt(as.character(obs$time_utc), 'UTC', 
                                 format = '%Y%m%d%H%M%S')
        if (qfTF) obs = obs %>% filter(qa >= tropomi_qa)
    }

    # TROPOMI HCHO
    if ( grepl('TROPOMI', sensor) & sensor_gas == 'HCHO' ) {      
        obs = grab.tropomi.hcho(tropomi.path = sensor_path, timestr, lon_lat, 
                                tropomi.fn = sensor_fn) 
        if ( is.null(obs) ) { cat(err_message); return() }
        if ( nrow(obs) == 0) { cat(err_message); return() }
        obs$datestr = as.POSIXlt(as.character(obs$time_utc), 'UTC', 
                                 format = '%Y%m%d%H%M%S')
        if (qfTF) obs = obs %>% filter(qa >= tropomi_qa)
    }


    if ( is.null(obs) ) { cat(err_message2); return() }
    if ( nrow(obs) == 0 ) { cat(err_message2); return() }

    # initiate swath number as 1 and add row number, DW, 03/09/2021
    obs = obs %>% mutate(swath = 1, row = seq(1, nrow(obs)))  
    if ( grepl('OCO', sensor) ) { 

        # assign polygon indx to OCO-3 SAM
        zero.indx   = which(obs$vertices == 1)
        obs$polygon = findInterval(obs$row, zero.indx)
     
        # if SAM, differentiate between swaths, DW, 08/04/2020
        if ('Land_SAM' %in% obs$mode ) {

            # figure out the different tracks using time and latitude difference
            diff.lat   = abs(diff((obs %>% arrange(time))$lat))
            track.indx = c(0, which(diff.lat > 0.5)) + 1
            track.indx = track.indx[track.indx != nrow(obs)]
            obs = obs %>% mutate(swath = findInterval(obs$row, track.indx))
        }   # end SAM

    } else if ( grepl('TROPOMI', sensor) ) {        

        # assign fake polygon indx
        zero.indx   = which(obs$corner == 0)
        obs$polygon = findInterval(obs$row, zero.indx)
    }   # end if

    obs = obs %>% dplyr::select(-row)
    return(obs)
}
