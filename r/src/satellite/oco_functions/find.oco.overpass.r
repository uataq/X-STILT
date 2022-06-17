# ---------------------------------------------------------------------------- #
# function that search for any OCO-2/3 overpases for a given region, DW, 05/15/2017
#' @param date.range date range, c(YYYYMMDD_1, YYYYMMDD_2)
#' @param oco.path default in lin-group
#' @param oco.ver needs to contain version number
#' @param lon.lat a data frame containining c(minlat, maxlat, minlon, maxlon), 
#'                please follow THIS ORDER!!!
#' @param nfTF for searching soundings near near-field region, DW, 06/15/2018
# ---------------------------------------------------------------------------- #
# add count for good quality data, DW, 12/20/2017
# update for v9 data, DW, 10/19/2018 
# drop the scientific notation for sounding ID, DR, DW, 09/04/2019
# update for OCO-3 Vearly data, DW, 06/29/2020 
# remove warn levels, DW, 07/01/2020 
# ---------------------------------------------------------------------------- #

find.oco.overpass = function(date.range, lon.lat, oco.ver = 'V10r', oco.path, 
                             nfTF = F, nf.dlon = 0.5, nf.dlat = 0.5){ 
  
  # find.overpass() in r/src/satellite/oco_functions
  library(geosphere); library(ncdf4); library(dplyr)

  # path and filename for storing OCO-2 info
  all.fns = list.files(pattern = 'LtCO2_', path = oco.path)
  if (length(all.fns) == 0) stop('find.oco.overpass(): NO OCO data found...check file path\n')

  # get rid of some characters
  file.info = gsub('.nc4', '', all.fns)
  file.info = strsplit.to.df(file.info)
  all.timestr = file.info$V3
  if (grepl('QTS', oco.path)) all.timestr = file.info$V5 

  if (grepl('7', oco.ver)) {   # for version 7r
    all.timestr[nchar(all.timestr) == 6] = 
      paste0('20', all.timestr[nchar(all.timestr) == 6])

    oco.file = all.fns[all.timestr >= date.range[1] & 
                       all.timestr <= date.range[2]]

    timestr = all.timestr[all.timestr >= date.range[1] & 
                          all.timestr <= date.range[2]]

  } else if (!grepl('7', oco.ver)) {   # for version 8, 9r or even for OCO-3
    SEL.day = all.timestr >= substr(date.range[1], 3, 8) &
              all.timestr <= substr(date.range[2], 3, 8)
    oco.file = all.fns[SEL.day]
    timestr  = paste0('20', all.timestr[SEL.day])
  } 

  # loop over each overpass
  result = NULL
  for (f in 1 : length(oco.file)) {
    
    if (f %% 25 == 0)
    cat(paste('#--', signif(f / length(oco.file), 3) * 100, '% SEARCHED --#\n'))
    
    dat = nc_open(file.path(oco.path, oco.file[f]))

    ## grabbing OCO-2 levels, lat, lon
    oco.lev = ncvar_get(dat, 'levels')
    oco.lat = ncvar_get(dat, 'latitude')
    oco.lon = ncvar_get(dat, 'longitude')
    xco2 = ncvar_get(dat, 'xco2'); xco2[xco2 == -999999] = NA
    qf = ncvar_get(dat, 'xco2_quality_flag')

    # drop the scientific notation for sounding ID, DR, DW, 09/04/2019
    id = format(ncvar_get(dat, 'sounding_id'), scientific = F)
    orbit = ncvar_get(dat, 'Sounding/orbit')

    # 0=Nadir, 1=Glint, 2=Target, 3=Transition, 4=Snapshot Area Map
    mode = ncvar_get(dat, 'Sounding/operation_mode')
    obs = data.frame(lat = as.numeric(oco.lat), lon = as.numeric(oco.lon), 
                     qf = as.numeric(qf), xco2 = as.numeric(xco2), 
                     orbit = as.numeric(orbit), mode = mode,
                     timestr = as.numeric(substr(id, 1, 10)), 
                     stringsAsFactors = F) 
    
    obs = obs %>% filter(lat >= lon.lat$minlat & lat <= lon.lat$maxlat &
                         lon >= lon.lat$minlon & lon <= lon.lat$maxlon) 
    nc_close(dat)
    
    # store results if there are soundings over
    if (nrow(obs) > 0) {
      
      # if we have multiple overpass hours in a day, 
      # group by orbit # and get the minimal time string per orbit, DW, 03/12/2021
      # dplyr::count(x, ..., wt = NULL (the default), ...)
      # counts the number of rows in each group.
      tot.count = obs %>% count(orbit)
      qf.count  = obs %>% filter(qf == 0) %>% count(orbit)
      sam.count = obs %>% filter(mode == 4) %>% count(orbit)
      uni.timestr = obs %>% group_by(orbit) %>% 
                    dplyr::summarise(timestr = as.character(min(timestr)), .groups = 'drop')

      # merge all df, DW, 03/12/2021
      tmp = full_join(uni.timestr, 
                      tot.count %>% rename(tot.count = n), by = 'orbit') %>% 
            left_join(sam.count %>% rename(sam.count = n), by = 'orbit') %>% 
            left_join(qf.count %>% rename(qf.count = n), by = 'orbit')

      # also search for soundings near city center,
      if (nfTF) {
        nf.dat = obs %>% filter(lon >= (lon.lat$site_lon - nf.dlon),
                                lon <= (lon.lat$site_lon + nf.dlon),
                                lat >= (lon.lat$site_lat - nf.dlat),
                                lat <= (lon.lat$site_lat + nf.dlat))
        
        tot.count2 = nf.dat %>% count(orbit)
        qf.count2  = nf.dat %>% filter(qf == 0) %>% count(orbit)
        sam.count2 = nf.dat %>% filter(mode == 4) %>% count(orbit)
        uni.timestr2 = nf.dat %>% group_by(orbit) %>% 
                       dplyr::summarise(timestr = as.character(min(timestr)), .groups = 'drop')

        # merge all df, DW, 03/12/2021
        tmp.nf = uni.timestr2 %>% 
                 left_join(tot.count2 %>% rename(tot.nf.count = n), 
                           by = 'orbit') %>%  
                 left_join(qf.count2 %>% rename(qf.nf.count = n), 
                           by = 'orbit') 

        # combine 
        tmp = tmp %>% left_join(tmp.nf, by = c('timestr', 'orbit'))
      } # end if nfTF

      tmp = tmp %>% rename(orbit.id = orbit) %>% 
            mutate(sam.count = ifelse(is.na(sam.count), 0, sam.count))

      result = rbind(result, tmp)
    } else next  # if no data, jump to next file
    # end if tot.count
  }  # end for f

  result = result[order(result$timestr),]
  return(result)
}

