### subroutine to compute SIGUVERR,
# require txt file with modeled vs. observed wind data
# returning met-raob comparisons and wind stat for trans error
# DW

# pre-requirement: wind error statistics generated in txtfile *****
# generalize the code for other cities by introducing 'lon.lat', DW, 07/19/2018
# nhrs: selecting ? days backward for computing wind err, DW, 07/19/2018
# add agl for released heights in mAGL, DW, 07/31/2018
# add errors in wind direction, DW, 08/29/2018
# update variables according to output from cal.met.wind(), DW, 08/31/2018

get.siguverr <- function(met.raob, nfTF = F, forwardTF = F, lon.lat = NULL,
                         nhrs = NULL, agl = c(0, 3000), timestr = NULL){

  # get max time, ie., receptor time
  max.time <- as.numeric(substr(timestr, 1, 8))

  # filter pressure levels
  met.raob <- met.raob %>% na.omit() %>% filter(pres >= 300) %>%
              mutate(time = as.numeric(substr(timestr, 1, 8)))

  # select nearfield radiosondes
  # also true for forward-time definition of background
  if (nfTF | forwardTF) {

    # 1. select wind comparisons within 2x2 deg around the city center
    # select 1day comparisons
    nf <- met.raob %>% filter(time >= max.time,
                              lon >= (lon.lat$minlon - 1), 
                              lon <= (lon.lat$minlon + 1),
                              lat >= (lon.lat$minlat - 1), 
                              lat <= (lon.lat$minlat + 1))

    # 2. if there are no data given above selection,
    # expanding time range first, then spatial domain
    if (nrow(nf) == 0) {
      cat('***** Missing nf raob, Extending time range *****\n')

      # expanding time ranges
      uni.time <- unique(met.raob$time)
      for (l in 1 : length(uni.time)) {
        nf <- met.raob %>% filter(time >= uni.time[l],
                                  lon >= (lon.lat$minlon - 1), 
                                  lon <= (lon.lat$minlon + 1),
                                  lat >= (lon.lat$minlat - 1), 
                                  lat <= (lon.lat$minlat + 1))
        if (nrow(nf) == 0) { nf <- NA; next }
      } # end loop l

      # expanding horizontal lat lon ranges
      if (is.na(nf)) {
        cat('***** Missing raob, Extending spatial domain *****\n')

        # expanding spatial domain, now 8x8 deg around city center
        nf <- met.raob %>% filter(time >= max.time,
            lon >= (lon.lat$minlon - 4), lon <= (lon.lat$minlon + 4),
            lat >= (lon.lat$minlat - 4), lat <= (lon.lat$minlat + 4))
        if (nrow(nf) == 0) {
          cat('*** Very POOR RAOB network for this event; returning NA ***\n')
          return()
        }  # end if
      } # end if nf
    } # end if

    sel.met.raob <- nf

  } else {    #### if not nearfield wind statistics

    # we have 5 days statistics, only grab first 3 days backward
    # if backward, ndays < 0, DW, 01/23/2019
    date1 <- as.POSIXct(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
    date2 <- date1 + nhrs * 60 * 60
    min.date <- min(date1, date2)
    max.date <- max(date1, date2)

    # select met.raob based on time and mAGL as well, DW, 07/31/2018
    sel.met.raob <- met.raob %>% 
                    mutate(date = as.POSIXct(as.character(timestr), 'UTC', 
                                            format = '%Y%m%d%H'), 

                          # calculate the distance between raob and city, in km
                          dist = as.numeric(fields::rdist.earth(
                                              x1 = lon.lat[c('citylon', 'citylat')], 
                                              x2 = met.raob[c('lon', 'lat')], 
                                              miles = F))) %>% 

                    # remove wind error outlier, abs() > 40 m/s and far away raob
                    # select wind errors based on nhrs and recep hgts
                    filter(abs(u.err) <= 40, abs(v.err) <= 40, dist <= 2000, 
                           date >= min.date, date <= max.date, 
                           hgt >= min(unlist(agl)), hgt <= 3000) %>% #max(unlist(agl))) %>% 
                    na.omit()

  } # end if 

  # finally calculate meawind errors and biases
  if (length(sel.met.raob) > 0) {
    siguverr <- sqrt(mean(c(sel.met.raob$u.err^2, sel.met.raob$v.err^2)))
    u.bias   <- mean(sel.met.raob$u.err)
    v.bias   <- mean(sel.met.raob$v.err)
    ws.bias  <- mean(sel.met.raob$ws.err)
    wd.bias  <- mean(sel.met.raob$wd.err)

    stat <- data.frame(siguverr, u.bias, v.bias, ws.bias, wd.bias)
    stat
  } else return()

}  # end of subroutine
