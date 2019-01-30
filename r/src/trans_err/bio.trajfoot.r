#' subroutine to readin CarbonTracker-NeatRealTime for OCO-2 project
#' and calculate the dCO2 for each trajec for both biospheric signal)
#' @author: Dien Wu, 04/11/2017

#' @variables:
#' default ct.hr = 0, 3, 6, 9, 12, 15, 18, 21 UTC, left-edge hr, not centered hr
#' default variable name when reading as raster, 'bio_flux_opt' used in CT

#' @updates:
#' fix a bug, (all footprint columns are zero), 05/11/2017
#' add dmassTF for weighting footprint if violate mass conservation, DW, 10/20/2017
#' generalizt to merge with Ben's code, DW, 07/23/2018
#' read bio fluxes as raster layer that is more efficient, DW, 07/23/2018
#' add stilt.ver, if v = 1, convert colnames--
#' remove dmass correction
#' remove zero footprint to save time and space, DW, 01/29/2019 

bio.trajfoot <- function(trajdat, timestr, ctflux.path, ct.hr = seq(0, 21, 3),
                         varname = 'bio_flux_opt'){

  library(raster)

  # compute the total indx before any operations
  tot.p <- data.frame(indx = unique(trajdat$indx))

  # remove zero footprint to save time and space, DW, 01/29/2019 
  trajdat <- trajdat %>% filter(foot > 0)

  #### NOW grab a CO2.bio fluxes for each selected trajec ###
  # match 3-houly footprint and 3-houly biosperhic fluxes
  # first locate the date and hour for each backwards hours of STILT footprint
  recp.date <- as.POSIXlt(as.character(timestr), format = '%Y%m%d%H', tz = 'UTC')
  traj.date <- recp.date + trajdat$time * 60 # convert to seconds

  # get day, mon, year for trajec time
  trajdat <- trajdat %>%
             mutate(timestr = as.numeric(gsub('-', '', substr(traj.date, 1, 10))),
                    yr  = as.numeric(substr(timestr, 1, 4)),
                    mon = as.numeric(substr(timestr, 5, 6)),
                    day = as.numeric(substr(timestr, 7, 8)),
                    hr  = as.numeric(substr(traj.date, 12, 13)))
  uni.date <- unique(trajdat$timestr)

  # then match time to 3 hourly ct and get ct timestr for each particle,
  # find the correct hr and date.hr string given 3 hourly ct
  trajdat <- trajdat %>% mutate(match.hr = ct.hr[findInterval(hr, ct.hr)],
                                match.date.hr = paste0(timestr, 
                                                formatC(match.hr, width = 2, 
                                                        flag = 0)))
  uni.date.hr <- unique(trajdat$match.date.hr)
  uni.hr <- as.numeric(substr(uni.date.hr, 9, 10))

  # raster band, e.g., 1 for left-edge hr 00UTC (or centered time 1:30UTC in CT)
  bands <- findInterval(uni.hr, ct.hr)

  # loop over all combinations
  new.trajdat <- NULL
  for (d in 1 : length(uni.date.hr)) {

    # open the daily file just once
    ctfile <- list.files(path = ctflux.path, pattern = substr(uni.date.hr[d], 1, 8))
    if (length(ctfile) == 0) 
      stop('bio.trajfoot(): NO CT fluxes available for this overpass...\n')

    # read as raster, find the correct band from 'bands'
    bio <- raster(file.path(ctflux.path, ctfile), band = bands[d], 
                  varname = varname)

    # find corresponding 1x1deg CT grid for biospheric fluxes
    sel.trajdat <- trajdat %>% filter(match.date.hr == uni.date.hr[d])
    trajcor <- sel.trajdat %>% dplyr::select(long, lati)

    # get bio fluxes (convert to umol/m2/s) and finally get co2.bio
    sel.trajdat <- sel.trajdat %>% 
                   mutate(find.bio = raster::extract(x = bio * 1E6, y = trajcor), 
                          co2 = find.bio * foot)

    # store trajdat
    new.trajdat <- rbind(new.trajdat, sel.trajdat)
  } # end for d

  # also, compute total dCO2 for each traj over all backwards hours
  sum.trajdat <- new.trajdat %>% group_by(indx) %>%
                                 dplyr::summarize(bio.sum = sum(co2))

  # since we removed particles with zero footprint, we need to add zero to bio.sum 
  # to let `sum.trajdat` have the same amount of total particles as `tot.p`, 
  # DW, 01/29/2019 
  sum.trajdat <- sum.trajdat %>% right_join(tot.p, by = 'indx') 

  # NA will show up in above `sum.trajdat` for particles with foot = 0
  # thus, replace NA with 0 
  sum.trajdat$bio.sum[is.na(sum.trajdat$bio.sum)] <- 0 

  return(sum.trajdat)
} # end of subroutine
