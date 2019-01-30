# subroutine to readin ODIAC emissions, v2015a and then match with footprint
# column from STILT trajs, use ak*pw weighted traj for footprint,
# find out the anthro emissions at the gridcell where trajdaticle falls into
# written by DIEN WU, 01/08/2017

# add ODIACv2016, 03/15/2017
# fix a bug, (all footprint columns are zero), 05/11/2017
# add dmassTF for weighting footprint if violate mass conservation,
#    but probably do not need to turn on, as dmass fix can over-correct foot,
#    10% violation is a good value, DW, 10/20/2017

# cut particles beyond emission grid, fix bug for summertime track, DW, 12/04/2017
# generalizt to merge with Ben's code, DW, 07/23/2018
# change the form of 'emiss' as raster, DW, 07/23/2018
# remove dmass correction
# remove zero footprint to save time and space, DW, 01/29/2019 

#########
ff.trajfoot <- function(trajdat, emiss){

  library(dplyr)

  # compute the total indx before any operations
  tot.p <- data.frame(indx = unique(trajdat$indx))

  # cut particles beyond emission grid
  # remove zero footprint to save time and space, DW, 01/29/2019 
  trajdat <- trajdat %>% filter(long >= extent(emiss)@xmin,
                                long <  extent(emiss)@xmax, 
                                lati >= extent(emiss)@ymin,
                                lati <  extent(emiss)@ymax, 
                                foot > 0)

  # find emissions and calculate CO2
  trajcor <- trajdat %>% dplyr::select(long, lati) 
  coordinates(trajcor) <- c('long', 'lati')
  trajdat <- trajdat %>% mutate(find.emiss = raster::extract(x = emiss, y = trajcor),
                                co2 = find.emiss * foot)

  # checking --
  if (F) {
    m1 <- ggplot.map(map = 'ggmap', center.lat = lon.lat$citylat + 0.5,
                     center.lon = lon.lat$citylon + 0.1, zoom = 7)[[1]]
    c1 <- m1 + geom_point(data = trajdat[trajdat$indx < 2000, ],
                          aes(long, lati, colour = co2)) +
               scale_colour_gradient(low = 'yellow', high = 'red', trans = 'log10')
    e1 <- m1 + geom_point(data = trajdat[trajdat$indx < 2000, ],
                          aes(long, lati, colour = find.emiss)) +
               scale_colour_gradient(low = 'yellow', high = 'red', trans = 'log10')
  }

  # also, compute total dCO2 for each traj over all backwards hours
  sum.trajdat <- trajdat %>% group_by(indx) %>% na.omit() %>%
                             dplyr::summarize(ff.sum = sum(co2)) %>% ungroup()

  # since we removed particles with zero footprint, we need to add zero to ff.sum 
  # to let `sum.trajdat` have the same amount of total particles as `tot.p`, 
  # DW, 01/29/2019 
  sum.trajdat <- sum.trajdat %>% right_join(tot.p, by = 'indx') 

  # NA will show up in above `sum.trajdat` for particles with foot = 0
  # thus, replace NA with 0 
  sum.trajdat$ff.sum[is.na(sum.trajdat$ff.sum)] <- 0 

  return(sum.trajdat)

} # end of subroutine
