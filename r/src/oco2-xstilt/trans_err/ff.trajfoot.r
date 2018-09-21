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
# change to read 'emiss' in the form of raster, DW, 07/23/2018
# remove dmass correction

#########
ff.trajfoot <- function(trajdat, emiss){

  library(dplyr)

  # cut particles beyond emission grid
  trajdat <- trajdat %>% filter(
    long >= extent(emiss)@xmin & long < extent(emiss)@xmax &
    lati >= extent(emiss)@ymin & lati < extent(emiss)@ymax)

  # find emissions and calculate CO2
  trajcor <- trajdat %>% dplyr::select(long, lati)
  trajdat <- trajdat %>%
    mutate(find.emiss = extract(x = emiss, y = trajcor), co2 = find.emiss * foot)

  # verify--
  if (F) {
    m1 <- ggplot.map(map = 'ggmap', center.lat = lon.lat$citylat,
      center.lon = lon.lat$citylon + 0.1, zoom = 8)[[1]]
    c1 <- m1 + geom_point(data = trajdat[trajdat$indx < 2000, ],
        aes(long, lati, colour = co2)) +
        scale_colour_gradient(low = 'yellow', high = 'red', trans = 'log10')
    e1 <- m1 + geom_point(data = trajdat[trajdat$indx < 2000, ],
        aes(long, lati, colour = find.emiss)) +
        scale_colour_gradient(low = 'yellow', high = 'red', trans = 'log10')
  }

  # also, compute total dCO2 for each traj over all backwards hours
  sum.trajdat <- trajdat %>% group_by(indx) %>%
    dplyr::summarize(ff.sum = sum(co2))

  return(sum.trajdat)

} # end of subroutine
