# subroutine for loading map using ggplot2
# add ggmap as well, DW, 10/25/2017

####### ####### ####### ####### Input variables ####### ####### ####### #######
# map: loading maps with black background and blue ocean, OR google map
#
# --> if map == 'black':
# must-have variables: 'minlat', 'maxlat', 'minlon', 'maxlon' (lat lon boundary)
# variables for color control: 'ocean.col', 'land.col', 'land.outline',
#
# --> if map == 'ggmap':
# must-have variables: 'maptype', 'center.lat', 'center.lon', 'zoom'
#
# maptype: c('terrain', 'terrain-background', 'satellite', 'roadmap',
#            'hybrid', 'toner', 'watercolor', 'terrain-labels', 'terrain-lines',
#            'toner-2010', 'toner-2011', 'toner-background', 'toner-hybrid',
#            'toner-labels', 'toner-lines', 'toner-lite')
#
####### ####### ####### ####### ####### ####### ####### ####### ####### #######

ggplot.map <- function(map = c('black', 'ggmap'), maptype = 'roadmap', 
                       shape.file = './data/TM_WORLD_BORDERS-0.3.shp', 
                       minlat, maxlat, minlon, maxlon, land.col = 'black',
                       land.outline = 'gray30', us.land.outline = 'gray30',
                       ocean.col = 'lightsteelblue2', ocean.outline = 'grey', 
                       center.lat, center.lon, zoom = 8, color = c('color', 'bw')[1], 
                       source = c('google', 'osm', 'stamen')[1]){

  ## if load map with black background for land and blue for ocean (default col)
  if (map == 'black') {

    # load map, from https://susanejohnston.wordpress.com/
    # 2012/07/03/creating-a-large-scale-map-using-ggplot2-a-step-by-step-guide/
    library(ggplot2); library(maptools); gpclibPermit(); library(maps)

    usstates <- map_data('state')
    worldmap <- readShapeSpatial(shape.file)
    worldmap <- fortify(worldmap)

    # plot 2D map first
    latlimits <- c(minlat, maxlat)
    lonlimits <- c(minlon, maxlon)

    ticks  <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 20, 50, 100)
    seplat <- ticks[findInterval((maxlat-minlat)/7, ticks)]
    seplon <- ticks[findInterval((maxlon-minlon)/7, ticks)]

    latrange <- seq(minlat, maxlat, seplat)
    lonrange <- seq(minlon, maxlon, seplon)

    ylabels <- paste0(latrange, 'ºN')
    xlabels <- paste0(lonrange, 'ºE')

    xlabels[lonrange < 0] <- paste0(abs(lonrange[lonrange < 0]), 'ºW')
    ylabels[latrange < 0] <- paste0(abs(latrange[latrange < 0]), 'ºS')

    m1 <- ggplot() + labs(x = 'LONGITUDE', y = 'LATITUDE') +
          geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = land.col) +
          geom_path(data = usstates, aes(x = long, y = lat, group = group), colour = us.land.outline) +
          geom_path(data = worldmap, aes(x = long, y = lat, group = group), colour = land.outline) +
          coord_cartesian(xlim = lonlimits, ylim = latlimits) +
          scale_y_continuous(breaks = latrange, labels = ylabels) +
          scale_x_continuous(breaks = lonrange, labels = xlabels) +
          theme(panel.background = element_rect(fill = ocean.col, colour = ocean.outline),
                panel.grid.major = element_line(colour = NA),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(size = 12, vjust = 0),
                axis.text.y = element_text(size = 12, hjust = 1.2))

    return(m1)
  }  # end if black map


  ## if loading google map, with centered lat/lon
  if (map == 'ggmap') {

    library(ggmap);library(ggplot2)

    # load google map
    font.size = rel(1.2)
    sitemap <- get_map(location = c(lon = center.lon, lat = center.lat),
                       zoom = zoom, maptype = maptype, color = color, 
                       source = source, crop = FALSE)

    ### important to deal with the shift, grab center lat lon from ggmap
    box <- attr(sitemap, 'bb')
    ggmap.lat <- (box$ll.lat + box$ur.lat) / 2
    ggmap.lon <- (box$ll.lon + box$ur.lon) / 2
    shift.lat <- ggmap.lat - center.lat
    shift.lon <- ggmap.lon - center.lon

    m1 <- ggmap(sitemap) 
    return(list(m1, shift.lat, shift.lon))
  } # end if ggmap

}# end of subroutine
