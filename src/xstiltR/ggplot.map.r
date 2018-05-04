# subroutine for loading map using ggplot2
# add ggmap as well, DW, 10/25/2017

ggplot.map<-function(map=NULL,minlat=NULL,maxlat=NULL,minlon=NULL,maxlon=NULL, center.lat=NULL, center.lon=NULL, zoom=8, shiftlat=0, shiftlon=0, ocean.col="lightsteelblue2",land.col="black",land.outline="gray30",maptype='roadmap'){

  if(map=="black"){
    # load map, from https://susanejohnston.wordpress.com/2012/07/03/creating-a-large-scale-map-using-ggplot2-a-step-by-step-guide/
    library(ggplot2);library(maptools);gpclibPermit();library(reshape)
    worldmap<-readShapeSpatial("/uufs/chpc.utah.edu/common/home/lin-group1/wde/TM_WORLD_BORDERS-0.3.shp")
    worldmap<-fortify(worldmap)

    # plot 2D map first
    latlimits<-c(minlat+shiftlat, maxlat-shiftlat)
    lonlimits<-c(minlon+shiftlon, maxlon-shiftlon)

    ticks<-c(0.1,0.2,0.5,1,2,5,10,20,50,100)
    seplat<-ticks[findInterval((maxlat-minlat)/5,ticks)]
    seplon<-ticks[findInterval((maxlon-minlon)/5,ticks)]

    latrange<-seq(minlat,maxlat,seplat)
    lonrange<-seq(minlon,maxlon,seplon)

    ylabels<-paste(latrange,"ºN",sep="")
    xlabels<-paste(lonrange,"ºE",sep="")
    if(length(unique(lonrange<0))==2){
      #lonrange[lonrange<0]<-abs(lonrange[lonrange<0])
      xlabels<-paste(lonrange,"ºE",sep="")
      xlabels[lonrange<0]<-paste(abs(lonrange[lonrange<0]),"ºW",sep="")
    }

    m1<-ggplot()+
        geom_polygon(data=worldmap,aes(x=long,y=lat,group=group),fill=land.col)+
        geom_path(data=worldmap,aes(x=long,y=lat,group=group),colour=land.outline)+
        coord_cartesian(xlim=lonlimits,ylim=latlimits)+
        scale_y_continuous(breaks=seq(minlat,maxlat,seplat),labels=ylabels)+
        scale_x_continuous(breaks=seq(minlon,maxlon,seplon),labels=xlabels)+
        theme(panel.background=element_rect(fill=ocean.col,colour="grey"),panel.grid.major=element_line(colour=NA),panel.grid.minor=element_blank(),
              axis.ticks=element_blank(),axis.text.x=element_text(size=12,vjust=0),axis.text.y=element_text(size=12,hjust=1.2))

    return(m1)
  }  # end if black map

  if(map=="ggmap"){

      library(ggmap)
      # load google map
      font.size=rel(1.3);alpha=1
      sitemap<-get_map(location=c(lon=center.lon,lat=center.lat),zoom=zoom,maptype=maptype,crop=FALSE)

      ### important to deal with the shift, grab center lat lon from ggmap
      box<-attr(sitemap,"bb")
      ggmap.lat<-(box$ll.lat+box$ur.lat)/2
      ggmap.lon<-(box$ll.lon+box$ur.lon)/2
      shift.lat<-ggmap.lat-center.lat
      shift.lon<-ggmap.lon-center.lon
      #maptype c("terrain", "terrain-background", "satellite", "roadmap", "hybrid", "toner", "watercolor", "terrain-labels", "terrain-lines", "toner-2010", "toner-2011", "toner-background", "toner-hybrid", "toner-labels", "toner-lines", "toner-lite")
      m1<-ggmap(sitemap)
      return(list(m1,shift.lat,shift.lon))
  } # end if ggmap

}
