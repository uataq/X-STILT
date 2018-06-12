# script to find OCO-2 files given wanted lat, lon ranges
# DW, 05/15/2017

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- "/uufs/chpc.utah.edu/common/home"
workdir <- file.path(homedir, "lin-group1/wde/github/XSTILT")
source(file.path(workdir, "src/sourceall.r"))

library(ggplot2);library(ggrepel);library(scales);library(ncdf4)
library(ggmap);library(RColorBrewer);library(ggpubr);  library(geosphere)
library(cowplot)

####
index <- 14
site <- c("Riyadh", "Medina", "Mecca", "Jerusalem", "Cairo",
          "PRD", "Beijing", "Xian", "Lanzhou", "Mumbai", "Java",
          "Indianapolis", "Phoenix", "SLC")[index]

### spatial domain
# lat.lon <- c(lat.lon[1], lat.lon[2], lat.lon[3], lat.lon[4], city.lat, city.lon)
if(site=="Riyadh")lat.lon <- c(23.5, 26, 46, 48, 24.63+0.1, 46.72+0.1) # Riyadh
if(site=="Medina")lat.lon <- c(23.5, 25.5, 38, 41, 24.46, 39.60)  # Medina, 24.46N, 39.60E
if(site=="Mecca")lat.lon <- c(20.5, 22.5, 38.5, 40.5, 21.42, 39.82+0.2)  # Mecca, 21.42N, 39.82E
if(site=="Jerusalem")lat.lon <- c(31, 33, 34, 36, 31.78, 35.22)  # Jerusalem, 31.78N, 35.22E
if(site=="Cairo")lat.lon <- c(29, 32, 30, 32, 30.05+0.2, 31.23+0.1)  # Cairo

if(site=="PRD")lat.lon <- c(22, 24, 112.5, 115, 22.40+0.4, 114.11-0.4) # PRD
if(site=="Beijing")lat.lon <- c(39, 41, 115.5, 117.5, 39.90, 116.41+0.3) # Beijing, 39.9N, 116.41E
if(site=="Xian")lat.lon <- c(33, 35.5, 107.5, 110.5, 34.27, 108.90)  # Xi'an, 34.27N, 108.9E
if(site=="Lanzhou")lat.lon <- c(35, 37.5, 102.5, 105, 36.037+0.2, 103.8) # Lanzhou, 36.037N, 103.8E
if(site=="Mumbai")lat.lon <- c(17.5, 20, 72, 74.5, 18.98, 72.83+0.4) # Mumbai, 18.98N, 72.83E

if(site=="Indianapolis")lat.lon <- c(39, 41, -87, -85.5, 39.77, -86.15) # Indianapolis
if(site=="Phoenix")lat.lon <- c(32, 35, -113, -110.5, 33.45, -112.07+0.5) # Phoenix, 33.45N, 112.07W
if(site=="SLC")lat.lon <- c(37, 43, -114, -111, 40.75, -111.88) # SLC, 40.75N, 111.88W
target.region <- lat.lon[1:4]

# grabbing OCO-2 info
searchTF <- T
oco2.version <- c("b7rb", "b8r")[2]
ocopath <- paste(homedir, "/lin-group1/wde/STILT_input/OCO-2/OCO2_lite_",
                 oco2.version, "/", sep="")
filenm <- paste("./count_overpass_", site, "_", oco2.version, ".txt", sep="")
if(searchTF){
  date <- c(20140901, 20180228)
  find.info <- find.overpass(date=date, targe.region=target.region,
                             oco2.version=oco2.version)
  write.table(find.info, file=filenm, quote=F, row.names=F, sep=",")
}else{
  find.info <- read.table(file=filenm, header=T, sep=",")
}

## if no min distance, calculate min distance from the city center
if(length(grep("dist", colnames(find.info)))==0){
  find.info$dist <- rep(NA, nrow(find.info))

  for(d in 1:length(find.info$timestr)){
    obs <- grab.oco2(ocopath=ocopath, timestr=find.info$timestr[d], lat.lon=lat.lon)
    uni.foot <- unique(obs$foot)
    center.obs <- obs[obs$foot == uni.foot[1], c("lat", "lon")]

    dist <- NULL
    for(l in 1:nrow(center.obs)){
      dist <- c(dist, distCosine(c(lat.lon[6], lat.lon[5]),
                                 c(center.obs$lon[l], center.obs$lat[l])))
    }
    find.info$dist[d] <- min(dist)
  }
  write.table(find.info,file=filenm, quote=F, row.names=F, sep=",")
}

#### plot overpass info
gg <- 100  # number of screened soundings criteria
dd <- 50   # min distance in km
tt <- 500  # total number of soundings

val <- c("1"=cols[1], "2"=cols[2], "3"=cols[3], "4"=cols[4], "5"=cols[5])
cols<-c("gray70", "gray25", "orange", "brown", "deepskyblue2")

find.info <- find.info[find.info$tot.count > tt, ]
find.info$col <- cols[1]
find.info$x <- seq(1,nrow(find.info),1)
find.info[find.info$qf.count > gg, "col"] <- cols[2]
find.info[find.info$dist/1E3 < dd & find.info$qf.count > gg, "col"] <- cols[3]
lab <- c("# of total OBS",
         "# of screened OBS (QF=0)",
         "Min distance\nfrom city center [km]",
         "Mean RMSE [m/s]\nfor wind errors below 3km")

font.size=rel(1.0)
o1 <- ggplot()+theme_bw()
o1 <- o1+geom_bar(data=find.info, aes(x=x, y=tot.count, fill=as.factor(1)),
                  width=0.95, stat = "identity")
o1 <- o1+scale_x_continuous(name="Overpass Date [YYYY-MM-DD]", breaks=find.info$x,
                            labels=paste(substr(find.info$timestr,1,4), "-",
                                        substr(find.info$timestr,5,6), "-",
                                        substr(find.info$timestr,7,8), sep=""))
o1 <- o1+scale_y_continuous(breaks=seq(100,1000,100), labels=seq(100,1000,100),
                            name="Sounding numbers",
                            sec.axis=sec_axis(~., name="Min Distance from city center [km]"))
o2 <- o1+theme(legend.position="bottom",legend.key.width=unit(4, "cm"),
               legend.key.height=unit(0.3, "cm"),
               legend.text=element_text(size=font.size),
               legend.key = element_blank(),
               panel.grid.minor=element_blank(),
               axis.title.y=element_text(size=font.size,angle = 90),
               axis.title.y.right=element_text(color=cols[3]),
               axis.title.x=element_text(size=font.size,angle=0),
               axis.text=element_text(size=font.size),
               axis.ticks=element_line(size=font.size),title=element_text(size=font.size),
               axis.text.x = element_text(angle=45, hjust = 1,color=find.info$col),
               axis.text.y.right = element_text(color=cols[3]))

o2 <- o2+geom_bar(data=find.info, aes(x=x, y=qf.count, fill=as.factor(2)),
                  width=0.95, stat = "identity")
o2 <- o2+geom_hline(yintercept=dd, linetype=2, color=cols[3], size=1.1)+
         geom_hline(yintercept=gg, linetype=2, color=cols[2], size=1.3)

o3 <- o2+geom_bar(data=find.info[find.info$qf.count > gg,],
                  aes(x=x, y=dist/1E3, fill=as.factor(3)),
                  width=0.5, stat = "identity")

o4 <- o3+scale_fill_manual(name=NULL,values=val[1:3],labels=lab[1:3])+
         scale_colour_manual(name=NULL,values=val[3],labels=lab[3])

o4 <- o4+labs(title=paste("Overpass for", site, lat.lon[1], "-", lat.lon[2],
                           "N;  ", lat.lon[3],"-", lat.lon[4],"E"))
picname <- paste("tot_screened_dist_overpass_", site,".png",sep="")
ggsave(o4,filename=picname, width=20, height=8)


### plot overpass
filenm <- paste("./count_overpass_", site, "_", oco2.version, ".txt", sep="")
find.info <- read.table(file=filenm, header=T, sep=",")
#find.info <- find.info[find.info$tot.count> 1200, ]

#ram.time15 <- seq(as.Date("2015/05/17"), as.Date("2015/08/17"), by="day")  # 06/17-07/17/2015
#ram.time16 <- seq(as.Date("2016/05/06"), as.Date("2016/08/05"), by="day")
#ram.time17 <- seq(as.Date("2017/04/26"), as.Date("2017/07/24"), by="day")  # allow for more month
#ram.time <- c(ram.time15, ram.time16, ram.time17)
#find.info$date <- as.Date(as.character(find.info$timestr), "%Y%m%d")
#ram.info <- find.info[!is.na(match(find.info$date, ram.time)),]
#non.ram.info <- find.info[is.na(match(find.info$date, ram.time)),]
#print(ram.info)

alpha <- 1
tot <- nrow(find.info)

#font.size=rel(0.9)
#limit<-c(400,410)
col <- c('black','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8',
         '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')[-1]
m1 <- ggmap(get_map(location=c(lon=lat.lon[6],lat=lat.lon[5]),zoom=7,maptype='roadmap'))
z1 <- ggmap(get_map(location=c(lon=lat.lon[6],lat=lat.lon[5]-0.1),zoom=9,maptype='roadmap'))

picpath <- file.path(homedir, "u0947337/public_html/LAIR_group/OCO-2/overpass")

for (i in 1:nrow(find.info)){

  obs <- grab.oco2(ocopath=ocopath, timestr=find.info$timestr[i], lat.lon=lat.lon)
  urban.obs <- obs[obs$lat >= 40 & obs$lat <= 41 & obs$lon > -112.5 & obs$lon < -111.5, ]
  if(nrow(urban.obs)< 50){
    cat(find.info$timestr[i])
    next
  }
  print(range(urban.obs$xco2))

  limit <- c(388,406)
  if(find.info$timestr[i] >= 20150101)limit <- c(388,410) # 20 ppm
  if(find.info$timestr[i] >= 20160101)limit <- c(392,415) # 25 ppm
  if(find.info$timestr[i] >= 20170101)limit <- c(393,420) # 25 ppm
  if(find.info$timestr[i] >= 20180101)limit <- c(399,421) # 22 ppm

  font.size=rel(1.2)
  title <- paste("Map of XCO2 concentrations of OCO-2 overpasses ",
                 find.info$timestr[i],"\nover Utah [",lat.lon[1], "ºN-",
                 lat.lon[2],"ºN;",lat.lon[3],"ºE-",lat.lon[4],"ºE]", sep="")
  m2 <- m1+geom_point(data=obs,aes(x=lon,y=lat,colour=xco2),alpha=0.7)+
           labs(title=title,x="LONGITUDE [ºW]",y="LATITUDE [ºN]")+theme_bw()
  m2 <- m2+theme(legend.position="right",legend.text=element_text(size=font.size),
                 legend.key = element_blank(), axis.title.y=element_text(size=font.size,angle = 90),
                 axis.title.x=element_text(size=font.size,angle=0),axis.text=element_text(size=font.size),
                 axis.ticks=element_line(size=font.size),title=element_text(size=font.size))
  m3 <- m2+scale_colour_gradientn(name="XCO2:",colours=col,breaks=seq(380,430,2), limit=limit)+
           theme(legend.key.width=unit(0.5, "cm"),legend.key.height=unit(4, "cm"))
  #ggsave(m3,filename=paste("ggmap_xco2_",find.info$timestr[i],"_",site,"_",oco2.version,".png",sep=""),width=12,height=10)
  #map.all[[i]]<-p3

  # zoomed map
  title <- paste("Map of XCO2 concentrations of OCO-2 overpasses ", find.info$timestr[i],
                 "\nover ",site," [40ºN-41ºN;112.5ºW-112ºW]", sep="")
  z2 <- z1+geom_point(data=urban.obs,aes(x=lon,y=lat,colour=xco2),alpha=alpha, size=2)+
           labs(title=title,x="LONGITUDE [ºW]",y="LATITUDE [ºN]")+theme_bw()
  z2 <- z2+theme(legend.position="right",legend.text=element_text(size=font.size),
                 legend.key = element_blank(), axis.title.y=element_text(size=font.size,angle = 90),
                 axis.title.x=element_text(size=font.size,angle=0),axis.text=element_text(size=font.size),
                 axis.ticks=element_line(size=font.size),title=element_text(size=font.size))
  z3 <- z2+scale_colour_gradientn(name="XCO2:",colours=col,breaks=seq(380,430,2), limit=limit)+
           theme(legend.key.width=unit(0.5, "cm"),legend.key.height=unit(4, "cm"))

  ### plot on Latitude series
  title <- paste("Latitude series of XCO2 concentrations of OCO-2 overpasses ",
                find.info$timestr[i]," over ",site, sep="")
  p1 <- ggplot()+theme_bw()+geom_point(data=obs,aes(x=lat,y=xco2),colour="gray",shape=17,size=4)
  p1 <- p1+geom_point(data=obs[obs$qf==0,],aes(x=lat,y=xco2),shape=17,size=4)+
           labs(title=title,x="LATITUDE",y="XCO2 [ppm]")
  p2 <- p1+theme(legend.position="right",legend.key.width=unit(0.5, "cm"),
                 legend.key.height=unit(4, "cm"),legend.text=element_text(size=font.size),
                 legend.key = element_blank(),axis.title.y=element_text(size=font.size,angle = 90),
                 axis.title.x=element_text(size=font.size,angle=0),axis.text=element_text(size=font.size),
                 axis.ticks=element_line(size=font.size),title=element_text(size=font.size))
  p3 <- p2+scale_y_continuous(limit=limit,breaks=seq(380,430,2),labels=seq(380,430,2))
  #  highlight urban enhancements in red
  p4 <- p3+geom_point(data=urban.obs[urban.obs$qf==0,],aes(x=lat,y=xco2),shape=17,size=3,colour="orange")

  #ggsave(p3,filename=paste("LS_xco2_",find.info$timestr[i],"_",site,"_",oco2.version,".png",sep=""),width=12,height=10)
  #ls.all[[i]]<-p3

  mm <- ggdraw()+draw_plot(m3, x=0, y=.5, width=.5, height=.5) +
                draw_plot(z3, x=.5, y=.5, width=.5, height=.5) +
                draw_plot(p4, x=0, y=0, width =1, height=0.5) +
                draw_plot_label(label = c("A","B","C"), size = 15,
                                x = c(0,0.5,0), y = c(1,1,0.5))

  #mm <- ggarrange(plotlist=list(m3,z3,p3),labels=c("A","B", "C"),ncol=2,nrow=2,heights=c(1.7,1))  #,font.label=list(size = 20, face = "bold")
  picname <- paste("merge_xco2_",find.info$timestr[i],"_",site,"_",oco2.version,".png",sep="")
  picname <- file.path(picpath, picname)
  ggsave(mm, filename=picname, width=20, height=20)

}
