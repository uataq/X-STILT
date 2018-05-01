# script to define background region with forward-time run interpolated polygons
# by DW, 11/02/2017

# try box receptors/sources,"dxyp" in trajecmulti(), DW, 11/08/2017
# add time windows for releasing particles, DW, 11/09/2017
# clear things up and add more sites, DW, 04/16/2018

sourcepath<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/stiltR/"
source(paste(sourcepath,"sourceall.r",sep=""))
library(ggplot2);library(reshape)

####
index<-1  # for choosing sites
met<-c("1km","gdas","gdas0p5")[3] # met files
site<-c("Riyadh","Cairo","PRD","Lanzhou")[index]
oco.hr<-c("10","10","05","06")[index]

#columnTF<-F
# if using multiple receptors, or box of receptors or sources, turn it on,
# then call updated Trajecmulti() instead of Trajec()
multiTF<-T    # whether to release particles from a box or a column
windowTF<-T   # whether to continously release particles
nummodel<-15;if(multiTF)nummodel=997  # multiTF for releasing from a box or column
uncertTF<-T 	# sim using trajec with transport error components
#offset.hr<-6  # first releasing hours, # of hours ahead of OCO-2 overpassing time, e.g., 1000UTC overpass time, first release at 0400UTC.
offset.hr<-10   # used in paper for Riyadh

str<-"_forward"
if(multiTF)str<-"_box"
if(uncertTF==TRUE)str<-paste(str,"_uncert",sep="")
aglstr<-"fixed_agl"

# generate forward trajec or directly grab existing trajec
trajpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/Traj/",site,"/",met,"/",aglstr,"/orig_traj",str,"/",sep="")

# overpass time, YYYYMMDDHH
#sim.timestr<-c("2014122710","2014122910","2015012810","2015081710","2015121610","2016011510","2016021610","2016072510","2016103110") # sim for riyadh
#sim.timestr<-c("2015022810","2015031810","2016022411") # sim for cairo
#sim.timestr<-"2015011505"  # sim for PRD
sim.timestr<-"0"

find.track<-read.table(paste("/uufs/chpc.utah.edu/common/home/lin-group1/wde/OCO-2_proj/count_overpass_",site,"_b8r.txt",sep=""),sep=",",header=T)
if(site=="Riyadh"){
  # analysis on Ramandan signals, overpass dates fall during 06/17-07/17/2015, 06/06-07/05/2016, 05/26-06/24/2017
  ram.time15<-seq(as.Date("2015/06/17"),as.Date("2015/07/17"),by="day")
  ram.time16<-seq(as.Date("2016/06/06"),as.Date("2016/07/05"),by="day")
  ram.time17<-seq(as.Date("2017/05/26"),as.Date("2017/06/24"),by="day")
  ram.time<-c(ram.time15,ram.time16,ram.time17)
  find.track$date<-as.Date(as.character(find.track$timestr),"%Y%m%d")
  ram.info<-find.track[!is.na(match(find.track$date,ram.time)),]
  track.timestr<-paste(ram.info$timestr[ram.info$wl.count>0],oco.hr,sep="")

}

#find.track<-find.track[find.track$tot.count>=400 & find.track$wl.count>100,]
find.track<-find.track[find.track$tot.count>300 & find.track$wl.count>10,]
sel.timestr<-paste(find.track$timestr,oco.hr,sep="")
find.track$wod<-weekdays(as.Date(sel.timestr,"%Y%m%d"))
#track.timestr<-sel.timestr[-match(sim.timestr,sel.timestr)]
track.timestr<-sel.timestr

#nhrs<-as.numeric(substr(track.timestr,9,10))+2  # allow for 2 hours for upper bound
track.timestr<-as.character(as.numeric(track.timestr)-offset.hr)
filestr<-paste(substr(track.timestr,1,4),"x",substr(track.timestr,5,6),"x",substr(track.timestr,7,8),"x",substr(track.timestr,9,10),sep="")

### whether to generate trajectories
overwriteTF<-T
if(overwriteTF){
  rundir<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/STILT_Exe/"  # Where to run STILT, where Copy lies
  metpath<-"/uufs/chpc.utah.edu/common/home/u0947337/"    	# path for the ARL format of WRF and GDAS # Meteological paths
  #metpath<-"/uufs/chpc.utah.edu/common/home/lin-group1/group_data/ARL_model_data/GDAS0p5/"

  ### looping over all timestr
  for (i in 1:length(track.timestr)){

    timestr<-track.timestr[i] # current overpass time

    # get met files
    # e.g., ARL naming convection like this wrfout_d04_cairo_201503.arl
    # always remember to move the WRF files in home directory
    if(met=="1km")metfile<-rev(list.files(path=metpath, pattern=substr(timestr,1,6)))  # grab met files based on year and month
    if(met=="gdas")metfile<-paste("gdas_",substr(timestr,1,6),sep="")  # grab met files based on year and month
    if(met=="gdas0p5")metfile<-list.files(path=metpath,pattern=paste("gdas0p5_",substr(timestr,1,8),sep=""))
    #if(met=="gdas0p5")metfile<-list.files(path=metpath,pattern=substr(timestr,1,8))
    varstrajec<-c("time","index","lat","lon","agl","grdht","foot","sampt","dmass","zi","pres")  # output variables
    print(metfile)

    #### obtaining wind errors, transport error component
    siguverr<-NULL;TLuverr<-NULL;horcoruverr<-NULL;zcoruverr<-NULL  # SD in wind errors, correlation timescale, horizontal and vertical lengthscales
    #;sigzierr=NULL;TLzierr=NULL;horcorzierr=NULL

    #### whether add additional wind error component, Lin & Gerbig, 2005
    if(uncertTF){
      if(met=="gdas0p5"){errpath<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/err_analysis/gdas_interpolate/temporal/riyadh/gdas_0p5deg/";TLuverr<-1*60;zcoruverr<-600;horcoruverr<-40}
      if(met=="gdas"){errpath<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/err_analysis/gdas_interpolate/temporal/riyadh/gdas_1deg/";TLuverr<-(2.39+2.45)/2*60;zcoruverr<-c(770,625,613)[index];horcoruverr<-c(109,85,98)[index]}

      if(F){
        # add errors, mainly siguverr, create a subroutine to compute siguverr from model-data wind comparisons
        err.info<-get.SIGUVERR(site=site,timestr=timestr,gdaspath=errpath,nfTF=FALSE,forwardTF=F)
        if(length(err.info)!=2)next
        met.rad<-err.info[[1]]
        siguverr<-as.numeric(err.info[[2]][1])
        u.bias<-as.numeric(err.info[[2]][2])
        v.bias<-as.numeric(err.info[[2]][3])

        cat(paste("SIGUVERR:",signif(siguverr,3),"\nu.bias:",signif(u.bias,3),"\nv.bias:",signif(v.bias,3),"\n"))
      }else{
        siguverr<-1.8   # make a very conservative assumption about the wind error
        #siguverr<-2
      }
    } # end if uncertTF

    #### lat, lon, agl, nhrs for receptors/sources for Riyadh
    if(site=="Riyadh"){recp.lat<-24.71;recp.lon<-46.75}     # center of city Riyadh
    if(site=="Lanzhou"){recp.lat<-c(36.03,36.56);recp.lon<-c(103.80,104.21)}     # center of city Lanzhou and Baiyin

    for(l in 1:length(recp.lat)){
    #l=2
    cat(paste("Generating forward trajec for recp.lat=",recp.lat[l]),"N\n")
    tmp.lat<-recp.lat[l];tmp.lon<-recp.lon[l]
    agl<-10;nlevel<-1
    recp.yr <-as.numeric(substr(timestr,1,4));recp.mon<-as.numeric(substr(timestr,5,6))
    recp.day<-as.numeric(substr(timestr,7,8));recp.hr <-as.numeric(substr(timestr,9,10))

    #### try box receptors/sources, DW, 11/08/2017
    # fortran code from Thomas Nehrkorn
    if(multiTF){

      # try automatically placing receptors
      dxyp<-0.2*2   # degree of lat, lon, thus latlon range will be recp.lat+/-dxyp, recp.lon+/-dxyp
      npar<-10000   # increase particle number, since we have box of receptors
      mn<-0   # minutes
      delt<-2   # record time step
      ident<-paste(recp.yr,"x",formatC(recp.mon,width=2,flag=0),"x",formatC(recp.day,width=2,flag=0),"x",formatC(recp.hr,width=2,flag=0),"x",
               formatC(signif(tmp.lat,6),width=5,format='f',digits=4,flag=0),"Nx",formatC(signif(tmp.lon,7),width=5,format='f',digits=4,flag=0),"Ex",
               formatC(agl,width=5,flag=0),"x",formatC(npar,width=5,flag=0),"Px",dxyp,"deg",sep="")

      #### whether to continously release particles
      if(windowTF){
        npar<-1000  # # of particles per run
        dtime<-seq(0,offset.hr,0.5)   # offset hours, e.g., release every half an hour
        dhr<-floor(dtime)       # in hours
        dmin<-dtime*60-dhr*60   # in minutes

        # update hours or days
        update.hr<-recp.hr+dhr  # same time zone as recp.hr, UTC
        update.day<-rep(recp.day,length(update.hr))

        forw.index<-update.hr>=24
        update.hr [forw.index]<-update.hr [forw.index]-24
        update.day[forw.index]<-update.day[forw.index]+1 # for forward run

        back.index<-update.hr<0
        update.hr [back.index]<-update.hr [back.index]+24
        update.day[back.index]<-update.day[back.index]-1 # for backward run

        #### update all receptor info, vectors for all input
        recp.yr<-rep(recp.yr,length(dtime));recp.mon<-rep(recp.mon,length(dtime))
        recp.day<-update.day;recp.hr<-update.hr;recp.min<-dmin   # mn in trajectmulti() has unit in minutes, simply is the receptor minutes...
        tmp.lat<-rep(tmp.lat,length(dtime));tmp.lon<-rep(tmp.lon,length(dtime));agl<-rep(agl,length(dtime))

        #### trajec names
        ident<-paste(recp.yr,"x",formatC(recp.mon,width=2,flag=0),"x",formatC(recp.day,width=2,flag=0),"x",formatC(recp.hr,width=2,flag=0),"x",formatC(recp.min,width=2,flag=0),"x",
                    formatC(signif(tmp.lat,6),width=5,format='f',digits=4,flag=0),"Nx",formatC(signif(tmp.lon,7),width=5,format='f',digits=4,flag=0),"Ex",
                    formatC(agl,width=5,flag=0),"x",formatC(npar,width=5,flag=0),"Px",dxyp,"deg",sep="")

      }  # endif windowTF

      #### the updated Trajecmulti() will randomly place receptors according to dxyp
      # call trajecmulti() to generate trajec, it takes time...
      nhrs<-as.numeric(oco.hr)+2
      info<-Trajecmulti(yr=recp.yr-2000,mon=recp.mon,day=recp.day,hr=recp.hr,mn=recp.min,dxyp=dxyp,dzp=0,
                        lat=tmp.lat,lon=tmp.lon,agl=agl,outname=ident,nhrs=rep(nhrs,length(recp.yr)),numpar=npar,nummodel=nummodel,
                        metd=c("fnl","awrf"),metfile=metfile,metlib=metpath,conv=F,doublefiles=T,overwrite=TRUE,
                        outpath=trajpath,varsout=varstrajec,rundir=rundir,setup.list=list(DELT=delt,VEGHT=0.5),
                        siguverr=siguverr,TLuverr=TLuverr, zcoruverr=zcoruverr, horcoruverr=horcoruverr)

      #info<-ident.to.info(ident=ident)[[1]]
      # check trajec
      #trajdat<-getr(xname=ident,path=trajpath)
      #range(trajdat[,"time"])
      #sel.traj<-trajdat[trajdat[,"time"]==2,]
    }
  }

    }else{          # multiTF==FALSE

      ### particles released from one fixed location
      npar<-2000
      ident<-paste(recp.yr[1],"x",formatC(recp.mon[1],width=2,flag=0),"x",formatC(recp.day[1],width=2,flag=0),"x",formatC(recp.hr[1],width=2,flag=0),"x",
               formatC(signif(recp.lat[l],6),width=5,format='f',digits=4,flag=0),"Nx",formatC(signif(recp.lon[l],7),width=5,format='f',digits=4,flag=0),"Ex",
               formatC(agl,width=5,flag=0),"x",formatC(npar,width=5,flag=0),sep="")

      info<-Trajec(yr=recp.yr[t]-2000,mon=recp.mon[t],day=recp.day[t],hr=recp.hr[t],
                  lat=recp.lat[l],lon=recp.lon[t],agl=agl,outname=ident[t],nhrs=nhrs,numpar=npar,nummodel=nummodel,
                  delt=2,metd=c("fnl","awrf"),mgmin=2000,veght=0.5,metfile=metfile,metlib=metpath,conv=F,doublefiles=T,
                  overwrite=TRUE,outpath=trajpath,varsout=varstrajec,rundir=rundir,
                  siguverr=siguverr,TLuverr=TLuverr, zcoruverr=zcoruverr, horcoruverr=horcoruverr)

    }  # end if multiTF
    #} # end for looping over all traj, loop l
  } # end for looping over all timestr, loop i
}   # end if overwriteTF

### for plotting particle distributions
plotTF<-F
if(plotTF){

  index<-1
  met<-c("gdas","gdas0p5")[2]
  site<-c("Riyadh","Cairo","PRD","Lanzhou")[index]
  oco.hr<-c(10,10,5,6)[index]
  uncertTF<-T
  multiTF<-T
  str<-"_forward";if(multiTF)str<-"_box"
  if(uncertTF==TRUE)str<-paste(str,"_uncert",sep="")

  # generate forward trajec or directly grab existing trajec
  trajpath<-paste("/uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/Traj/",site,"/",met,"/fixed_agl/orig_traj",str,"/",sep="")

  trajfile<-list.files(path=trajpath,pattern=".RData")
  if(uncertTF==FALSE)trajfile<-list.files(path=trajpath,pattern="46.7500E")
  trajinfo<-ident.to.info(ident=trajfile)[[1]]
  uni.lat<-unique(trajinfo$recp.lat)
  uni.lon<-unique(trajinfo$recp.lon)
  filestr<-unique(substr(trajfile,1,10))
  track.timestr<-paste(substr(filestr,1,4),substr(filestr,6,7),substr(filestr,9,10),sep="")

  # output polluted lat ranges in txt file
  #header<-c("timestr","north.bg","south.bg","final.bg","min.lat","max.lat")
  #write(header,file=paste("bg_lat_range_",site,".txt",sep=""),ncolumns=length(header),append=F,sep=",")
  tot.plume<-NULL

  ### looping over all timestr
  for (i in 1:length(track.timestr)){

    if(i==3 & site=="Riyadh")next
    timestr<-track.timestr[i]
    ident<-trajfile[substr(trajfile,1,10)==filestr[i]]
    #ident<-ident[7:length(ident)]

    info<-ident.to.info(ident=ident)[[1]]
    if(F){s=1;ident<-ident[info$recp.lat==uni.lat[s]]}
    ident<-substr(ident, 1, nchar(ident)-6)
    dxyp<-info$dxyp/2

    #### read and merge all particles, when we release particles continously
    merge.trajdat<-NULL
    for(ii in 1:length(ident)){
      trajdat<-getr(xname=ident[ii],path=trajpath)
      merge.trajdat<-rbind(merge.trajdat,trajdat)
    }

    # select particles, cut by 11 UTC, recptor hour is 4UTC
    # no need to minus 1, recp.hr, xtime=0, recp.hr+1, xtime=60min
    xtime<-(oco.hr-info$recp.hour[1]+2)*60 # in minutes
    sel.trajdat<-merge.trajdat[merge.trajdat[,"time"]<= xtime,]

    #### prepare for plotting
    zoom<-8#;  stat<-c("kernel_density","stat_density")[1] # 2D kernel density (5% of the max) for define city plume
    sourcepath<-"/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/stiltR/"
    source(paste(sourcepath,"sourceall.r",sep=""))

    #### call ggplot.forward.polygon() to plot particle distribution with kernel density and OCO-2 observed XCO2
    pp<-ggplot.forward.polygon(zoom=zoom,ident=ident,trajpath=trajpath,trajdat=sel.trajdat,multiTF=T,site=site,timestr=timestr,map="ggmap",oco2.version="b8r",pch.size=zoom/11,alpha=0.4,hr=1)

    #library(ggpubr)
    #mm7<-ggarrange(plotlist=list(pp1,pp2,pp3), labels = paste(seq(1,6,1),")",sep=""),ncol=2,nrow=2)
    #ggsave(mm7,filename=paste(stat,"_bg_merge_",timestr,"_",dxyp[1],"deg_err_window.png",sep=""),width=20,height=20)

    picname<-paste(stat,"_bg_",site,"_",timestr,"_",dxyp[1],"deg_err_window.png",sep="")
    #picname<-paste(stat,"_bg_Baiyin_",timestr,"_",dxyp[1],"deg_err_window.png",sep="")
    if(uncertTF==FALSE)picname<-paste(stat,"_bg_",site,"_",timestr,"_",dxyp[1],"deg_window.png",sep="")
    ggsave(pp[[1]],filename=picname,width=12,height=11)

    #pp.uncert<-pp[[1]]
    #pp.all<-pp.uncert+geom_polygon(data=bound.traj,aes(x=X,y=Y),colour="gray30",linetype=2,fill=NA,size=0.7,alpha=0.7)#+geom_polygon(data=bg.bound,aes(x=X,y=Y),colour="gray50",fill="gray40",alpha=0.5)
    #ggsave(pp.all,filename=picname,width=12,height=11) # save plot
    #ggsave(pp[[1]],filename=picname,width=12,height=11) # save plot

    if(length(pp)>1){
      ggsave(pp[[1]],filename=picname,width=12,height=11)

      bg.lat<-as.numeric(pp[[2]])
      result<-c(timestr,bg.lat)
      #write(result,file=paste("bg_lat_range_",site,".txt",sep=""),ncolumns=length(result),append=T,sep=",")

      ## also get city plume
      tmp.plume<-pp[[3]]
      tmp.plume$timestr<-rep(timestr,nrow(tmp.plume))
      tmp.plume$fac<-rep(i,nrow(tmp.plume))
      tot.plume<-rbind(tot.plume,tmp.plume)
    }

    cat("# ----- ",signif(i/length(track.timestr))*100,"% done ----- #\n\n")
    #gc()
  } # end for i

  #store boundary
  #assignr(xname=paste("city_plume_",site,sep=""),value=tot.plume,path="./")

  # plot latitude range
  library(ggplot2);font.size=rel(1.1);site<-"Riyadh"
  bg.lat<-read.table(paste("bg_lat_range_",site,".txt",sep=""),sep=",",header=T)
  bg.lat<-bg.lat[bg.lat$timestr!="20160803",]
  l1<-ggplot()+theme_bw()+geom_segment(data=bg.lat,aes(x=as.character(timestr),xend=as.character(timestr), y=min.lat, yend=max.lat),size=1.2)
  l2<-l1+labs(x="Date [YYYYMMDD]",y="Latitude [degN]",title=paste("Latitude range for city plume vs. OCO-2 tracks for",site))+theme(legend.position="bottom",legend.key.width=unit(0.5, "cm"),legend.key.height=unit(3, "cm"),axis.text.x = element_text(angle = 45, hjust = 1),
              legend.text=element_text(size=font.size),legend.key = element_blank(),panel.grid.minor=element_blank(),axis.title.y=element_text(size=font.size,angle = 90),axis.title.x=element_text(size=font.size,angle=0),
              axis.text=element_text(size=font.size),axis.ticks=element_line(size=font.size),title=element_text(size=font.size))
  #ggsave(l2,filename=paste("lat_range_timestr_",site,".png",sep=""),width=10,height=5)

} # enfif plotTF
