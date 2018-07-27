# subroutine to search OCO-2 files for any overpases for a given region
# DW, 05/15/2017
# add one count for good quality data, DW, 12/20/2017


# target.region contains c(min.lat, max.lat, min.lon, max.lon), NEED THIS ORDER!!!
# default city center is for Riyadh
find.overpass <- function(date, targe.region, oco2.version = c("b7rb", "b8r")){

  library(geosphere)

  # path and filename for storing OCO-2 info
  oco2.path<-paste("/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_input/OCO-2/OCO2_lite_",oco2.version,"/",sep="")
  allfile<-list.files(pattern="oco2_LtCO2_",path=oco2.path)

  #start<-paste("oco2_LtCO2_",substr(date[1],1,8),"_B7000rb_COoffline.nc4",sep="")
  #end<-paste("oco2_LtCO2_",substr(date[length(date)],1,8),"_B7000rb_COoffline.nc4",sep="")

  # there can be missing files
  file.info<-matrix(unlist(strsplit(allfile,"_")),ncol=5,byrow=T)
  ocodate<-file.info[,3]

  if(oco2.version=="b8r"){
    ocofile<-allfile[ocodate >= substr(date[1],3,8) & ocodate <= substr(date[2],3,8)]
    timestr<-substr(ocofile,12,17)
    timestr<-paste("20",timestr,sep="")
  }else if(oco2.version=="b7rb"){
    ocodate[file.info[,4]=="B7305Br"]<-paste("20",ocodate[file.info[,4]=="B7305Br"],sep="")
    ocofile<-allfile[ocodate >= date[1] & ocodate <= date[2]]
    timestr<-substr(ocofile,12,19)
    timestr[substr(timestr,8,8)=="B"]<-paste("20",substr(timestr[substr(timestr,8,8)=="B"],1,6),sep="")
    #timestr<-sort(timestr)
  }

  # loop over each overpass
  result<-NULL
  for (f in 1:length(ocofile)){

    if(f%%20==0)cat(paste("#------ ",signif(f/length(ocofile),3)*100,"% SEARCHED ------#\n"))
    ocodat<-nc_open(paste(oco2.path,ocofile[f],sep=""))
    oco.lat<-ncvar_get(ocodat,"latitude");oco.lon<-ncvar_get(ocodat,"longitude")      # grabbing OCO-2 levels, lat, lon
    oco.qf<-ncvar_get(ocodat,"xco2_quality_flag")
    oco.wl<-ncvar_get(ocodat,"warn_level")

    # determine whether there are any overpasses for the target region
    flag<-oco.lat >= target.region[1] & oco.lat <= target.region[2] & oco.lon >= target.region[3] & oco.lon <= target.region[4]
    tot.count<-length(which(flag==T))
    qf.count<-length(which(flag==T & oco.qf ==0))
    if(oco2.version=="b8r")wl.count<-length(which(flag==T & oco.wl ==0))
    if(oco2.version=="b7rb")wl.count<-length(which(flag==T & oco.wl <= 15))

    if(tot.count>0){
      tmp<-cbind(as.numeric(timestr[f]),as.numeric(tot.count),as.numeric(qf.count),as.numeric(wl.count))
      result<-rbind(result,tmp)
    }
    nc_close(ocodat)
  }

  result<-as.data.frame(result)
  colnames(result)<-c("timestr","tot.count","qf.count","wl.count")
  result<-result[order(result$timestr),]

  return(result)
}
