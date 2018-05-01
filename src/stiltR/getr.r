getr<-function(xname, path="./",gz=FALSE){
#see also assignr() and existsr()
#to avoid large databases that take too long to load memory (as R tries...)
#similar to get()
#gets object from stored location (file path/xname.Rdata)
#name of object is xname if file was writen with assignr()
#2/27/04 by CHG
#modified to gunzip gzipped files (flag gz)
#edited by DW on 01/16/2017 (change .RDataxname to xname.RData, no more hidden files)

#  $Id: getr.r,v 1.5 2008-03-26 19:19:17 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------

     xxname<-paste(xname,".RData",sep="")
     if(gz){
       fls<-dir(path,pattern=xxname,all.files=TRUE)
       if(xxname%in%fls) {
         print("found unzipped version, use this")
         delflag<-F
       }else{
         unix(paste("gunzip -c ",path,xname,".RData.gz > ",path,xname,".RData",sep=""))
         delflag<-T
       }
     }
     attach(paste(path,xname,".RData",sep=""),pos=2)
     dat<-get(xname,pos=2)
     detach(2)
     if(gz){
       if(delflag){#delete unzipped file, wasn't there before, avoid build up ...
         unix(paste("rm ",path,xname,".RData",sep=""))
       }
     }
     return(dat)
}
