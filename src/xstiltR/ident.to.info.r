# subroutine to get info from trajec name
# by DW, 07/10/2017

# ident can be a vector
ident.to.info<-function(ident=ident){

  if(grepl(".RData",ident[1]))ident<-substr(ident,1,nchar(ident)-6) # get rid of .RData
  library(stringr)
  bin<-str_count(ident,"x")+1

  all.recp<-NULL;all.agl<-NULL
  for (i in 1:length(ident)){

    info<-matrix(unlist(strsplit(ident[i],"x")), ncol=bin[i], byrow=TRUE)
    if(bin[i]==10)colnames(info)<-c("year","mon","day","hour","min","lat","lon","agl","numpar","dxyp")
    if(bin[i]==9)colnames(info)<-c("year","mon","day","hour","lat","lon","agl","numpar","dxyp")
    if(bin[i]==8)colnames(info)<-c("year","mon","day","hour","lat","lon","agl","numpar")
    if(bin[i]==7)colnames(info)<-c("year","mon","day","hour","lat","lon","agl")
    aglinfo<-as.character(info[,"agl"])
    recp.year<-as.numeric(info[,"year"]);recp.mon<-as.numeric(info[,"mon"])
    recp.day<-as.numeric(info[,"day"]);recp.hour<-as.numeric(info[,"hour"])
    recp.lat<-as.numeric(substr(info[,"lat"],1,nchar(info[,"lat"])-1))
    recp.lon<-as.numeric(substr(info[,"lon"],1,nchar(info[,"lon"])-1))
    timestr<-as.numeric(paste(info[,"year"],formatC(info[,"mon"],width=2, flag=0),formatC(info[,"day"],width=2, flag=0),formatC(info[,"hour"],width=2, flag=0),sep=""))

    # grab height info, column or fixed heights
    if(grepl("by",aglinfo)){
    	find.agl<-unlist(strsplit(aglinfo,"&"))
      columnTF<-TRUE

    	# lower AGL levels if unequal dh, dw (02/08/2017)
    	# for summer time, we have three different agl
    	nagl<-length(find.agl)
    	stilt.agl<-NULL
    	for(hh in 1:nagl){
    		max<-as.numeric(substr(find.agl[hh],7,11))
    		min<-as.numeric(substr(find.agl[hh],1,5))
    		dh <-as.numeric(substr(find.agl[hh],nchar(find.agl[hh])-4, nchar(find.agl[hh])))
    		agl<-seq(min, max, dh)
    		stilt.agl<-c(stilt.agl,agl)
    	}

    }else{		# for const AGL
      columnTF<-FALSE
    	#maxagl<-as.numeric(substr(aglinfo,7,11))
    	#minagl<-as.numeric(substr(aglinfo,1,5))
    	#dh<-as.numeric(substr(aglinfo,14,18))
    	stilt.agl<-as.numeric(aglinfo)
    }

    stilt.nlevel<-length(stilt.agl)
    stilt.npar<-as.numeric(gsub("P","",info[,"numpar"]))

    # storing results...
    each.par<-NA;dxyp<-NA;recp.min<-NA
    if(columnTF)each.par<-stilt.npar/stilt.nlevel	# for all trajs
    if(bin[i]>= 9)dxyp<-as.numeric(gsub("deg","",info[,"dxyp"]))
    if(bin[i]==10)recp.min<-as.numeric(info[,"min"])

    tmp<-data.frame(recp.year,recp.mon,recp.day,recp.hour,recp.min,recp.lat,recp.lon,dxyp,timestr,stilt.npar,each.par,stilt.nlevel,aglinfo)
    all.agl<-rbind(all.agl,stilt.agl)
    all.recp<-rbind(all.recp,tmp)
  }

  rownames(all.agl)<-seq(1,length(ident))
  rownames(all.recp)<-seq(1,length(ident))

  return(list(all.recp,all.agl))
}
