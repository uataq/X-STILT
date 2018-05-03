getr <- function(xname, path="./", gz=FALSE){

#see also assignr() and existsr()
#to avoid large databases that take too long to load memory (as R tries...)
#similar to get()
#gets object from stored location (file path/xname.Rdata)
#name of object is xname if file was writen with assignr()
#2/27/04 by CHG
#modified to gunzip gzipped files (flag gz)

# change .RDataxname to xname.RData, no more hidden files, DW on 01/16/2017
# change paste() to file.path(), DW, 05/02/2018

#  $Id: getr.r,v 1.5 2008-03-26 19:19:17 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------

     xxname <- paste(xname, ".RData", sep="")

     # change paste() to file.path(), DW, 05/02/2018
     rfile <- file.path(path, xxname)

     if(gz){
       fls <- dir(path, pattern=xxname, all.files=TRUE)
       if(xxname %in% fls){
         print("found unzipped version, use this")
         delflag <- F
       }else{
         # unzip gz file
         unix(paste("gunzip -c ", rfile, ".gz > ", rfile, sep=""))
         delflag <- T
       }  # end if find file
     }  # end if gz

     #attach(paste(path,xname,".RData",sep=""),pos=2)
     attach(rfile, pos=2)
     dat <- get(xname, pos=2)
     detach(2)

     if(gz){
       if(delflag){ #delete unzipped file, wasn't there before, avoid build up ...
         unix(paste("rm ", rfile, sep=""))
       }
     }
     return(dat)
}
