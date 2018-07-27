assignr <- function(xname, value, path="", printTF=FALSE, gz=F){

#see also getr() and existsr()
#to avoid large databases that take too long to load memory (as R tries...)
#similar to assign(), but object is removed after call to assignr
#assigns object to name xname, saves it with save() under path/xname.RData
#and REMOVES it from local database (search()[1])
#example: assignr("test",temp,path="/mydir/") creates file "/mydir/test.RData"
#that can be attached and contains a single object with name "test"
#printTF can be set to TRUE to print a statement that the object has been assigned
#2/27/04 by CHG
#modified to allow compression of files (flag gz)


# change .RDataxname to xname.RData, no more hidden files, DW on 01/16/2017
# change paste() to file.path(), DW, 05/02/2018

#  $Id: assignr.r,v 1.4 2008-08-12 08:51:22 skoerner Exp $
#---------------------------------------------------------------------------------------------------

     assign(x=xname, value=value, pos=1)

     # change paste() to file.path(), DW, 05/02/2018
     rfile <- file.path(path, paste(xname, ".RData", sep=""))
     #save(list=xname, file=paste(path,xname,".RData",sep=""))
     save(list=xname, file=rfile)
     remove(list=xname, pos=1)

     if(gz){
       unix(paste("gzip ", rfile, sep=""))
       xname<-paste(xname,".gz", sep="")
     }

     if (printTF) cat(rfile, " created\n", sep="")

}
