
# subroutinr to select trajec
# assign level index to original trajdat-- group particles, sort traj files by "index"
# add one more column indicating which starting level that those particles belong to

bootstrap.level<-function(orig.trajdat=NULL,orig.trajpath=orig.trajpath,orig.outname=orig.outname,orig.dpar=orig.dpar,orig.agl=orig.agl,sel.agl=sel.agl,assignTF=T){

  sel.maxagl<-max(sel.agl)
  sel.dh<-unique(diff(sel.agl))
  sel.numpar<-length(sel.agl)*orig.dpar
  xname<-paste(substr(orig.outname,1,32),"00000-",formatC(sel.maxagl,width=5,flag=0),"by",formatC(sel.dh,width=5,flag=0),"x",formatC(sel.numpar,width=5,flag=0),"P",sep="")

  # figure out the level index from base traj for the newly selected runs
  agl.index<-match(sel.agl, orig.agl)
  # for checking...
  #print(diff.agl<-orig.agl[agl.index]-sel.agl)

  # first try to find selected trajec, if not, readin orig.trajdat and select
  findTF<-list.files(pattern=xname,path=orig.trajpath)

  if(length(findTF)>0){

    sel.level.trajdat<-getr(xname=xname,path=orig.trajpath)

  }else{

    if(length(orig.trajdat)==0){
      cat("Reading base run traj...it takes time...\n")
      orig.outname<-substr(orig.outname,1, nchar(orig.outname)-6)
      orig.trajdat<-getr(path=orig.trajpath, xname=orig.outname)	# readin base run traj

      cat("bootstrap.level(): assigning LEVELS to original trajdat & matching selected LEVELS...takes a while\n")
      level<-orig.trajdat[,"index"]%/% orig.dpar
      orig.trajdat<-cbind(orig.trajdat,level)
      orig.trajdat[orig.trajdat[,"index"]%%orig.dpar!=0,"level"]<-orig.trajdat[orig.trajdat[,"index"]%%orig.dpar!=0,"level"]+1
    }

    # now subset the whole traj at multiple selected levels
    # find out whether trajdat's level ==agl.index
    cat("bootstrap.level(): subsetting according to selected levels...takes a while\n")
    level.flag<-orig.trajdat[,"level"] %in% agl.index	# if true, those particles have the same level
    sel.level.trajdat<-orig.trajdat[level.flag, ]

    if(assignTF)assignr(value=sel.level.trajdat,xname=xname,path=orig.trajpath)
  }

  return(list(sel.level.trajdat,agl.index))

}
