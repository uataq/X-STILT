# subroutine to readin ODIAC emissions, v2015a and then match with footprint column from STILT trajs
# use ak*pw weighted traj for footprint, find out the anthro emissions at the gridcell where trajdaticle falls into
# written by DIEN WU, 01/08/2017

# add ODIACv2016, 03/15/2017
# fix a bug, (all footprint columns are zero), 05/11/2017
# add dmassTF for weighting footprint if violate mass conservation, DW, 10/20/2017
# cut particles beyond emission grid, bug shows up for summertime track, DW, 12/04/2017

#########
# input from main script
debugTF<-FALSE
if(debugTF){

  trajdat=upper.traj
  odiac.dat=odiac.co2
  dmassTF=F
}

odiac.trajfoot<-function(trajdat=trajdat,odiac.dat=odiac.co2,dmassTF=F,max.foot=1E-6){

  library(reshape)

  # determine the grid cell that each traj falls into at each time step
  # hard to deal with 0.083333..., so divided by 1/120 to convert to integer
  # if for PRD, cut longitude by 60E and 130E
  odiac.lat<-as.numeric(rownames(odiac.co2))  #[lat,lon]
  odiac.lon<-as.numeric(colnames(odiac.co2))
  odiac.res<-odiac.lat[2]-odiac.lat[1]

  max.index<-max(trajdat[,"index"])
  full.index<-seq(1,max.index,1)

  # mass conservation correction, grabbed from Trajecfoot()
  if (dmassTF){

    #remove trajdaticles with too strong dmass violation
    ind<-unique(trajdat[trajdat[,"dmass"]>1E3|trajdat[,"dmass"]<1/1E3,"index"])

    if (length(ind) >= length(unique(trajdat[, "index"]))/2){
      message("Trajecvprm(): ", length(ind), ' of ', length(unique(trajdat[, "index"])), ' trajdaticles have mass defect; returning NA')
      return(NA)
    }
    trajdat<-trajdat[!trajdat[,"index"]%in%ind,]
    #trajdat<-cbind(trajdat,hour=trajdat[,"time"]/60)

    # get average dmass to "correct correction" (allow multiplication w/ dmass without changing total mass)
    # i.e. correction for average mass loss of trajdaticles, since they get attracted to areas of mass destruction
    mean.dmass<-rev(tapply(trajdat[, "dmass"], trajdat[, "time"], mean)) # this gives for each btime a mean dmass

    # DMM, To account for situations where mean.dmass is zero (mass violation total), need to avoid division by zero to avoid downstream problems.
    mean.dmass[which(mean.dmass == 0)] <- 0.00001

    # need to "merge" this with trajdat; can't use array since not all paticles are left at very large time
    nparleft<-rle(trajdat[, "time"])$length # number of times the same time is repeated
    mean.dmass<-rep(mean.dmass, nparleft) # long vector of mean dmass

    # need to link this info to each trajdaticles dmass: normalize individual dmass by mean.dmass
    trajdat[, "dmass"]<-trajdat[, "dmass"]/mean.dmass             # Dan Matross gets problems with that

    # also multiplied by dmass (accumulated weight of trajdaticles due to mass violation, normalized by average dmass to conserve total mass over time)
    trajdat[,"foot"]<-trajdat[,"foot"]*trajdat[,"dmass"]
  }

  # cut particles beyond emission grid [LAT, LON], bug shows up for summertime track, DW, 12/04/2017
  max.lon<-max(odiac.lon)+odiac.res;max.lat<-max(odiac.lat)+odiac.res
  min.lon<-min(odiac.lon);min.lat<-min(odiac.lat)

  # end particles (mostly higher particles, foot=0)
  trajdat<-trajdat[trajdat[,"lon"]< max.lon & trajdat[,"lon"]>= min.lon,]
  trajdat<-trajdat[trajdat[,"lat"]< max.lat & trajdat[,"lat"]>= min.lat,]

  #lat.index<-(trunc(trajdat[,"lat"]/odiac.res)-odiac.lat[1]/odiac.res)+1
  #lon.index<-(trunc(trajdat[,"lon"]/odiac.res)-odiac.lon[1]/odiac.res)+1
  lat.index<-findInterval(trajdat[,"lat"],odiac.lat)
  lon.index<-findInterval(trajdat[,"lon"],odiac.lon)
  combine.index<-paste(lat.index, lon.index)
  #diff.lat<-odiac.lat[lat.index]-trajdat[,"lat"]

  # select traj where foot is non-zero
  foot.flag<-trajdat[,"foot"]>max.foot
  sel.trajdat<-trajdat[foot.flag,]

  # initialize with NA for emissions
  emiss<-rep(NA,nrow(trajdat))
  trajdat<-cbind(trajdat,emiss)

  if(length(sel.trajdat)!=0){

    # conpute the lat lon indices, only for non-zero foot
    #sel.lat.index<-(trunc(sel.trajdat[,"lat"]/odiac.res)-odiac.lat[1]/odiac.res)+1
    #sel.lon.index<-(trunc(sel.trajdat[,"lon"]/odiac.res)-odiac.lon[1]/odiac.res)+1
    sel.lat.index<-findInterval(sel.trajdat[,"lat"],odiac.lat)
    sel.lon.index<-findInterval(sel.trajdat[,"lon"],odiac.lon)
    sel.combine.index<-paste(sel.lat.index, sel.lon.index)  #488541

    # find the unique gridcells with footprint > max.foot and
    uni.index<-unique(sel.combine.index)  # 166973
    uni.str<-matrix(as.numeric(unlist(strsplit(uni.index, " "))), ncol=2, byrow=TRUE)
    colnames(uni.str)<-c("lat","lon")

    # now find whether all combine.index are within the uni.index
    # if yes, return the number in "uni.index"; if no, return NA
    match.index<-match(combine.index, uni.index)

    # in order to save time, only grab CO2 from ODIAC ONLY ONCE for all uni.index
    uni.co2<-NULL
    cat("odiac.trajfoot(): Grabbing ODAIC emission. It takes time...\n")
    for(i in 1:length(uni.index))uni.co2<-c(uni.co2,odiac.co2[uni.str[i,"lat"],uni.str[i,"lon"]])

    # now put emission in this large matrix
    # if match.index==NA, no emission assigned; if not, find the emission from vector "uni.co2"
    assign.flag<-!is.na(match.index)  # assign only when match.index is non-NA
    nonzero.match.index<-match.index[assign.flag] # 1916649
    trajdat[assign.flag,"emiss"]<-uni.co2[nonzero.match.index]

  }

  # if emission is NA or all footprint == 0 , replace with zero
  trajdat[is.na(trajdat[,"emiss"]),"emiss"]<-0

  # NOW, foot and emiss have the same dim, multiple to get a dCO2
  dco2<-trajdat[,"emiss"]*trajdat[,"foot"]
  trajdat<-cbind(trajdat,dco2)

  # output trajdat in csv format for CESIUM PLOTTING
  writeTF<-FALSE
  if(writeTF){
    asl<-trajdat[,"agl"]+trajdat[,"grdht"]
    trajdat<-cbind(trajdat,asl)
    output<-trajdat[,c("time","index","lat","lon","asl","dco2")]

    filename<-"traj_dco2_riyadh.csv"
    headers = c("time","index","lat","lon","asl","dco2")
    write(headers, file = filename, ncolumns = length(headers), append = FALSE, sep = ",")

    result<-c(output[l,"time"],output[l,"index"],output[l,"lat"],output[l,"lon"],output[l,"asl"],output[l,"dco2"])
    write.table(result, file=filename, ncolumns = length(headers), append=TRUE,sep=",")
    # append TRUE, keep adding lines without overwriting

  }# end if writeTF

  # also, compute total dCO2 for each traj over all backwards hours
  traj.dco2.anthro<-tapply(trajdat[,"dco2"], trajdat[,"index"], sum)
  melt.dco2.anthro<-melt(traj.dco2.anthro);colnames(melt.dco2.anthro)<-c("index","dco2.anthro")

  # fill with all index, if dmassTF==T, particles being removed, fill with 0
  if(dmassTF){
    sel.index<-sort(unique(trajdat[,"index"]))
    rm.index<-full.index[full.index%in%sel.index==FALSE]
    rm.dco2.anthro<-data.frame(index=rm.index,dco2.anthro=rep(0,length(rm.index)))

    all.dco2.anthro<-rbind(melt.dco2.anthro,rm.dco2.anthro)
    sort.dco2.anthro<-all.dco2.anthro[order(all.dco2.anthro$index),]
  }else{
    sort.dco2.anthro<-melt.dco2.anthro[order(melt.dco2.anthro$index),]
  }

  return(sort.dco2.anthro)
} # end of subroutine
