# subroutine to weight the column of footprint for each particle
# create new profile based on OCO-2 pressure weighting and averaging kernel
# and then apply to STILT footprint, by weighting the footprints of particles based on different releasing heights
# OCO-2 only provides the PW, AK and a priori at 20 levels, we might wanna use linear interpolation to "approx" values at given STILT releasing levels

# add 1, for obtaining interpolated weighting function, call "get.weight.func()"
# fix 2, 05/22/2017, trajec with level column is not passed to weight.trajecfoot(), add now...DW

if(F){
	ident=ident
	trajdat=orig.trajdat
	each.par=100
	oco.ak.norm=oco2.profiles$ak.norm
	oco.pw=oco2.profiles$pw
	oco.pres=oco2.profiles$pres
	oco.apriori=oco2.profiles$apriori
	recp.grdhgt=recp.grdhgt
	new.trajpath="./"
	ak.weight=TRUE
	pw.weight=TRUE
	overwrite=FALSE
}
weight.trajecfoot<-function(ident=outname, trajdat=trajdat, each.par=each.par, oco.ak.norm=sel.ak.norm, oco.pw=sel.pw, oco.pres=sel.pres, oco.apriori=sel.apriori, recp.grdhgt=recp.grdhgt, new.trajpath=new.trajpath, ak.weight=TRUE, pw.weight=TRUE, overwrite=FALSE){

# HERE, ak.weight and pw.weight is passed on for weighting trajec
combine.profile<-get.weight.func(ident=ident, trajdat=trajdat, each.par=each.par,
																oco.ak.norm=oco.ak.norm, oco.pw=oco.pw, oco.pres=oco.pres, oco.apriori=oco.apriori,
																recp.grdhgt=recp.grdhgt, ak.weight=ak.weight, pw.weight=pw.weight)

### STARTing weighting trajec based on profiles
if(ak.weight==F & pw.weight==F){

		# if ak.weight==F && pw.weight==F, return trajec with original footprint, BE CAREFULLY!!!
		# no longer need any following weighting...
		# !!! but still need to return weighting functions and other info
		cat("weight.trajecfoot(): no weighting turn on, return profiles and original trajec\n")
		result<-list(combine.profile,trajdat)

}else{

	#### ---------------------- START WEIGHTING FOOTPRINT COLUMN FROM .RData FILE ---------------------- ####
	# now if particle has a level of one, adjust footprint column by multiplying footprint with AK*PW weighting
	# previously, particles have already been group by trajdat$level
	# first of all, initialize newfoot with original foot, find whether weighted traj existed
	find.flag<-TRUE
	if (overwrite==FALSE){
		find.traj<-list.files(path=new.trajpath, pattern=paste(ident,".RData",sep=""), all.files=TRUE)
		if(length(find.traj)==0)find.flag<-FALSE
		if(find.flag==TRUE){
			cat("weight.trajecfoot(): FOUND weighted traj..");cat("\n")
			result<-getr(ident,path=new.trajpath)	# get weighting profiles and weighted trajec
			#new.traj<-result[[2]]
			#combine.profile<-result[[1]]
		}	# end if
	}  # end overwrite ==FALSE

	# if overwrite == TRUE OR there is no file found, weight the trajec
	if(overwrite==TRUE | find.flag==FALSE){

		# group particles, sort traj files by "index", 05/22/2017
		# add one more column indicating which starting level that those particles belong to
		level<-trajdat[,"index"]%/% each.par
		trajdat<-cbind(trajdat,level)
		trajdat[trajdat[,"index"]%%each.par!=0,"level"]<-trajdat[trajdat[,"index"]%%each.par!=0,"level"]+1
		stilt.nlevel<-max(trajdat[,"level"])

		# initialize weighted foot column with normal footprint
		newfoot<-trajdat[,"foot"]
		trajdat<-cbind(trajdat,newfoot)

		#cat("weight.trajecfoot(): Weighted-traj NOT FOUND, weighting foot columns for each group of particles...");cat("\n")
		# weighting newfoot by multipling AK and PW profiles from "combine.profile", along with number of STILT levels
		stilt.profile<-combine.profile[combine.profile$stiltTF==TRUE,]

		# DW, 04/20/2017, add pw.weight flag too
		# only weight footprint in trajec if one of the two flags/or both are TRUE
		if(ak.weight==T & pw.weight==T){	# weighted by AK*PW and total STILT levels !!!
			cat("weight.trajecfoot(): weight trajec by both AK & PW profiles\n")
			for (l in 1:stilt.nlevel){trajdat[which(trajdat[,"level"]==l),"newfoot"]<-trajdat[which(trajdat[,"level"]==l),"foot"] * stilt.profile$ak.pw[l]*stilt.nlevel}

		}else if(ak.weight==F & pw.weight==T){				# weighted by PW and total STILT levels !!!
			cat("weight.trajecfoot(): weight trajec only by PW profiles\n")
			for (l in 1:stilt.nlevel){trajdat[which(trajdat[,"level"]==l),"newfoot"]<-trajdat[which(trajdat[,"level"]==l),"foot"] * stilt.profile$pw[l]*stilt.nlevel}

		}else if(ak.weight==T & pw.weight==F){					# weighted by only AK, for calculating dCO2 for each trajec in doing transport error !!!
			cat("weight.trajecfoot(): weight trajec only by AK profiles\n")
			for (l in 1:stilt.nlevel){trajdat[which(trajdat[,"level"]==l),"newfoot"]<-trajdat[which(trajdat[,"level"]==l),"foot"] * stilt.profile$ak.norm[l]}
		} # end if flag, ak.weight

		# for testing, store two sets of trajdat
		# one weighting over AK.norm*PW, newfoot are much smaller than original foot
		newtraj<-trajdat[,-which(colnames(trajdat)=="foot")]
		colnames(newtraj)[colnames(newtraj)=="newfoot"]<-"foot"

		# use assignr.r to store in new path,
		# DW -- 11/28/2016
		# add interpolated profiles in RData files as well, DW, 04/19/2017
		result<-list(combine.profile, newtraj)
		assignr(xname=ident, value=result, path=new.trajpath, printTF=TRUE)
	} # end if overwrite/find.flag, weighting foot column in traj

}	# end if for weighting traj

return(result)	# return both weighting profiles and weighted trajec

}  # end of subroutine
