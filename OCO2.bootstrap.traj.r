# subroutine for bootstrapping the STILT trajectory
# it aims for a sensitivity test, trying to evaluate the impact of STILT column receptors on simulated XCO2
# we consider three factors -- 1) the max agl heights receptors go; 2) vertical spacing dh; and 3) the particles at each level

# the base run has Xreceptors set from 0m to 10,000m, with a dh=25m, and 200 particles released at each level
# thus, total 401 levels, total numpar = 401 * 200 particles

# Inputs--trajdat for base run, selected levels, selected particles
#		for selected levels, dh should be 25N, in order to use the 25m base run
# Outputs--create new names (with new agl, numpar), store into .RData file, and return the name of the .RData
# Written by Dien Wu, 10/31/2016

# move selecting levels to main script, 10/23/2017,DW

debugTF<-FALSE
if(debugTF){

	trajdat<-sel.trajdat
	sel.par=sel.par
	boot.trajpath=boot.trajpath
	boot.outname=boot.outname
	storeTF=FALSE

}

bootstrap.traj<-function(trajdat=trajdat, agl.index=agl.index, sel.par=sel.par, boot.trajpath=boot.trajpath, boot.outname=boot.outname, storeTF=FALSE){

# put the bootstrapped traj into /uufs/chpc.utah.edu/common/home/lin-group4/wde/STILT_output/OCO-2/Traj/Riyadh/1km/multi_agl/agl_test/orig_traj, for "weighted.footprint.r" to interpolate the AK and PW profiles

# apply Random Samples and Permutations to particles at each level
# testing...
#sel.par<-50 	# numpar for selected run for each level

# randomly sampling the traj index
# if replace ==TRUE, allows for repetitionin sampling, meaning potentially we can sample same traj for multiple times
#new.index<-sample(orig.par, size=sel.par, replace=TRUE)

cat("bootstrap.traj(): bootstrapping indices and generating random indices...\n")
bootstrap<-function(x, size=sel.par, replace=TRUE){sample(x, size, replace)}
new.index<-tapply(trajdat[,"index"], trajdat[,"level"], bootstrap)	# return lists of index numbers for each selected levels

# now grab the particle according to "new.index" for each selected level, from trajdat
# e.g., for selected level 1, grab the particles whose index name is the same as new.index[[1]]

# define a function to return the row number of y that matches each x, x can be repeated
# and also rename the selected particle indices
find.row.rename<-function(x, y, row.num=NULL, rename.str=NULL, rename.index=rename.index){
	for(i in 1:length(x)){
		tmp.row.num<-which(y==x[i])		# current row.number
		row.num<-c(row.num, tmp.row.num)		# combine

		# also count how many particles that matches
		rep.time<-length(tmp.row.num)

		#create a string for renaming particle index as a next step, combining
		rename.str<-c(rename.str, rep(rename.index[i], rep.time))	# should have the same dimension as row.num
	}

	results<-data.frame(row.num, rename.str)
	colnames(results)<-list("row.number","rename.string")
	return(results)	# results is a list, use "$" to access each list
}

# for debugging...
#df<-by(trajdat, trajdat[,"level"], function(df){return(df)})	# split the entire data frame by level
#df<-df[[2]]

#define a function to subset the trajdat according to new.index after bootstrapping
# df is the input for getSubset, but the output of by (), see below
# df has same type as trajdat
getSubset<-function(df, rename.level=NULL, rename.index=NULL){

	#return(df)	# split by level, return lists of sub.data.frame equal to "tmp.level.trajdat"

	level<-as.character(df[1,"level"])	# grab the level as a string
	tmp.index<-get(level, new.index)	# get(attr,list), attr must be a string

	# rename original levels (base run) to new levels (selected levels)
	rename.level<-which(agl.index==level)

	# rename those index from 1 to sel.par, inputs for "find.row.rename"
	# depending on new levels, as 1, 2, 3, ..., length(sel.agl)
	min.index<-(rename.level-1) * sel.par + 1
	max.index<-rename.level * sel.par
	rename.index<-seq(min.index, max.index, 1)

	# call the find.row.rename function for row numbers and rename strings
	# results will have two colnames, one for selected row.number, one for renaming the traj[,"index"]
	# !!! remember to pass "rename.index" in the find.row.rename() function, EVERYTIME !!!
	results <- find.row.rename(y=df[,"index"], x=tmp.index, row.num=NULL, rename.str=NULL, rename.index=rename.index)

	# now grab the row for the subset
	sub.df<-df[results$row.number,]

	# add another column for rename particle indices and levels after bootstrapping
	sub.df[,"rename.index"]<-results[,"rename.string"]
	# for checking...
	#sub.df[,"rename.level"]<-rep(rename.level, nrow(sub.df))

	return(sub.df)
}

# use by() to subset particles at all selected levels
cat("bootstrap.traj(): subsetting according to selected particles...takes a while\n")
sel.par.trajdat<-by(trajdat, trajdat[,"level"], getSubset, simplify=TRUE)
sel.par.trajdat<-do.call(rbind, sel.par.trajdat)	# convert and merge the "by" object into dataframe
sel.par.trajdat<-as.matrix(sel.par.trajdat)	# then convert to matrix
rownames(sel.par.trajdat)<-NULL

# checking...
#level1<-sel.par.trajdat[sel.par.trajdat[,"level"]==1,]
#table(level1[,"index"])

# need to reorder the subset, by time (-2 min to -4320min), and index ranging from 1-sel.par*length(sel.agl)
new.seltraj<-sel.par.trajdat[order(-sel.par.trajdat[,"time"]),]

# get rid of old level and index and rename new colnames, and turn it back to matrix
rm.col<-which(colnames(new.seltraj) %in% c("level","index"))
new.seltraj2<-new.seltraj[,-rm.col]
colnames(new.seltraj2)[colnames(new.seltraj2)=="rename.index"]<-"index"

# store the resultant traj into .RData files
if(storeTF==TRUE)assignr(value=new.seltraj2, xname=boot.outname, path=boot.trajpath, printTF=TRUE)

# return the subset to main script
return(new.seltraj2)	# still in matrix form

}	# end of subroutine
