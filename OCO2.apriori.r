# Subroutine to calculate the CO2 contribution from a priori CO2
# what needed are the combined a priori CO2 profiles and the adjusted combine AK *PW profiles from "combine.profile" from main script
# Written by Dien WU, 09/16/2016

oco2.apriori<-function(combine.profile=combine.profile){

ak.norm<-combine.profile$ak.norm
pw<-combine.profile$pw

# the CO2 contribution from a priori portion is calculated as the product of (I-AK)* CO2.apriori
# grab normalized AK profiles and apriori CO2 profile
apriori<-combine.profile$oco2.prior
co2.apriori<-(rep(1,length(ak.norm))-ak.norm)*apriori
xco2.apriori<-sum(co2.apriori*pw)


### for plotting modeled boundary conditions
if(F){
	### plot AK.norm
  cex=4.5
	png(filename="int_boundary_condition.tiff",width=2000, height=2400)
	par(mar=c(12, 13, 10, 5), mgp=c(8,3,0))
	plot(oco2.info[[1]]$apriori,oco2.info[[1]]$pres,ylim=c(1000,0),xlim=c(377,407),pch=19,cex=cex*0.6,main="Modeled CO2 boundary conditions vs. OCO-2 prior CO2 profile",xlab="CO2 [ppm]",ylab="Pressure [hpa]",cex.main=cex,cex.lab=cex,cex.axis=cex)
	lines(oco2.info[[1]]$apriori,oco2.info[[1]]$pres,lty=2,lwd=cex-1)

	points(combine.profile[combine.profile$stiltTF==FALSE,"ctnrt.bg"],combine.profile[combine.profile$stiltTF==FALSE,"pres"],col="blue",cex=cex,lwd=cex)
	points(combine.profile[combine.profile$stiltTF,"ctnrt.bg"],combine.profile[combine.profile$stiltTF,"pres"],col="red",cex=cex,lwd=cex)
	text(380,475,"MAXAGL",col="red",cex=cex)
	abline(h=c(upper.oco.pres[1],upper.oco.pres[length(upper.oco.pres)]),col="blue",lty=2,lwd=cex)
	abline(h=c(stilt.pres[1],stilt.pres[length(stilt.pres)]),col="red",lty=2,lwd=cex)
	legend(377,750,c("OCO-2 prior","Boundary conditions at upper model levels","Boundary conditions (endpoint-trajec)\n at lower model levels"),col=c("black","blue","red"),text.col=c("black","blue","red"),bty="n",lty=c(2,NA,NA),pch=c(19,1,1),cex=cex,lwd=c(cex,NA,NA))
	dev.off()
}


return(xco2.apriori)

} # end of subroutine
