# update, use SD instead of CI for smooth spline, DW 11/16/2017
# add posterior error for observed XCO2, DW, 01/17/2018
# add the particle uncertainties on simulations, 8%, DW, 02/21/2018

ggplot.lat.series <- function(site, timestr, met, oco2.ver,
  sim.gdas=NULL, all.bg, sim.wrf=NULL, all.obs, lon.lat,
  plotTF=F, rmTF=T, rawTF=F, covTF=T, obs.uncertTF=F){

  library(ggplot2)
  library(zoo)
  library(geosphere)

  sim.gdas <- sim.gdas[sim.gdas$timestr == timestr, ]

  trans.errTF <- 'trans.err' %in% colnames(sim.gdas) & length(sim.gdas$trans.err)>0
  emiss.errTF <- 'emiss.err' %in% colnames(sim.gdas) & length(sim.gdas$emiss.err)>0
  addTF  <- 'add.pp' %in% colnames(sim.gdas)
  gdasTF <- !is.na(match('gdas', met)) | !is.na(match('gdas0p5', met))
  wrfTF  <- !is.na(match('1km', met))

  print(trans.errTF)
  print(emiss.errTF)
  print(addTF)

  ### plot legend and labels
  obs.span <- 0.2; sim.span <- obs.span # few points for simulated, thus larger span
  col.name <- 'XCO2\nConcentrations:'
  fill.name <- 'Uncertainty\nSources:'
  col.labels <- c('1'='WRF-XSTILT modeled XCO2 with smooth spline',
                  '2'='GDAS-XSTILT modeled XCO2 with smooth spline',
                  '3'='GDAS-XSTILT modeled XCO2 with PP10+PP14',
                  '4'='Model-defined background XCO2 (ocean+bio+endpts)',
                  '5'=NULL,
                  '6'='Observed XCO2 (QF=0)', '7'=NULL, '8'=NULL,
                  '9'='Bin-averaged observed XCO2 (QF=0)')

  fill.labels <- c('1'='WRF-XSTILT transport errors',
                   '2'='GDAS-XSTILT transport errors',
                   '3'='Statistical background XCO2 uncertainty',
                   '4'=NULL,
                   '5'='Overpass-specific background uncertainty',
                   '6'='Observed XCO2 (QF=0) uncertainties',
                   '8'='XCO2 uncertainty (due to prior emissions)',
                   '9'='Observed XCO2 (QF=0) uncertainties')

  values <- c('1'='deepskyblue2', '2'='purple', '3'='hotpink', '5'='limegreen',
              '6'='gray70', '8'='orange', '9'='gray20')

  lt.values <- c('1'=2, '2'=2, '3'=2, '4'=2, '5'=2, '6'=NULL, '7'=2, '8'=1, '9'=1)
  pch.values <- c('1'=19, '2'=19, '3'=19, '4'=NULL, '5'=NULL,
                  '6'=17, '7'=NULL, '8'=NULL, '9'=17)

  ## calculate total obs errors

  for(n in 1:nrow(avg.obs)){avg.obs$obs.uncert[n]<-sqrt(sum(avg.obs[n,grepl('sd',colnames(avg.obs))]^2))}
  obs.sort<-avg.obs[order(avg.obs$lat),]

  ylimit <- min.max(x=obs, y=sim.gdas, fixTF=fixTF, trans.errTF=trans.errTF,
                  emiss.errTF=emiss.errTF, addTF=addTF)
  minxco2 <- ylimit[1]; maxxco2 <- ylimit[2]

  #### select background
  dy <- -1*(maxxco2-minxco2)/40
  sel.bg <- data.frame(bg=all.bg$final.bg, bg.uncert=all.bg$poly.bg.sd)

  sel.bg$x.pos <- maxlat+0.1
  sel.bg$y.pos <- sel.bg$bg + dy
  sel.bg$line.col <- 'darkgreen'
  sel.bg$line.lab <- 'XCO2.bg'

  # add a second x-axis for km, DW, 11/16/2017
  # edit distance, d=0 for the least distance from city center to OCO-2 track, DW, 11/21/2017
  if(site=='Riyadh'){center.lat<-24.71;center.lon<-46.74}

  # use foot=2 to get the angle of OCO2 track
  center.obs <- obs[obs$foot == 2,c('lat','lon')]
  d <- NULL
  for(l in 1:nrow(center.obs)){
    d<-c(d,distCosine(c(center.lon,center.lat),c(center.obs$lon[l],center.obs$lat[l])))
  }
  min.obs<-center.obs[d==min(d),]   # determine the center lattitude on OCO-2 track

  # then use the distCosine() to calculate the distance of each sounding from the min.obs, DW, 11/21/2017
  dist <- NULL; slope <- NULL  # use km
  for (l in 1:nrow(center.obs)){dist <- c(dist,distCosine(c(min.obs$lon,min.obs$lat),c(center.obs$lon[l],center.obs$lat[l]))/1000)}
  for (l in 1:nrow(center.obs)){slope<- c(slope,dist[l]/(center.obs$lat[l]-center.obs$lat[which(d==min(d))]))}
  mean.slope <- mean(abs(slope), na.rm=T)

  ### start plotting--
  lwd<-0.9;font.size=rel(1.2)
  p1<-ggplot()+theme_bw()#+geom_vline(xintercept=min.obs$lat,colour='gray60',linetype=2)

  # go with normal stat, format(round(1.1234, 2), nsmall = 2)
  # apply spline to mean values
  if(fixTF)area.sim<-sum(diff(sim.gdas$lat)*rollmean(sim.gdas$sim.total.xco2.fix-sel.bg$bg,2))
  if(fixTF==FALSE)area.sim<-sum(diff(sim.gdas$lat)*rollmean(sim.gdas$sim.total.xco2-sel.bg$bg,2))

  ### whether bin up
  area.obs<-sum(diff(avg.obs$lat)*rollmean(avg.obs$xco2-sel.bg$bg,2))
  #area.obs<-sum(diff(sim.gdas$lat))*rollmean(sim.gdas$obs.total.xco2-sel.bg$bg,2)

  area.trans.err<-NA;  area.emiss.err<-NA
  area.trans.err.exc<-NA;  area.emiss.err.exc<-NA
  area.part.err<-NA

  #### ------------------------------------------------------ calculate ensemble errors ------------------------------------------------------ ####
  # not includ covariances
  sim.gdas$part.err<-sim.gdas$dxco2.anthro*0.08
  if(fixTF==TRUE)frac.part.err<-sqrt(sum(sim.gdas$part.err^2,na.rm=T))/sum(sim.gdas$sim.total.xco2.fix-sel.bg$bg)
  if(fixTF==FALSE)frac.part.err<-sqrt(sum(sim.gdas$part.err^2,na.rm=T))/sum(sim.gdas$sim.total.xco2-sel.bg$bg)
  area.part.err<-frac.part.err*area.sim
  cat(paste('area for ensemble errors:',signif(area.part.err,3),'ppm\n'))
  gdas.part.err<-as.character(paste('+/-',signif(area.part.err,3),sep=''))

  #### ------------------------------------------------------ calculate emission errors ------------------------------------------------------ ####
  if(emiss.errTF){
  	### add emission uncertainties
    p1<-p1+geom_ribbon(data=sim.gdas,aes(x=lat,ymin=sim.total.xco2.fix-dxco2.anthro.uncert,ymax=sim.total.xco2.fix+dxco2.anthro.uncert,fill=as.factor(8)),alpha=0.3)
  }

  #### ------------------------------------------------------  calculate trans err ------------------------------------------------------- ####
  if(trans.errTF){
    #### add WRF & GDAS transport error uncertainties
    p1<-p1+geom_ribbon(data=sim.gdas,aes(x=lat,ymin=sim.total.xco2.fix-trans.err,ymax=sim.total.xco2.fix+trans.err,fill=as.factor(2)),alpha=0.3)
  }

  #### ------------------------------------------------------  add title, labs ------------------------------------------------------ ##
  title<-paste('Latitude Series of Simulated vs. Observed XCO2 from Lite',oco2.version,'[ppm] for',site,'for track on',timestr)
  if(addTF)title<-paste(title,'with PP10+PP14')
  title.bg<-paste('Modeled background XCO2:',signif(all.bg$sim.bg,5),'ppm; M2H background (daily median):',signif(all.bg$hakka.bg,5),'ppm; M2S background (mean-sd):',signif(all.bg$silva.bg,5),'ppm\n Out-of-plume defined background:',signif(all.bg$poly.bg,5),'ppm')
  title.sig<-paste('Simulated urban signal (AUC):',signif(all.bg$sim.sig,3),'ppm for sim-GDAS; Observed urban signal (AUC):',signif(all.bg$obs.sig,3),'+/-',signif(all.bg$obs.err),'ppm')
  if(length(met)==2)title.sig<-paste('Simulated urban signal (AUC):',signif(all.bg$sim.sig,3),'ppm for sim-GDAS and',signif(all.bg$wrf.sig,3),'ppm for sim-WRF; Observed urban signal (AUC):',signif(all.bg$obs.sig,3),'ppm')
  if(trans.errTF)title.err<-paste('Transport Uncertainty (AUC):',gdas.trans.err,'ppm for sim-GDAS')
  if(emiss.errTF)title.err<-paste('Emission Uncertainty (AUC):',gdas.emiss.err,'ppm')
  if(trans.errTF&emiss.errTF)title.err<-paste('Transport Uncertainty (AUC):',gdas.trans.err,'ppm for sim-GDAS; Emission Uncertainty (AUC):',gdas.emiss.err,'ppm')
  title<-paste(title,'\n',title.bg,'\n',title.sig)
  if(trans.errTF|emiss.errTF)title<-paste(title,'\n',title.err)

  #title<- NULL
  p1<-p1+labs(title=title, x='Latitude', y='XCO2 [ppm]')

  # add background uncertainty
  bg.sd<-data.frame(lat=seq(minlat,maxlat,0.1))
  bg.sd$sd<-rep(sel.bg$bg.uncert,length(bg.sd$lat))
  bg.sd$mean<-rep(sel.bg$bg, length(bg.sd$lat))
  p2<-p1+geom_ribbon(data=bg.sd, aes(x=lat,ymin=mean-sd,ymax=mean+sd,fill=as.factor(5)),alpha=0.2)

  ### ------------------------------------------------------  add avg obs and SD ------------------------------------------------------ ##
  if(obs.uncertTF)p2<-p2+geom_ribbon(data=avg.obs,aes(x=lat,ymin=xco2-obs.uncert, ymax=xco2+obs.uncert,fill=as.factor(9)),alpha=0.2)
  p2<-p2+geom_point(data=obs,aes(x=lat,y=xco2,color=as.factor(6),shape=as.factor(6)),size=3)#+geom_smooth(data=avg.obs,aes(x=lat,y=xco2,color=as.factor(6)),span=sim.span,se=F,method='loess')
  p2<-p2+geom_point(data=avg.obs,aes(x=lat,y=xco2,color=as.factor(9),shape=as.factor(9)),size=3)
  p2<-p2+geom_smooth(data=avg.obs,aes(x=lat,y=xco2,color=as.factor(9),linetype=as.factor(9)),span=sim.span,se=F,method='loess')
  #p2<-p2+geom_line(data=avg.obs,aes(x=lat,y=xco2,color=as.factor(9),linetype=as.factor(9)))

  #p2<-p1+geom_ribbon(data=sim.gdas,aes(x=lat,ymin=obs.total.xco2-obs.uncert, ymax=obs.total.xco2+obs.uncert,fill=as.factor(9)),alpha=0.2)
  #p2<-p2+geom_point(data=sim.gdas,aes(x=lat,y=obs.total.xco2,color=as.factor(9),shape=as.factor(9)),size=3)+geom_smooth(data=sim.gdas,aes(x=lat,y=obs.total.xco2,color=as.factor(9)),span=sim.span,se=F,method='loess')

  ### ------------------------------------------------------ add WRF-sim ------------------------------------------------------ ###
  p3<-p2
  if(wrfTF & fixTF==FALSE)p3<-p3+geom_point(data=sim.wrf,aes(x=lat,y=sim.total.xco2,colour=as.factor(1),shape=as.factor(1)),size=3)+geom_smooth(data=sim.wrf,aes(x=lat,y=sim.total.xco2,colour=as.factor(1),linetype=as.factor(1)),span=sim.span,se=F,method='loess')#+geom_line(data=sim.wrf,aes(x=lat,y=sim.total.xco2,colour=as.factor(1),linetype=as.factor(1)))#
  if(wrfTF & fixTF==TRUE)p3<-p3+geom_point(data=sim.wrf,aes(x=lat,y=sim.total.xco2.fix,colour=as.factor(1),shape=as.factor(1)),size=3)+geom_smooth(data=sim.wrf,aes(x=lat,y=sim.total.xco2.fix,colour=as.factor(1),linetype=as.factor(1)),span=sim.span,se=F,method='loess')#+geom_line(data=sim.wrf,aes(x=lat,y=sim.total.xco2.fix,colour=as.factor(1),linetype=as.factor(1)))#

  ### ------------------------------------------------------ add GDAS-sim ------------------------------------------------------ ###
  ### add additional power plant
  p4<-p3
  if(gdasTF & fixTF==FALSE & addTF)p4<-p4+geom_point(data=sim.gdas,aes(x=lat,y=sim.total.xco2+add.pp,colour=as.factor(3),shape=as.factor(3)),size=3,alpha=0.7)#+geom_line(data=sim.gdas,aes(x=lat,y=sim.total.xco2+add.pp,colour=as.factor(3),linetype=as.factor(3)),size=lwd,alpha=0.7)
  if(gdasTF & fixTF==TRUE  & addTF)p4<-p4+geom_point(data=sim.gdas,aes(x=lat,y=sim.total.xco2.fix+add.pp,colour=as.factor(3),shape=as.factor(3)),size=3,alpha=0.7)#+geom_line(data=sim.gdas,aes(x=lat,y=sim.total.xco2.fix+add.pp,colour=as.factor(3),linetype=as.factor(3)),size=lwd,alpha=0.7)
  if(gdasTF & fixTF==FALSE)p4<-p4+geom_point(data=sim.gdas,aes(x=lat,y=sim.total.xco2,colour=as.factor(2),shape=as.factor(2)),size=3,alpha=0.7)+geom_smooth(data=sim.gdas,aes(x=lat,y=sim.total.xco2,colour=as.factor(2),linetype=as.factor(2)),span=sim.span,se=F)#+geom_line(data=sim.gdas,aes(x=lat,y=sim.total.xco2,colour=as.factor(2),linetype=as.factor(2)))#
  if(gdasTF & fixTF==TRUE)p4<-p4+geom_point(data=sim.gdas,aes(x=lat,y=sim.total.xco2.fix,colour=as.factor(2),shape=as.factor(2)),size=3,alpha=0.7)+geom_smooth(data=sim.gdas,aes(x=lat,y=sim.total.xco2.fix,colour=as.factor(2),linetype=as.factor(2)),span=sim.span,se=F)#+geom_line(data=sim.gdas,aes(x=lat,y=sim.total.xco2.fix,colour=as.factor(2),linetype=as.factor(2)))#

  ### ------------------------------------------------------  only plot final background, based on polygon bg ------------------------------------------------------ ###
	p5<-p4+geom_hline(size=1,data=sel.bg,aes(yintercept=bg),linetype=2,colour=sel.bg$line.col)
  #+geom_hline(size=1,aes(yintercept=all.bg$obs.bg),colour='darkgreen',linetype=2)
  p5<-p5+annotate('text',x=sel.bg$x.pos,y=sel.bg$y.pos,size=5,label=sel.bg$line.lab,colour=sel.bg$line.col)

	#### ----------------  add styles for colour, fill, shape, linetype, add theme style ------------------------- ###
  p6<-p5
  p6<-p6+scale_shape_manual(name=col.name,values=pch.values,labels=col.labels)
  p6<-p6+scale_linetype_manual(name=col.name,values=lt.values,labels=col.labels)
	p6<-p6+scale_fill_manual(name=fill.name,values=values,labels=fill.labels)
  p6<-p6+scale_colour_manual(name=col.name,values=values,labels=col.labels)

  p7<-p6+theme(legend.position='bottom',legend.text=element_text(size=font.size),
               legend.key = element_blank(),panel.grid.minor=element_blank(),
               axis.title.y=element_text(size=font.size,angle = 90),
               axis.title.x=element_text(size=font.size,angle=0),
               axis.text=element_text(size=font.size),axis.ticks=element_line(size=font.size),
               title=element_text(size=font.size))
	p7<-p7+guides(colour=guide_legend(ncol=1,order=1),fill=guide_legend(ncol=1,order=2),
                shape=guide_legend(ncol=1,order=1),linetype=FALSE)
  p7<-p7+scale_y_continuous(limits=c(minxco2,maxxco2),breaks=seq(minxco2,maxxco2,1))
  p7<-p7+scale_x_continuous(limits=c(minlat,maxlat),breaks=seq(22,28,seplat),
                            labels=paste(seq(22,28,seplat),'ºN',sep=''),
                            sec.axis=sec_axis(~(.-min.obs$lat)*mean.slope,name='Distance [km]',
                            breaks=seq(-160,160,20),labels=seq(-160,160,20)))

  #### whether plot it and store it
  if(plotTF){
    picname<-paste('LS_',site,'_',timestr,'_',oco2.version,'_fix',substr(fixTF,1,1),'_raw',substr(rawTF,1,1),'.png',sep='')

    if(addTF)picname<-paste('LS_',site,'_',timestr,'_',oco2.version,'_fix',substr(fixTF,1,1),'_raw',substr(rawTF,1,1),'_PP.png',sep='')
    if(trans.errTF)picname<-paste('LS_uncert_',site,'_',timestr,'_',oco2.version,'_fix',substr(fixTF,1,1),'_raw',substr(rawTF,1,1),'.png',sep='')
    if(trans.errTF&addTF)picname<-paste('LS_uncert_',site,'_',timestr,'_',oco2.version,'_fix',substr(fixTF,1,1),'_raw',substr(rawTF,1,1),'_PP.png',sep='')
    if(trans.errTF&rmTF)picname<-paste('LS_uncert_',site,'_',timestr,'_',oco2.version,'_fix',substr(fixTF,1,1),'_raw',substr(rawTF,1,1),'_rm.png',sep='')
    if(trans.errTF&rmTF&addTF)picname<-paste('LS_uncert_',site,'_',timestr,'_',oco2.version,'_fix',substr(fixTF,1,1),'_raw',substr(rawTF,1,1),'_rm_PP.png',sep='')
    print(picname)

    height=9;if(length(title)==0)height=8
    ggsave(p7,filename=picname,width=16,height=height)
  }

  # return plot and uncertainties at overpass level
  uncert<-list(p7,c(trans.err=area.trans.err,emiss.err=area.emiss.err,
                    part.err=area.part.err,trans.err.exc=area.trans.err.exc,
                    emiss.err.exc=area.emiss.err.exc))
  return(uncert)
}
#### end of script














### now create data frame for covariance, spatial cov E
if(F){
  # calculate the distance between points first
  num.recp<-seq(1,nrow(sim.gdas))  # number of soundings/receptors
  grid.dist<-array(0,dim=c(max(num.recp),max(num.recp)),dimnames=list(num.recp,num.recp))
  for(n in 1:length(num.recp)){
    point1<-data.frame(lon=sim.gdas[n,'lon'],lat=sim.gdas[n,'lat'])
    point2<-data.frame(lon=sim.gdas[,'lon'], lat=sim.gdas[,'lat'])
    tmp.dist<-distCosine(p1=point1,p2=point2) # in m
    grid.dist[n,]<-tmp.dist
  }
  sp.cov.emiss<-exp(-grid.dist/1000/Lx.emiss)

  # add covariance as well
  var.emiss.err<-array(0,dim=c(max(num.recp),max(num.recp)),dimnames=list(num.recp,num.recp))
  diag(var.emiss.err)<-sim.gdas$emiss.err  # SD in ppm

  prod.emiss<-var.emiss.err%*%sp.cov.emiss%*%t(var.emiss.err)
  pp<-plot.matrix(x=prod.emiss,timestr=timestr,sim.gdas=sim.gdas,Lx=Lx.emiss)
  #ggsave(pp,filename=paste('sp_cov_emiss_err_nxn_',site,'_',timestr,'.png',sep=''),width=10,height=11)

  tot.sd<-sqrt(sum(prod.emiss,na.rm=T)) # total sd
  if(fixTF)frac.emiss.err<-tot.sd/sum(sim.gdas$sim.total.xco2.fix-sel.bg$bg)
  if(fixTF==FALSE)frac.emiss.err<-tot.sd/sum(sim.gdas$sim.total.xco2-sel.bg$bg)
  area.emiss.err<-frac.emiss.err*area.sim # in ppm


  ### now create data frame for covariance, spatial cov E
  # calculate the distance between points first
  num.recp<-seq(1,nrow(sim.gdas))  # number of soundings/receptors
  grid.dist<-array(0,dim=c(max(num.recp),max(num.recp)),dimnames=list(num.recp,num.recp))
  for(n in 1:length(num.recp)){
    point1<-data.frame(lon=sim.gdas[n,'lon'],lat=sim.gdas[n,'lat'])
    point2<-data.frame(lon=sim.gdas[,'lon'], lat=sim.gdas[,'lat'])
    tmp.dist<-distCosine(p1=point1,p2=point2) # in m
    grid.dist[n,]<-tmp.dist
  }
  sp.cov.trans<-exp(-grid.dist/1000/Lx.trans)   # fix lat first, looping over long

  # add covariance as well
  #tot.var<-sqrt(sum(sim.gdas$trans.err^2,na.rm=T))
  var.trans.err<-array(0,dim=c(max(num.recp),max(num.recp)),dimnames=list(num.recp,num.recp))
  diag(var.trans.err)<-sim.gdas$trans.err  # SD in ppm

  # variance * correlation * variance^T
  plot.matrix<-function(x=x,timestr=timestr,sim.gdas=sim.gdas,Lx=Lx.trans){
    font.size<-rel(1.3)
    melt.x<-melt(x)[melt(x)$value>0,]
    colnames(melt.x)<-list('x','y','value')
    # add latitude for soundings
    add.lat<-data.frame(x=seq(1,dim(x)[1]), y=seq(1,dim(x)[2]), lat=sim.gdas$lat)
    sel.lat<-add.lat[add.lat$x%%5==0,]

    c1<-ggplot()+theme_bw()+geom_raster(data=melt.x,aes(x=x,y=y,fill=value))+scale_fill_gradient(low='gray95',high='purple',name='Variance-covariance\nof transport errors [ppm^2]')
    c2<-c1+labs(title=paste('Spatial covariance matrix [ppm^2] of errors for Riyadh for',timestr,'\nusing correlation length scale of',Lx,'km'),x=paste('# of soundings/receptors (',dim(x)[1],')'),y=paste('# of soundings/receptors (',dim(x)[2],')'))
    c2<-c2+theme(legend.position='bottom',legend.key.width=unit(3.5, 'cm'),legend.key.height=unit(0.5, 'cm'),legend.text=element_text(size=font.size),legend.key = element_blank(),panel.grid.minor=element_blank(),
                axis.title.y=element_text(size=font.size,angle = 90),axis.title.x=element_text(size=font.size,angle=0),axis.text=element_text(size=font.size),axis.ticks=element_line(size=font.size),title=element_text(size=font.size))
    c3<-c2+geom_abline(intercept=0,slope=1,linetype=2,colour='gray30')+geom_text(data=sel.lat,aes(x=x,y=y,label=lat),colour='gray30',size=5)
    return(c3)
  }

  prod.trans<-var.trans.err%*%sp.cov.trans%*%t(var.trans.err)
  pp<-plot.matrix(x=prod.trans,timestr=timestr,sim.gdas=sim.gdas,Lx=Lx.trans)
  ggsave(pp,filename=paste('sp_cov_trans_err_nxn_',site,'_',timestr,'.png',sep=''),width=10,height=11)

  tot.sd<-sqrt(sum(prod.trans,na.rm=T))
  if(fixTF)frac.trans.err<-tot.sd/sum(sim.gdas$sim.total.xco2.fix-sel.bg$bg)
  if(fixTF==FALSE)frac.trans.err<-tot.sd/sum(sim.gdas$sim.total.xco2-sel.bg$bg)
  area.trans.err<-frac.trans.err*area.sim # in ppm

}
