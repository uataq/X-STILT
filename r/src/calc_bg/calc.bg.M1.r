#' subroutine to compute M1 background in Wu et al., 2018
#' @author: Dien Wu, 09/16/2018

#' require: 
#' 1x1 deg hourly footprint, column weighted trajectories  
#' (from 'create_namelist_oco2-xsilt.r')
#' paths and data for CarbonTracker fluxes and mole fractions

calc.bg.M1 <- function(all.timestr, foot.path, traj.path, ct.ver, flux.path, 
                       mf.path, output.path, oco2.ver, txtfile, writeTF = F, 
                       nhrs = -72) {

  bg <- NULL 
  for (tt in 1:length(all.timestr)) {
    
    timestr <- all.timestr[tt]
    cat(paste('Working on overpass on ', timestr, '...\n'))

    fn <- file.path(output.path, 
                    paste0(timestr, '_', site, '_XCO2_CT_', oco2.ver, '.txt'))

    # check whether simulation has been done
    if (file.exists(fn)) {
      xco2.bg <- read.table(file = fn, sep = ',', header = T)

    } else {
      cat('NO existing txtfile found, start simulation...\n')

      # grab weighted trajec info and 1x1 deg footprint info
      # if no foot, use 'create_namelist_oco2-xstilt.r' to generate
      foot.file <- list.files(path = foot.path[tt], pattern = '_1x1_foot.nc', 
                              recursive = T, full.names = T)
      foot.ext <- extent(raster(foot.file[1]))
      
      # paths for weighted trajec
      traj.file <- list.files(path = traj.path[tt], pattern = 'X_wgttraj.rds', 
                              recursive = T, full.names = T)    

      # call function to get XCO2.bio and XCO2.ocean, traj-edpt and aprior
      cat('Couple footprint with CT fluxes...\n')
      xco2.ct <- foot.ctnrt(foot.file, ct.ver[tt], flux.path[tt], timestr, nhrs, 
        writeTF = F)
      
      cat('Retrieve CT CO2 values at trajec-endpoints...\n')
      xco2.bound <- trajec.edpt.ctnrt(traj.file, foot.ext, ct.ver[tt], 
        mf.path[tt], timestr, writeTF = F)

      # merge two results and store in txt file
      xco2.bg <- full_join(xco2.ct, xco2.bound, by = c('timestr', 'lon', 'lat')) %>% 
                 mutate(bg = xco2.bio + xco2.ocn + xco2.ap + xco2.bound)

      write.table(x = xco2.bg, file = fn, sep = ',', row.names = F, quote = F)
    } # end if file.exists
    
    mean.xco2.bg <- mean(xco2.bg$bg, na.rm = T)
    bg <- c(bg, mean.xco2.bg)
  }  # end for tt

  bg.info <- data.frame(timestr = all.timestr, ct.bg = bg)
  if (writeTF) 
    write.table(x = bg.info, file = txtfile, sep = ',', row.names = F, quote = F)
  
  return(bg.info)
}