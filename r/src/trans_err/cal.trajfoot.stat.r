#' script to obtain the dCO2 for each trajectories and variances of dCO2
#' need existing traj, need ak, pw profiles from oco2
#' please refer to Wu et al., 2018 GMD for modifications on column trans error
#' @author: Dien Wu, 01/11/2017

#' @variables:
# pct = 0.99: only select lower 99th percentile values, 
#   remove upper 1st percentile for large value

#' @updates:
# merge two traject (2 error stat), 03/15/2017
# add biospheric and background portions, 04/11/2017
# use un-weighted trajec for estimating CO2.true (co2.ff, co2.bio and co2.edp)
#   profiles at each level, DW, 04/20/1027
# add more tracks, DW, 11/25/2017

# generalize as a subroutine, DW, 07/19/2018
# path1, path2 for paths that store original trajec vs. trajec with error component
# add vertical integration of errors, DW
# set an upper limit of sd in ppm, e.g., max.sd.trans <- 200
# distribute each receptor to a core for calculating trans err, DW, 10/21/2018 
# get rid of the checking portion, DW, 01/27/2019 
# use non-weighted trajec-level foot for getting errors, DW, 01/28/2019 

cal.trajfoot.stat <- function(workdir, outdir, emiss.file, met, ct.ver, 
                              ctflux.path, ctmole.path, r_run_time, r_lati, 
                              r_long, r_zagl, pct = 0.99, max.sd.trans = 200) {

    print(r_run_time)
        
    # Ensure dependencies are loaded for current node/process
    setwd(workdir); source('r/dependencies.r'); library(raster); library(dplyr)

    # 1. grab all trajec info
    # new version of STILT stored both trajec before and after perturbations into
    # one single rds file
    timestr <- substr(format(r_run_time, "%Y%m%d%H%M%S"), 1, 10)
    traj.path <- file.path(outdir, 'by-id')
	wgt.traj.file <- list.files(traj.path, '_X_wgttraj.rds', recursive = T, 
                                full.names = T)
    traj.file <- list.files(traj.path, '_X_traj.rds', recursive = T, 
                            full.names = T)

    # get the correct traj file
    traj.file <- traj.file[grep(r_lati, traj.file)]
    wgt.traj.file <- wgt.traj.file[grep(r_lati, wgt.traj.file)]

    if (length(traj.file) == 0 | length(wgt.traj.file) == 0) {
        cat('cal.trajfoot.stat(): NO trajec file found...\n'); return()}
    emiss <- raster(emiss.file)   # get emissions

	#------------------------------------------------------------------------- #
    #### 2. directly read from exiting weighted trajs -------------------------
    cat(paste('cal.trajfoot.stat(): Reading in trajec for', r_lati, 'N...\n'))
    r.wgt <- readRDS(wgt.traj.file)
    r.all <- readRDS(traj.file)

    combine.prof <- r.wgt$wgt.prof
    p1 <- r.all$particle         # traj without error
    p2 <- r.all$particle_error   # traj with error
        
    # get column receptor info from particle lists
    xhgt <- unique(p1$xhgt)
    nlevel <- length(xhgt)
    numpar <- length(unique(p1$indx))
    dpar <- numpar / nlevel
    
    #### 3. USE non-weighted trajec (orig vs. err) ----------------------------
    # to calculate the dCO2.ff
    # return data frame with index number and corresponding modeled co2
    cat('cal.trajfoot.stat(): working on FFCO2 Contribution...\n')
    co2.ff1 <- ff.trajfoot(trajdat = p1, emiss)
    co2.ff2 <- ff.trajfoot(trajdat = p2, emiss)

    #### 4. USE non-weighted trajec (orig vs. err) ----------------------------
    # to calculate the dCO2.bio
    cat('cal.trajfoot.stat(): working on bio CO2 Contribution...\n')
    co2.bio1 <- bio.trajfoot(trajdat = p1, timestr, ctflux.path)
    co2.bio2 <- bio.trajfoot(trajdat = p2, timestr, ctflux.path)

    ### 5. USE non-weighted trajec (orig vs. err) -----------------------------
    # to grab background CO2 for STILT levels
    cat('cal.trajfoot.stat(): working on CO2 endpoints...\n')
    # !!! no need to weight background concentration with Averaging Kernel here
    # because error for each level will be further weighted in cal.trans.err()
    # bug fixed by DW, 01/28/2019
    co2.edp1 <- endpts.trajfoot(trajdat = p1, timestr, ctmole.path)
    co2.edp2 <- endpts.trajfoot(trajdat = p2, timestr, ctmole.path)

    ### 6. NOW, sum up all 3 contributions ------------------------------------
    # to calculate CO2.true profiles (w/wout errors) for each particle
    # (at this point, CO2 for each particle has not been weighted through AK,PWF)
    merge.co2 <- co2.ff1 %>% full_join(co2.ff2,  by = 'indx') %>%
                             full_join(co2.bio1, by = 'indx') %>% 
                             full_join(co2.bio2, by = 'indx') %>%
                             full_join(co2.edp1, by = 'indx') %>% 
                             full_join(co2.edp2, by = 'indx') %>%
                             mutate(tot1 = ff.sum.x + bio.sum.x + edp.x,
                                    tot2 = ff.sum.y + bio.sum.y + edp.y) %>%
                             na.omit() %>% 
                             # also, assign level index to co2 for each traj
                             mutate(level = rep(seq(1, nlevel, 1), each = dpar))

    ### 7. calculate CO2 variances before and after randomizations ------------
    stat.orig <- cal.varv2(x = merge.co2, colname = 'tot1', pct = pct)
    stat.err  <- cal.varv2(x = merge.co2, colname = 'tot2', pct = pct)
    colnames(stat.orig) <- c('mean.orig', 'sd.orig', 'var.orig', 'level')
    colnames(stat.err)  <- c('mean.err', 'sd.err', 'var.err', 'level')

    # merge statistics and calculate diff in variances, i.e., transort errors
    co2.stat <- full_join(stat.orig, stat.err, by = 'level') %>% na.omit() %>%
                mutate(dVAR = var.err - var.orig, hgt = unique(p1$xhgt))


    ### 8. An additional step ----------------------------------------------------
    # to remove negative trans error and output in txt file
    # scale trans errors based on weighted linear regression lines
    # based on Wu et al., 2018
    lr.stat.info <- scale.dvar(co2.stat)

    # store vertical profile of trans errors in the same dir as 'by-id'
    outname <- gsub('_X_traj.rds', '', basename(traj.file))
    stat.file <- file.path(dirname(traj.file), 
                           paste(outname, met, 'emiss_info.rds', sep = '_'))

    # add `combine.prof` for further error accumulation in calc.trans.err()
    all.stat.info <- list(merge.co2 = merge.co2, stat.info = lr.stat.info, 
                          combine.prof = combine.prof)
    saveRDS(all.stat.info, file = stat.file)
    cat(paste('cal.trajfoot.stat(): saving all error statistic info as in', 
               stat.file, '...\n'))
} # end of subroutine