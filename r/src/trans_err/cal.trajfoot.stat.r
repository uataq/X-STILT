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

# for debug..
if (F) {
    X = 1
    r_run_time = recp.info$run_time; r_lati = recp.info$lati
    r_long = recp.info$long; r_zagl = recp.info$zagl
}

cal.trajfoot.stat <- function(X, workdir, outdir, emiss.file, met, dpar, ct.ver, 
                              ctflux.path, ctmole.path, r_run_time, r_lati, 
                              r_long, r_zagl, pct = 0.99, max.sd.trans = 200) {

  try({
    print(r_run_time)
        
    # Ensure dependencies are loaded for current node/process
    setwd(workdir); source('r/dependencies.r')
    library(raster); library(dplyr)

    # 1. grab all trajec info
    # new version of STILT stored both trajec before and after perturbations into
    # one single rds file
    traj.path <- file.path(outdir, 'by-id')
    traj.patt <- '_X_wgttraj.rds'
	traj.file <- list.files(traj.path, traj.patt, recursive = T, full.names = T)

    if (length(r_run_time) > 1) {
        r_run_time <- r_run_time[X]
        r_lati  <- format(r_lati[X], digits = 4, nsmall = 4)
        r_long  <- format(r_long[X], digits = 4, nsmall = 4)
        r_zagl  <- r_zagl[X]
    }
    timestr <- substr(format(r_run_time, "%Y%m%d%H%M%S"), 1, 10)

    # get the correct traj file
    traj.file <- traj.file[grep(r_lati, traj.file)]
    if (length(traj.file) == 0) {
        cat('cal.trans.err(): NO trajec found...\n'); return()}
    emiss <- raster(emiss.file)   # get emissions


	#--------------------------------------------------------------------------- #
    #### 2. directly read from exiting weighted trajs -------------------------
    cat(paste('cal.trans.err(): Reading in trajec for', r_lati, 'N...\n'))
    r.all <- readRDS(traj.file)
    p1 <- r.all$particle         # traj without error
    p2 <- r.all$particle_error   # traj with error
    combine.prof <- r.all$wgt.prof


    #### 3. USE non-weighted trajec (orig vs. err) ----------------------------
    # to calculate the dCO2.anthro
    # return data frame with index number and corresponding modeled co2
    cat('working on FFCO2 Contribution...\n')
    co2.ff1 <- ff.trajfoot(trajdat = p1, emiss)
    co2.ff2 <- ff.trajfoot(trajdat = p2, emiss)

    # checking...
    check <- function(x, combine.prof) {
        weight <- data.frame(x, level = rep(seq(1, nrow(x)/100), each = 100))
        colnames(weight)[2] <- 'co2'
        sum.co2 <- sum(tapply(weight$co2, weight$level, mean) *
                   combine.prof[combine.prof$stiltTF, 'ak.pwf'])
        return(sum.co2)
    } # end of function

    xco2.ff1 <- check(x = co2.ff1, combine.prof)
    xco2.ff2 <- check(x = co2.ff2, combine.prof)
    cat(paste('Checking...co2.ff.orig =', signif(xco2.ff1, 4), '[ppm-CO2];',
                          'co2.ff.err =', signif(xco2.ff2, 4), '[ppm-CO2]\n'))


    #### 4. USE non-weighted trajec (orig vs. err) ----------------------------
    # to calculate the dCO2.bio
    cat('working on bio CO2 Contribution...\n')
    co2.bio1 <- bio.trajfoot(trajdat = p1, timestr = timestr, ctflux.path = ctflux.path)
    co2.bio2 <- bio.trajfoot(trajdat = p2, timestr = timestr, ctflux.path = ctflux.path)

    # checking...
    xco2.bio1 <- check(x = co2.bio1, combine.prof)
    xco2.bio2 <- check(x = co2.bio2, combine.prof)
    cat(paste('Checking...bio.orig =', signif(xco2.bio1, 4), '[ppm-CO2];',
                         'bio.err =', signif(xco2.bio2, 4), '[ppm-CO2]\n'))


    ### 5. USE non-weighted trajec --------------------------------------------
    # to grab background CO2 for STILT levels
    cat('working on CO2 endpoints...\n')
    # !!! remember to weight background concentration with Averaging Kernel
    co2.edp1 <- endpts.trajfoot(trajdat = p1, timestr, ctmole.path, combine.prof)
    co2.edp2 <- endpts.trajfoot(trajdat = p2, timestr, ctmole.path, combine.prof)


    ### 6. NOW, sum up all 3 contributions ------------------------------------
    # to calculate CO2.true profiles (w/wout errors) for each particle
    merge.co2 <- co2.ff1 %>% full_join(co2.ff2,  by = 'indx') %>%
                             full_join(co2.bio1, by = 'indx') %>% 
                             full_join(co2.bio2, by = 'indx') %>%
                             full_join(co2.edp1, by = 'indx') %>% 
                             full_join(co2.edp2, by = 'indx') %>%
                             mutate(tot1 = ff.sum.x + bio.sum.x + edp.x,
                                    tot2 = ff.sum.y + bio.sum.y + edp.y) %>%
                             na.omit()

    # assign level index to co2 for each traj
    merge.co2$level <- rep(seq(1, nrow(merge.co2) / dpar), each = dpar)


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
    lr.stat.info <- scale.dvar(co2.stat = co2.stat)

    # store vertical profile of trans errors in the same dir as 'by-id'
    outname <- gsub(traj.patt, '', basename(traj.file))
    stat.file <- file.path(dirname(traj.file), 
                           paste(outname, met, 'info.rds', sep = '_'))

    all.stat.info <- list(merge.co2 = merge.co2, stat.info = lr.stat.info, 
                          combine.prof = combine.prof)
    saveRDS(all.stat.info, file = stat.file)
    cat(paste('saving all error statistic info as in', stat.file, '...\n'))
  })  # try()
} # end of subroutine