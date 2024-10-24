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
# use un-weighted trajec for estimating CO2.true (co2_ff, co2_bio and co2_edp)
#   profiles at each level, DW, 04/20/1027
# add more tracks, DW, 11/25/2017

# generalize as a subroutine, DW, 07/19/2018
# path1, path2 for paths that store original trajec vs. trajec with error component
# add vertical integration of errors, DW
# set an upper limit of sd in ppm, e.g., max.sd.trans = 200
# distribute each receptor to a core for calculating trans err, DW, 10/21/2018 
# get rid of the checking portion, DW, 01/27/2019 
# use non-weighted trajec-level foot for getting errors, DW, 01/28/2019 

# due to momery limit, do not load weighted traj for AK and PW weighting prof, 
#   instead, call get.wgt.funcv3() to get `combine.prof` and
#   remove zero footprint to save time and space in ff/bio_trajfoot.r, DW, 01/29/2019 

cal.trajfoot.stat = function(workdir, output = NULL, outdir, met, emiss.file, 
                             combine.prof, ct.ver, ctflux.path, ctmole.path, 
                             r_run_time, r_lati, r_long, r_zagl, pct = 0.99, 
                             max.sd.trans = 200) {
  try({

    # 1. grab all trajec info
    # new version of STILT stored both trajec before/after perturbations into
    # one single rds file
    timestr = substr(format(r_run_time, "%Y%m%d%H%M%S"), 1, 10)
    traj_dir = file.path(outdir, 'by-id')
    traj_fn = list.files(traj_dir, '_X_traj.rds', recursive = T, full.names = T)

    # get the correct traj file
    traj_fn = traj_fn[grep(r_lati, traj_fn)]

    if (length(traj_fn) == 0) {
        cat('cal.trajfoot.stat(): NO trajec file found...\n'); return()}
    emiss = raster(emiss.file)   # get emissions

	  #------------------------------------------------------------------------- #
    #### 2. directly read from exiting weighted trajs -------------------------
    cat(paste('cal.trajfoot.stat(): Reading in trajec for', r_lati, 'N...\n'))
    if (is.null(output)) output = readRDS(traj_fn) 
    p1 = output$particle         # traj without error
    p2 = output$particle_error   # traj with error
    invisible(gc())

    # get column receptor info from particle lists
    xhgt   = unique(p1$xhgt)
    nlevel = length(xhgt)
    numpar = length(unique(p1$indx))
    dpar   = numpar / nlevel
    
    #### 3. USE non-weighted trajec (orig vs. err) ----------------------------
    # to calculate the dCO2.ff
    # return data frame with index number and corresponding modeled co2
    cat('cal.trajfoot.stat(): working on FFCO2 Contribution...\n')
    co2_ff1 = ff.trajfoot(trajdat = p1, emiss)
    co2_ff2 = ff.trajfoot(trajdat = p2, emiss)
    invisible(gc())

    #### 4. USE non-weighted trajec (orig vs. err) ----------------------------
    # to calculate the dCO2.bio
    cat('cal.trajfoot.stat(): working on bio CO2 Contribution...\n')
    co2_bio1 = bio.trajfoot(trajdat = p1, timestr, ctflux.path)
    co2_bio2 = bio.trajfoot(trajdat = p2, timestr, ctflux.path)
    invisible(gc())

    ### 5. USE non-weighted trajec (orig vs. err) -----------------------------
    # to grab background CO2 for STILT levels
    cat('cal.trajfoot.stat(): working on CO2 endpoints...\n')
    # !!! no need to weight background concentration with Averaging Kernel here
    # because error for each level will be further weighted in cal.trans.err()
    # bug fixed by DW, 01/28/2019
    co2_edp1 = endpts.trajfoot(trajdat = p1, timestr, ctmole.path)
    co2_edp2 = endpts.trajfoot(trajdat = p2, timestr, ctmole.path)
    invisible(gc())

    ### 6. NOW, sum up all 3 contributions ------------------------------------
    # to calculate CO2.true profiles (w/wout errors) for each particle
    # (at this point, CO2 for each particle has not been weighted through AK,PWF)
    co2_df = co2_ff1 %>% full_join(co2_ff2,  by = 'indx') %>%
                         full_join(co2_bio1, by = 'indx') %>% 
                         full_join(co2_bio2, by = 'indx') %>%
                         full_join(co2_edp1, by = 'indx') %>% 
                         full_join(co2_edp2, by = 'indx') %>%
                         mutate(tot1 = ff.sum.x + bio.sum.x + edp.x,
                                tot2 = ff.sum.y + bio.sum.y + edp.y) %>%
                         na.omit() %>% 
                         # also, assign level index to co2 for each traj
                         mutate(level = rep(seq(1, nlevel, 1), each = dpar))

    ### 7. calculate CO2 variances before and after randomizations ------------
    stat_orig = cal.varv2(x = co2_df, colname = 'tot1', pct = pct)
    stat_err  = cal.varv2(x = co2_df, colname = 'tot2', pct = pct)
    colnames(stat_orig) = c('mean.orig', 'sd.orig', 'var.orig', 'level')
    colnames(stat_err)  = c('mean.err', 'sd.err', 'var.err', 'level')

    # merge statistics and calculate diff in variances, i.e., transort errors
    co2_stat = full_join(stat_orig, stat_err, by = 'level') %>% na.omit() %>%
               mutate(dVAR = var.err - var.orig, hgt = unique(p1$xhgt))

    ### 8. An additional step -------------------------------------------------
    # to remove negative trans error and output in txt file
    # scale trans errors based on weighted linear regression lines
    # based on Wu et al., 2018
    lr_stat_info = scale.dvar(co2_stat)

    # store vertical profile of trans errors in the same dir as 'by-id'
    outname = gsub('_X_traj.rds', '', basename(traj_fn))
    stat_file = file.path(dirname(traj_fn), 
                          paste(outname, met, 'emiss_info.rds', sep = '_'))

    # add `combine.prof` for further error accumulation in calc.trans.err()
    all_stat_info = list(merge.co2 = co2_df, stat.info = lr_stat_info, 
                         combine.prof = combine.prof)
    saveRDS(all_stat_info, file = stat_file)
    cat(paste('cal.trajfoot.stat(): saving error statistic info for', r_lati, 'N\n'))
    invisible(gc())

    return(stat_file) # return the name of the stat_file
  })
  
} # end of subroutine