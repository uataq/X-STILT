# subroutine to get info from trajec name
# ident can be a vector
# by DW, 07/10/2017

# if aglTF, then interpolate AGLs from trajec
# add STILTv2 trajec name convention, DW, 07/31/2018

ident.to.info <- function(ident, stilt.ver, aglTF = T){

  if (stilt.ver == 1) {  # if using version 1 naming convention
    library(stringr)

    # get rid of .RData, if there is
    if (grepl('.RData', ident[1])) ident <- substr(ident, 1, nchar(ident) - 6)

    # get # of columns
    ncol <- str_count(ident, 'x') + 1
    colnms <- c('year', 'mon', 'day', 'hour', 'lat', 'lon', 'agl', 'numpar', 'dxyp')
    colnms.mm <- append(x = colnms, values = 'min', after = 4)  # if mins in recp time

    # each ident could have different format, thus, do a loop -->
    all.recp <- NULL; all.agl  <- NULL
    for (i in 1 : length(ident)){

      tmp.ncol  <- ncol[i]
      tmp.ident <- ident[i]
      tmp.info  <- unlist(strsplit(tmp.ident, 'x'))
      tmp.info  <- as.data.frame(matrix(tmp.info, ncol = tmp.ncol, byrow = TRUE),
                                 stringsAsFactors = F)

      ##### get receptor info, 'recp.info'
      dxyp <- NA; dpar <- NA; recp.min <- NA
      # if there is mins in recp time string
      if (!is.na(as.numeric(tmp.info$V5))){
        colnames(tmp.info) <- colnms.mm
        recp.min <- as.numeric(tmp.info$min)
        timestr <- paste0(tmp.info$year, tmp.info$mon, tmp.info$day,
                          tmp.info$hour, tmp.info$min)

      } else {
        colnames(tmp.info) <- colnms[1:tmp.ncol]
        timestr <- paste0(tmp.info$year, tmp.info$mon, tmp.info$day, tmp.info$hour)
      }
      if (tmp.ncol >= 9) dxyp <- as.numeric(tmp.info$dxyp)

      # get receptors time info
      recp.year <- as.numeric(tmp.info$year)
      recp.mon  <- as.numeric(tmp.info$mon)
      recp.day  <- as.numeric(tmp.info$day)
      recp.hour <- as.numeric(tmp.info$hour)

      # get rid of N and E for lat lon
      recp.lat <- as.numeric(substr(tmp.info$lat, 1, nchar(tmp.info$lat) - 1))
      recp.lon <- as.numeric(substr(tmp.info$lon, 1, nchar(tmp.info$lon) - 1))

      # if for western or southern hemisphere
      if (grepl('W', tmp.info$lon)) recp.lon <- recp.lon * (-1)
      if (grepl('S', tmp.info$lat)) recp.lat <- recp.lat * (-1)

      npar <- as.numeric(gsub('P', '', tmp.info$numpar))
      agl.info <- as.character(tmp.info$agl)  # get AGL info, 'agl.info'

      #### whether grab height info, column or fixed heights
      if (aglTF) {
        if (grepl('by', agl.info)) {
          if (grepl('&', agl.info)) find.agl <- unlist(strsplit(agl.info, '&'))

          # count different vertical spacing, by counting 'by'
          #count.dh <- nchar(find.agl) - nchar(gsub(agl.info,'','b'))
          # lower AGL levels if unequal dh, dw (02/08/2017)
          # for summer time, we have three different agl
          stilt.agl <- NULL

          for(h in 1 : length(find.agl)){
            # an example of 'find.agl'--'00000-03000by00100'
            ff.agl  <- unlist(strsplit(find.agl[h], '-'))
            ff.agl  <- unlist(strsplit(ff.agl, 'by'))
            max.agl <- as.numeric(ff.agl[2])
            min.agl <- as.numeric(ff.agl[1])
            dh      <- as.numeric(ff.agl[3])
            agl     <- seq(min.agl, max.agl, dh)
            stilt.agl <- c(stilt.agl, agl)
          }  # end for h

          stilt.nlevel <- length(stilt.agl)
          dpar <- npar/stilt.nlevel	# for all trajs

        }else{
          # for const AGL
          stilt.agl <- as.numeric(agl.info)
          stilt.nlevel <- length(stilt.agl)
        }  # end if 'by'

        # may have a warning,when agls are different
        # store AGL info
        all.agl  <- rbind(all.agl, stilt.agl)

        # storing recp.info results...
        recp.info <- data.frame(recp.year, recp.mon, recp.day, recp.hour,
                                recp.min, recp.lat, recp.lon, timestr, dxyp,
                                npar, dpar, stilt.nlevel, stringsAsFactors = F)
        all.recp <- rbind(all.recp, recp.info)

      }else{

        # storing recp.info results...
        recp.info <- data.frame(recp.year, recp.mon, recp.day, recp.hour,
                                recp.min, recp.lat, recp.lon, timestr, dxyp,
                                npar, dpar, stringsAsFactors = F)
        all.recp <- rbind(all.recp, recp.info)
      } # end if aglTF
    } # end loop i

    if(aglTF){
      rownames(all.agl) <- seq(1, length(ident))
      colnames(all.agl) <- seq(1, stilt.nlevel)
    } # end if aglTF

    rownames(all.recp)<- seq(1,length(ident))
    return(list(all.recp, all.agl))

  } else {  # foe STILT version 2 convention

    ident <- gsub('_X_traj.rds', '', ident)
    recp.info <- matrix(unlist(strsplit(ident, '_')), ncol = 3, byrow = T)
    recp.info <- as.data.frame(recp.info)
    colnames(recp.info) <- list('timestr', 'recp.lon', 'recp.lat')

    recp.info <- recp.info %>% mutate_all(funs(factor)) %>%
                               mutate_all(funs(as.character)) %>% 
                               mutate_all(funs(as.numeric))

    return(recp.info)
  } # end if stilt.ver

} # end of subroutine
