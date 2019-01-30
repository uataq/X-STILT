#' Main script to obtain the CO2 error for each receptor
#' need existing error statistics info, need ak, pw profiles from oco2
#' @author: Dien Wu, 01/11/2017

#' @variables: 
#' traj.path: path that stores trajec (that include a 'by-id' directory)
#' store.path: path that store the txtfile (result)
#' Lh: vertical error covariance length scale for transport errors, in meter

#' @updates:
#' add vertical integration of errors, DW
#' set an upper limit of sd in ppm, e.g., max.sd.trans <- 200, DW

cal.trans.err <- function(site, timestr, workdir, traj.path, store.path, met, 
                          max.sd.trans = 200, Lh = NULL) {

	# txt file name for outputting model error results
	txtfile <- file.path(store.path, 
	                     paste0('horerr_', timestr, '_', site,'_', met, '.txt'))

	traj.patt <- '_X_wgttraj.rds'
	traj.path <- file.path(outdir, 'by-id')
	traj.file <- list.files(traj.path, traj.patt, recursive = T, full.names = T)

	outname   <- gsub(traj.patt, '', basename(traj.file))
	recp.info <- strsplit.to.df(outname, colnms = c('timestr', 'lon', 'lat'))
	recp.info <- recp.info %>% mutate_all(funs(as.numeric), colnames(recp.info)) %>% 
                 mutate(xco2.trans = NA) %>% dplyr::arrange(lat)

	#------------------------------------------------------------------------- #
	### loop over all unique lat lon and merge two traj
	cat('cal.trans.err(): start vertically integrating trans error...\n')
	for (i in 1 : length(traj.file)) {

		# get error statistic info rds file in each by-id
		stat.file <- list.files(dirname(traj.file[i]), 
		                        paste0(met, '_emiss_info.rds'), full.names = T)

		if (length(stat.file) == 0) {
			cat(paste('NO rds file found for horizontal trans error stat for',
			    recp.info$lat[i], 'N, check `by-id` and skip..\n'))
		    next
		} 
		
		stat.all <- readRDS(file = stat.file)
		stat.info <- stat.all$stat.info
		co2.stat  <- stat.info[[1]]
		sd.trans  <- co2.stat$sd.trans

		combine.prof <- stat.all$combine.prof
		model.prof   <- combine.prof %>% filter(stiltTF == T) %>% 
										 mutate(sd.trans = sd.trans)
		
		### 9. Vertically integrate transport errors for each recp --------
		# up to now, we get a transport error per release level
		# need to vertically integrate errors with correlation length scale
		if (is.null(Lh)) Lh <- 356  # empirical mean length scale

		# w for total weighting function (AK*PWF)
		w <- model.prof[sd.trans < max.sd.trans, 'ak.pwf']
		x <- co2.stat[sd.trans < max.sd.trans, ]

		# call cal.var.cov() for vertical-integrated trans error variance in ppm^2
		prod.trans <- cal.var.cov(sd = sd.trans[sd.trans < max.sd.trans], 
								  L = Lh, w = w, x = x, type = "ver")

		# multiple grid.wgt with prod.trans
		# total SD for trans error for each receptor
		# aggregated over vertial levels
		tot.trans.sd <- sqrt(sum(prod.trans, na.rm=T))
		print(tot.trans.sd)		# now in ppm

		# store trans
		recp.info$xco2.trans[i] <- tot.trans.sd	# in ppm 
	} # end of lat loop i

	# write in txtfile
	recp.info <- recp.info %>% na.omit()
    write.table(recp.info, file = txtfile, sep = ',', quote = F, row.names = F)
	cat(paste('cal.trans.err(): receptor-level trans error has been written in', 
	          txtfile, '...\n'))

    return(recp.info)
} # end of subroutine