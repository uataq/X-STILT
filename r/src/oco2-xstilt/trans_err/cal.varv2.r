# calculate diff in variances (without and with errors) at each level for OCO2
# DW, 04/20/2017

# add filtering to remove outliers, for some case:
# anthro with error goes up to >1000ppm, which is impossible, DW, 01/08/2018

# generalize and clear up the code, DW, 07/23/2018

cal.varv2 <- function(x, colname, pct = 0.99, func = 'normal'){

	library(MASS)
	uni.level <- unique(x$level)
	var.tot <- NULL; mean.tot <- NULL; sd.tot <- NULL

	# decide to go with test III -- upper 1st percentile
	for (l in 1:length(uni.level)) {

		tmp <- x %>% filter(level == uni.level[l]) %>% dplyr::select(colname)
        tmp <- as.numeric(unlist(tmp))

		# add filtering to remove outliers, DW, 01/08/2018
		# for some case, anthro with error goes up to 2000ppm,
		# which is impossible and skews the normal distribution
		if(F){
			## test I -- some arbitary number, e.g., 700 ppm
			sel1 <- tmp[tmp <= 700]

			# test II -- upperest bin
			breaks <- hist(tmp, plot = F)$breaks
			crit2  <- breaks[length(breaks) - 1]
			sel2   <- tmp[tmp <= crit2]
		}

		# decide to go with test III -- remove upper 1st percentile per level
		crit3 <- quantile(sort(tmp), pct)
		sel3  <- tmp[tmp <= crit3]

		# calculate the mean and sd
		fit.tmp  <- fitdistr(sel3, densfun = func)
		mean.tmp <- as.numeric(fit.tmp$estimate[1])
		sd.tmp   <- as.numeric(fit.tmp$estimate[2])
		var.tmp  <- sd.tmp ^ 2

		# storing...
		mean.tot <- c(mean.tot, mean.tmp); sd.tot   <- c(sd.tot, sd.tmp)
		var.tot  <- c(var.tot, var.tmp)
	}

	result <- data.frame(mean = mean.tot, sd = sd.tot, var = var.tot, 
	                     level = uni.level)
	return(result)
}
