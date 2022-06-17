# calculate diff in variances (without and with errors) at each level for OCO2
# DW, 04/20/2017

# add filtering to remove outliers, for some case:
# anthro with error goes up to >1000ppm, which is impossible, DW, 01/08/2018

# generalize and clear up the code, DW, 07/23/2018

cal.varv2 = function(x, colname, pct = 0.99, func = 'normal'){

	library(MASS)
	uni_levels = unique(x$level)
	var_tot = NULL; mean_tot = NULL; sd_tot = NULL

	# decide to go with test III -- upper 1st percentile
	for (l in 1:length(uni_levels)) {

		tmp = x %>% filter(level == uni_levels[l]) %>% dplyr::select(colname)
        tmp = as.numeric(unlist(tmp))

		# add filtering to remove outliers, DW, 01/08/2018
		# for some case, anthro with error goes up to 2000ppm,
		# which is impossible and skews the normal distribution
		if(F){
			## test I -- some arbitary number, e.g., 700 ppm
			sel1 = tmp[tmp <= 700]

			# test II -- upperest bin
			breaks = hist(tmp, plot = F)$breaks
			crit2  = breaks[length(breaks) - 1]
			sel2   = tmp[tmp <= crit2]
		}

		# decide to go with test III -- remove upper 1st percentile per level
		crit3 = quantile(sort(tmp), pct)
		sel3  = tmp[tmp <= crit3]

		# calculate the mean and sd
		fit_tmp  = fitdistr(sel3, densfun = func)
		mean_tmp = as.numeric(fit_tmp$estimate[1])
		sd_tmp   = as.numeric(fit_tmp$estimate[2])
		var_tmp  = sd_tmp^2

		# storing...
		mean_tot = c(mean_tot, mean_tmp)
		sd_tot   = c(sd_tot, sd_tmp)
		var_tot  = c(var_tot, var_tmp)
	}

	result = data.frame(mean = mean_tot, sd = sd_tot, var = var_tot, 
	                    level = uni_levels)
	return(result)
}
