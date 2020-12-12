## subroutine to fit linear regression and get ggplots for linear regression plots
# outputs are linear regression line, CO2 statiatics, and ggplot
# see Appendix of Wu et al., 2018, GMD for details
# written by DW, 03/15/2018

# generalize and clear up codes, DW, 07/23/2018
# originate from 'get.ggplot.lr.r' in my stiltR,
# original methods c("all","1.5max","10ppm","wgt"), only use wgt
# change co2.stat to stat, DW, 01/29/2019 

scale.dvar <- function(stat){

	# M4, properly weighted by variances of error
	# select positive/negative error variances diff
	stat.pos <- stat[stat$dVAR > 0, ]
	stat.neg <- stat[stat$dVAR < 0, ]

	wgt.pos <- 1 / stat.pos$var.err  # positive err var
	wgt.neg <- 1 / stat.neg$var.err

    # fit two weighted linear regressions
	fit.lm.pos <- suppressMessages(
		lm(stat.pos[, "var.err"] ~ stat.pos[, "var.orig"], weights = wgt.pos))

	fit.lm.neg <- suppressMessages(
		lm(stat.neg[, "var.err"] ~ stat.neg[, "var.orig"], weights = wgt.neg))

	fit.lm <- data.frame(s = as.numeric(c(fit.lm.pos$coefficients[2], 
	                                      fit.lm.neg$coefficients[2])),
		                 i = as.numeric(c(fit.lm.pos$coefficients[1], 
						                  fit.lm.neg$coefficients[1])),
		                 validTF = c(TRUE, FALSE)) # pos for validTF == T

	# now scale variance with error component
	pos.lm <- fit.lm %>% filter(validTF == T)

	# scale all signals based on LR
    # calculate scaled dVAR, aka transport error in ppm
	stat <- stat %>% mutate(scaled.var.err = pos.lm$s * var.orig + pos.lm$i, 
							scaled.dvar = scaled.var.err - var.orig, 
							scaled.dsd = sqrt(scaled.dvar), 
							sd.trans = replace(scaled.dsd, is.na(scaled.dsd), 0))

	list(stat, fit.lm)
}
