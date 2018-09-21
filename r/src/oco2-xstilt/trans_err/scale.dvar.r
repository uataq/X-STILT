## subroutine to fit linear regression and get ggplots for linear regression plots
# outputs are linear regression line, CO2 statiatics, and ggplot
# written by DW, 03/15/2018

# generalize and clear up codes, DW, 07/23/2018
# originate from 'get.ggplot.lr.r' in my stiltR,
# original methods c("all","1.5max","10ppm","wgt"), only use wgt

scale.dvar <- function(co2.stat){

	# M4, properly weighted by variances of error
	# select positive/negative error variances diff
	stat.pos <- co2.stat[co2.stat$dVAR > 0, ]
	stat.neg <- co2.stat[co2.stat$dVAR < 0, ]

	wgt.pos <- 1/stat.pos$var.err  # positive err var
	wgt.neg <- 1/stat.neg$var.err

  # fit two weighted linear regressions
	fit.pos.lm <- suppressMessages(
		lm(stat.pos[, "var.err"] ~ stat.pos[, "var.orig"], weights = wgt.pos))

	fit.lm.neg <- suppressMessages(
		lm(stat.neg[, "var.err"] ~ stat.neg[, "var.orig"], weights = wgt.neg))

	fit.lm <- data.frame(
		s = as.numeric(c(fit.lm.true$coefficients[2], fit.lm.false$coefficients[2])),
		i = as.numeric(c(fit.lm.true$coefficients[1],fit.lm.false$coefficients[1])),
		validTF = c(TRUE, FALSE)
	)

	# now scale variance with error component
	pos.lm <- fit.lm %>% filter(validTF == T)

	# scale all signals based on LR
	co2.stat$scaled.var.err <- pos.lm$s * co2.stat$var.orig + pos.lm$i

	# calculate scaled dVAR, aka transport error in ppm
	co2.stat$scaled.dvar <- co2.stat$scaled.var.err - co2.stat$var.orig

	return(list(co2.stat, pos.lm))
}
