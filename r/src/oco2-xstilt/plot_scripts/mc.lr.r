# subroutine to carry out Monte Cario tests for fitting a linear regression 
# based on modeled vs. observed values with their uncertainties
# assuming normal distribution, DW, 09/17/2018

mc.lr <- function(x, y, x.sd, y.sd, fac, boottime = 5000) {

    cat('\n\nStart Monte Carlo tests for linear regression...\n')
    library(lmodel2)

    # fit linear regression for Orthogonal Nonlinear Least-Squares Regression
    count <- 100

    melt.norm <- NULL
    for (n in 1:length(x)) {
      x.norm <- rnorm(n = count, mean = x[n], sd = x.sd[n])
      y.norm <- rnorm(n = count, mean = y[n], sd = y.sd[n])
      tmp.norm <- data.frame(x = x.norm, y = y.norm, fac = rep(fac[n], count), 
        n = rep(n, count), count = seq(1, count))
      melt.norm <- rbind(melt.norm, tmp.norm)
    }

    if (F) {
      h = hist(x.norm[1,], plot = F)
      xfit <- seq(min(x.norm[1,]), max(x.norm[1,]), length = 40)
      yfit <- dnorm(xfit,mean = x[1], sd = x.sd[1])
      yfit <- yfit * diff(h$mids[1:2]) * length(x.norm[1,])
      hist(x.norm[1,], ylim = c(0,30))
      lines(xfit, yfit, col = 'blue', lwd = 2)
    }

    boot <- function(x){return(sample(x, size = 1))}

    slope.all <- NULL; int.all <- NULL; cor.all <- NULL
    for (c in 1:boottime) {

      if(c == 1 || c %% 1000 == 0)
        cat(paste(c / boottime * 100, '% Monte Carlo for linear regression...\n'))

      sel.norm <- data.frame(x = tapply(melt.norm$x, melt.norm$n, boot), 
                             y = tapply(melt.norm$y, melt.norm$n, boot))

      fit.lm <- suppressMessages(lmodel2(sel.norm$y ~ sel.norm$x)$regression.results)
      tmp.lm <- fit.lm[fit.lm$Method == 'SMA', ]
      slope.all <- c(slope.all, tmp.lm$Slope)
      int.all   <- c(int.all, tmp.lm$Intercept)
      cor.all   <- c(cor.all, cor(sel.norm$y, sel.norm$x))
    }

    stat.all <- data.frame(slope = slope.all, int = int.all, cor = cor.all)

    cat(paste('mean of r:', signif(mean(stat.all$cor), 3),
              ';SD of r:', signif(sd(stat.all$cor), 3), '\n\n'))
    cat(paste('mean of slope:', signif(mean(stat.all$slope), 3),
              ';sd of slope:', signif(sd(stat.all$slope), 3), '\n\n'))

    return(stat.all)
}