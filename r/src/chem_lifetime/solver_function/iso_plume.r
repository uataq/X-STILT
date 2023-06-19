iso_plume = function(var, pval = 0.1) {
     f.var = ecdf(var)
     f.inv = GoFKernel::inverse(f.var, lower = min(var, na.rm = T), 
                                        upper = max(var, na.rm = T))
     tno2_td = f.inv(1 - pval)
     zscore = (var - mean(var, na.rm = T)) / sd(var, na.rm = T)
     plumeTF = zscore >= tno2_td
     return(plumeTF)
}
