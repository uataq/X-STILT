# subroutine to split string and convert to data frame, DW, 12/20/2018

strsplit.to.df <- function(string, sep = '_', ncol = NULL, colnms = NULL) {

    library(stringr)
    ncol <- unique(str_count(string, sep) + 1)
    info <- matrix(unlist(strsplit(string, sep)), ncol = ncol, byrow = T)
    info.df <- data.frame(info, stringsAsFactors = F)

    if (is.null(colnms)) colnms <- paste0('V', seq(1, ncol, 1))
    colnames(info.df) <- colnms
    return(info.df)
}
