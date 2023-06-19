# subroutine to look for the correct file of chemical transport model 
# mimic find_met_file(), edited by DW 

find_ctm_files <- function(r_run_time, n_hours, ctm_file_format, ctm_path) {
    
    is_backward <- n_hours < 0

    request <- as.POSIXct(r_run_time, tz = 'UTC') %>%
        c(. + c(n_hours, is_backward * n_hours - 5) * 3600) %>%
        range() %>%
        (function(x) seq(x[1], x[2], by = 'hour')) %>%
        strftime(tz = 'UTC', format = ctm_file_format) %>% unique()
    
    available <- dir(ctm_path, full.names = T)
    available <- available[!grepl('\\.lock$', available)]
    
    idx <- do.call(c, lapply(request, function(pattern) {
                             grep(pattern = pattern, x = available) }))

    if (any(idx < 1)) return()
    
    fn <- unique(available[idx])
    if (length(fn) == 0) return() 

    if (length(fn) > 0) 
        cat('find_ctm_files(): Found CTM files.\n')
        return(fn)
}
