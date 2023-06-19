
load_epa_pp = function(epa_name, pp_fn = NULL) {
    if (is.null(pp_fn)) pp_fn = '/central/home/dienwu/postdoc_proj/NOx/pp/quarterly-emission-comparison-2020-vs-2021.xlsx'

    loc_df = readxl::read_excel(pp_fn, sheet = 3, skip = 3) %>% 
             filter(grepl(epa_name, `Facility Name`)) %>% 
             dplyr::select(lon = `Longitude (degrees)`, 
                           lat = `Latitude (degrees)`) %>% unique()
    
    if ( nrow(loc_df) == 0 ) loc_df = NULL
    return(loc_df)
}
