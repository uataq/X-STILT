
# if ODIAC or Vulcan, calc ENO based on non-EDGAR ECO2 and ERNO_EDGAR, DW, 03/03/2022
create_aq_non_edgar = function(eco2_non_edgar, eco2_edgar, eco_edgar, 
                               ech4_edgar, eno_edgar) {
    
    # ----------------------------------------------------------------------
    cat('create_aq_non_edgar(): ***ENOx is calc from ECO2 and ER from EDGAR...\n')
    
    # EDGAR based emission ratio
    erno_edgar = eno_edgar / eco2_edgar         # EDGAR based ER
    erco_edgar = eco_edgar / eco2_edgar 
    erch4_edgar = ech4_edgar / eco2_edgar 
    erno_prj = projectRaster(erno_edgar, eco2_non_edgar, method = 'ngb')
    erco_prj = projectRaster(erco_edgar, eco2_non_edgar, method = 'ngb')
    erch4_prj = projectRaster(erch4_edgar, eco2_non_edgar, method = 'ngb')

    eno_non_edgar = eco2_non_edgar * erno_prj  
    eco_non_edgar = eco2_non_edgar * erco_prj   
    ech4_non_edgar = eco2_non_edgar * erch4_prj   

    eno_non_edgar[eno_non_edgar < 0] = 0
    eco_non_edgar[eco_non_edgar < 0] = 0
    ech4_non_edgar[ech4_non_edgar < 0] = 0
    
    eno_rt = eno_non_edgar
    eco_rt = eco_non_edgar
    ech4_rt = ech4_non_edgar
    eco2_rt = eco2_non_edgar
    emis_stk = stack(eco2_rt, eco_rt, ech4_rt, eno_rt)
    names(emis_stk) = c('ECO2', 'ECO', 'ECH4', 'ENO')

    return(emis_stk)
}

