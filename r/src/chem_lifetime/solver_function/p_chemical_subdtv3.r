
# calculate [NOx] at sub-minute time scale
p_chemical_subdtv3 = function(p_subdt, p_before, delt = 1, tau_nox_df, OX_df, 
                              tau_chem = c('v20230503', 4, NA)[3]) {
    
    # default chucks of time stamp at sub-minute scale in second
    dt_bins = c(1, 2, 3, 4, 5, 6, 10, 12, 15, 20)    

    cat(paste('p_chemical_subscale(): reducing chem timestamp for', 
              nrow(p_subdt), 'particles\n'))

    # need to reduce delt, generate dynamic delt that << ts_nox * 60 
    # (at least 1/4 of ts_nox), with ts_* in hours
    dtt = abs(p_subdt$ts_nox) * 60 * 60 / 4                # in second
    dtt_crp = dt_bins[findInterval(dtt, dt_bins)]        # in second
    nloop = delt * 60 / dtt_crp     # # of loop for sub-min calc
    nindx = unique(p_subdt$indx)               # all indx # that needs fix
    
    # -------------------------------------------------------------------- #
    # find the air parcels whose delt >> net loss time scale
    p_corr = NULL
    for ( i in 1 : length(nindx) ) {

        tmp_loop = p_before %>% filter(indx == nindx[i]) 
        tmp_chem_total = 0 

        for ( l in 1 : nloop[i] ) {     # nloop for number of sub-min chucks
            
            # since tmp_nox and tmp_rto has been updated for the last iteration, we can remove some additional columns to avoid duplicated column names for the current iteration
            if (l > 1)  
                tmp_loop = tmp_loop[colnames(tmp_loop) %in% colnames(p_before)]

            # check out NOx NET change frequency again based on the newly updated [NOX] for each iteration         
            tmp_loop = lookup_nox_tendency(tmp_loop, tau_nox_df, OX_df, 
                                           delt = dtt[i] / 60, tau_chem)

            tmp_loop = p_chemical_dtv3(p_df = tmp_loop, tau_nox_df, 
                                       delt = dtt[i] / 60, num_chunk = nloop[i])

            # need to modify the dchem column to sum up all dchem from all iterations
            tmp_chem_total = tmp_chem_total + tmp_loop$dchem
            #print(i); print(tmp_chem_total)
        }   # end for l

        tmp_loop[tmp_loop$indx == nindx[i], 'dchem'] = tmp_chem_total
        p_corr = rbind(p_corr, tmp_loop)
    }   # end for i

    return(p_corr)
}

