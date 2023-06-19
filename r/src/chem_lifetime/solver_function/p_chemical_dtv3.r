

p_chemical_dtv3 = function(p_df, tau_nox_df, delt, num_chunk = 1) {
    
    # ------------------------------------
    # calculate NET changes in C, i.e., L - P = delt * [nox] / tau
    # positive for NET LOSS, p_noxa = p_noxb * (1 - delt / tau_net)
    # dchem is positve for net LOSS; is negative for net generation

    # *** need to check abs(ts_nox) [hr] versus delt [usually 1min]
    # if ts_nox is much smaller than delt = 1min, need to exp()
    # if ts_nox is much larger than delt = 1min, can do approx
    p_mini = p_df %>% 
             mutate(dchem = delt / ts_nox / 60 * tmp_nox, 
             
                    # since I defined ts_nox > 0 for net loss
                    # dchem > 0 for net production, < 0 for net loss
                    #' now directly modify @column tmp_nox
                    tmp_nox = tmp_nox - dchem + demiss_no / num_chunk + tmp_dentrain,       # [nox] at current time stamp now

                    # calculate [NO2], DW, 02/08/2022
                    # estimate OX at the current time, = OX_bg + slope * NOX
                    # NOx needs to be in ppb
                    #tmp_OX = slp * tmp_nox * 1e3 + int, 
                    #tmp_OX = ifelse(tmp_OX > 95, 95, tmp_OX), 
                    tmp_OX = 50, 

                    ### TWO WAYS TO ESTIMATE NO2 using [NO2] / [NO] = K_OH+O3 * [O3] / J(NO2), where O3 is prescribed as [OX] - [NO2]
                    # BUT either way, we need to calc J and K

                    # convert K from cm3 molecules-1 s-1 to ppm-1 s-1
                    # DW, 03/28/2022 -------------------------------------------
                    # requires air density (from STILT in kg m-3)
                    # * 1e-6 * 1e3 convert density from kg m-3 to g cm-3
                    # * 28.97 g mol-1 -> mol cm-3 air
                    # * 6.023e23 -> molecules cm-3 air 
                    # * 1e-9 for converting unitless to ppb
                    # thus the final unit of K is ppb-1 s-1 now
                    unit_conv = dens * 1E-6 * 1E3 / 28.97 * 6.023E23 * 1E-9, 
                    K = 2.07E-12 * exp(-1400 / temz) * unit_conv,   # ppb-1 s-1
                    J = 0.0167 * exp(-0.575 / cos(psza * pi / 180)), 
                    J = ifelse(psza > 85, 1e-16, J), 
                    #J = 0.0111 * exp(-58 / temp),      # in s-1

                    # Method 1 -----------------------------------------------
                    # the one [NO2] for [O3] is obtained from last time step
                    # calc ratio based on NO2 from the last time stamp, 
                    # i.e., r = 1 / ( 1 + [NO] / [NO2])
                    # before updating tmp_no2, this is the [NO2] from last time
                    #tmp_rto = 1 / (1 + J / (K * (tmp_OX - tmp_no2))),

                    # then update NO2 based on updated ratio and current NOx
                    #tmp_no2 = tmp_rto * tmp_nox, 
                    # END OF Method 1 ------------------------------------------

                    # Method 2 -----------------------------------------------
                    # solving quadratic solution for [NO2]
                    # J / k in [s-1 * ppb * s] = ppb
                    # tmp_nox in ppm, convert to ppb
                    b = -J / K - tmp_nox * 1e3 - tmp_OX,            # b [ppb]
                    c = tmp_nox * 1e3 * tmp_OX,                     # c [ppb^2]
                    
                    # When a, b, and c are integers, the solutions to ax2 + bx + c = 0 are rational numbers if and only if b2 - 4ac is a perfect square.
                    # the quadratic solution has unit of ppb-NO2
                    # overwrite [NO2] for the current time stamp in ppm now 
                    tmp_no2 = (-b - sqrt(b^2 - 4 * c)) / 2 / 1e3,  # back to ppm
                    tmp_rto = tmp_no2 / tmp_nox,        # ppm / ppm for ratio
                    # END OF Method 2 ------------------------------------------

                    # several checks 
                    tmp_rto = ifelse(tmp_rto > 1, 1, tmp_rto), 
                    tmp_rto = ifelse(tmp_rto < 0, 0, tmp_rto)
                    
                ) %>% dplyr::select(-c(K, J, b, c, int, slp, unit_conv))

    return(p_mini)

}
