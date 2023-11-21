
load_radm2_vocr = function() {
    
    # for K = A exp(-E / (RT)) and VOCR = K_OH+VOCi * [VOCi]
    # exceptions: K = T^2 * A * exp(-E/(RT)) for ethane ETH and PAN
    # list of A and E/R shown below for NMVOC lumped groups in RADM2
    vocs = c('eth', 'hc3', 'hc5', 'hc8', 'ol2', 'olt', 'oli', 'iso', # C-C, C=C
             'tol', 'csl', 'xyl',               # Aromatic
             'hcho', 'ald', 'ket', 'gly', 'mgly', 'dcb',    # Carbonyls
             'pan', 'onit', 'op1', 'op2', 'paa')   # Org N and peroxides

    A = c(1.37E-17, 1.59E-11, 1.73E-11, 3.64E-11,  # C-C
          2.15E-12, 5.32E-12, 1.07E-11, 2.55E-11,  # C=C
          2.10E-12, 4.00E-11, 1.89E-11,                # aromatics
          9.00E-12, 6.87E-12, 1.20E-11, 1.15E-11, 1.7E-11, 2.8E-11, 
          # C=O, carbonyls
        
         # organic nitrogen and peroxides (PAN, TPAN (x), ANs, OP1, OP2, PAA)
         6.85E-18, 1.55E-11, 1E-11, 1E-11, 1E-11)

    E_R = c( 444, 540, 380, 380, -411, -504, -549, -409, # C-C, C=C
            -322, NA, -116, NA, -256, 745, NA, NA, NA, 444, 540, NA, NA, NA)
    K_df = data.frame(name = vocs, A = A, E_R = E_R, stringsAsFactors = FALSE) 

    return(K_df)
}

load_radm2_ro2_no = function() {
    
}