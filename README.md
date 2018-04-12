# X-STILT
## Descriptions:
Main R functions for X-STILT, modified from the Stochastic Time-Inverted Lagrangian Transport (STILT) model [*Lin et al*., 2003] along with X-STILT error analysis, based on *Wu et al.* [in prep].

## Details:
- XCO2 modeling (./src/xco2_modeling, topmost subroutine is oco2.get.xco2() from "OCO2.get.xco2.r")
1. Generate column trajectories
2. Weight trajec-level footprint by satellite profiles (i.e., prior, averaging kernels and pressure weighting)
3. Generate weighted column footprint
4. Convolve footprint from 3. with different emission grids

- XCO2 error analysis (./src/xco2_uncert)
1. Estimate XCO2 errors due to column transport errors
2. Estimate XCO2 errors due to emission errors
