# X-STILT: X-Stochastic Time-Inverted Lagrangian Transport model 
## Descriptions:
Main R scripts and subroutines for X-STILT, modified from the Stochastic Time-Inverted Lagrangian Transport (STILT) model [*Lin et al*., 2003] along with X-STILT error analysis, based on *Wu et al.* [submitted].

## Details:
- X-STILT runs (column or fixed receptors)
1. Change parameters in "XSTILT_run_trajec.r" and store the namelist
2. Automatically call ./src/run.backward.trajec() OR ./src/run.forward.trajec()
   to generate backward or forward trajectories

- X-STILT XCO2 simulations (incorporate satellite profiles)
1. Change parameters in "XSTILT_sim_xco2.r" and store the namelist
2. Weight trajec-level footprint by satellite profiles (i.e., prior, averaging kernels and pressure weighting)
3. Generate/store weighted column footprint
4. Convolve column footprint (from above step) with different emission/flux grids

- XCO2 error analysis (./src)
1. Estimate XCO2 errors due to column transport errors
2. Estimate XCO2 errors due to emission errors
