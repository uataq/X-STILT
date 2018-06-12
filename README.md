# X-STILT: X-Stochastic Time-Inverted Lagrangian Transport model
## Descriptions:
Main R scripts and subroutines for X-STILT, modified from the Stochastic Time-Inverted Lagrangian Transport (STILT) model [*Lin et al*., 2003] based on the manuscript of *Wu et al.* (submitted to *Geoscientific Model Development*).

## Details:
- X-STILT runs (column or fixed receptors)
1. Change parameters in "XSTILT_run_trajec.r" and store the namelist
2. Automatically call subroutines to generate backward or forward trajectories
3. Allow for distributing receptors/soundings to several STILT copies and multiple program calculations (update on 05/10/2018)

- X-STILT XCO2 simulations (incorporate satellite profiles)
1. Change parameters in "XSTILT_sim_xco2.r" and store the namelist
2. Weight trajec-level footprint by satellite profiles (i.e., prior, averaging kernels and pressure weighting)
3. Generate/store weighted column footprint
4. Convolve column footprint (from above step) with different emission/flux grids

Model developments are ongoing. Please contact Dien if you have any questions/comments/suggestions. Thank you.
