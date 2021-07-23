# X-STILT: Column-Stochastic Time-Inverted Lagrangian Transport model

X-STILT is an atmospheric transport model that deals with vertically integrated column CO2 and potentially other trace gases. Subroutines for X-STILT are modified from the Stochastic Time-Inverted Lagrangian Transport (STILT) model [*Lin et al*., 2003](https://doi.org/10.1029/2002JD003161) and STILT-R version 2 [*Fasoli et al*., 2018](https://doi.org/10.5194/gmd-11-2813-2018). 

Since STILT-R version 2 now serves as a submodule of X-STILT repo, you will need to *"git clone --recursive"* the entire package with submodule. Required datasets and methodologies are described in [prerequisites](prerequisites) and [*Wu et al*., 2018](https://doi.org/10.5194/gmd-11-4843-2018). Please refer to the [reference](#reference) section if you plan to publish any study using the (X-)STILT model. 

Currently, X-STILT can work with OCO-2 and 3, as well as TROPOMI CO, CH4, and NO2. Model developments are ongoing and contributions are welcomed and appreciated. Please contact Dien (dienwu@caltech.edu) if you are interested in other column sensors/species or have any questions with the code. Thank you.

# X-STILT Features
## Table of Contents
- [**Latest commit**](#latest-commit-by-july-06-2021)
- [**Prerequisites**](#prerequisites)
- [**Obtain column footprint**](#obtain-column-footprint)
- [**Determine background XCO<sub>2</sub>**](#determine-background-xco2)
- [**Estimate horizontal and vertical transport errors**](#estimate-horizontal-and-vertical-transport-errors)
- [**Atmospheric inversion on XCO<sub>2</sub>**](#atmospheric-inversion-on-xco2)
- [**Generalization**](#generalization)
- [**Example figures of column footprints and XCO<sub>2.ff</sub>**](#example-figures-of-column-footprints-and-xco2ff)
- [**Reference**](#reference)


Latest commit by July 6, 2021
============
Changes to ['merge_stilt-hysplit' branch](https://github.com/uataq/X-STILT/tree/merge_stilt-hysplit) towards the release of X-STILT version 2 that can work with non-OCO sensors or without the dependence of satellite data:

:pushpin: X-STILT can now work with TROPOMI X-gas data, such as 
   * footprint weighting using TROPOMI species-specific averaging kernels;

   * functions to grab vertical column density and convert to vertical mixing ratio in ppb; 

   * need to specify TROPOMI info via `obs_sensor` and `obs_species` that are currently limited to CO, CH4, and NO2;

:pushpin: IDEAL simulations for column footprints/enhancements from any locations (no requirement of satellite data).
   * To use the ideal mode, `obs_sensor` needs to be set as NA; 
   * Column AK is simply treated as 1 from TOA to the surface (`ak_wgt` is forced as FALSE);
   * Column PWF is calculated using modeled variables/profiles (`pwf_wgt` = TRUE). 

   * Must provide a .csv/.txt file containing the latitudes (`lati`) and longitudes (`long`) of your desired receptor locations. Receptor time (in form of YYYYMMDDHH or YYYYMMDDHHss) can either be included in the csv/txt file as the 3rd column (with name of `time`) OR assigned to `timestr` in `run_xstilt.r`. See example of `receptor_demo.csv`. 
   
   * [USE WITH CAUTIONS] NOT suitable for comparisons with real satellite sensors that have variable vertical sensitivities.

:pushpin: **Built upon the latest HYSPLITv5**: Replace discrete release levels with *line source release* -- i.e., particles are now evenly distributed between `minagl` and `maxagl`, e.g., default is 3000 particles from 0 - 3km; 

:pushpin: PWF calculations in `r/src/recp_trajec_foot`

   * Instead of interpolate PWF from OCO-2 to each discrete release level :arrow_forward: calculate PWF based on dp and dry air column density (mol m-2) between particles. 

   * To calculate the dry air X-density, we extracted specific humidity and temperature profiles via `get.met.vars()`. 

   * If ak.wgt is turned off (`ak.wgt` == FALSE), the dependence on any satellite sensors will be removed. That being said, normalied AK is 1 and PWF is calculated based on MODELED surface pressures instead of RETRIEVED surface pressures. 

:pushpin: Other minor changes: 
   * Allow the calculation of another set of footprint with a different temporal resolution from the default run at the same time, e.g., hourly footprint implied by `time_integrate2 = F`;

   * Add interactive mode for choosing site, overpass time, etc. One can modify the code to turn it off;

   * Insert the google API as the first line of `insert_ggAPI.csv` that will be sourced by the model for geolocation and google map plotting; 

   * Get rid of the additional .rds files end with X_wgttraj.rds to reduce file sizes :arrow_forward: trajectory-level footprints after the AK and PWF weighting are now stored as additional columns in the initial *X_traj.rds files. Here are the descriptions of columns in output$particle from rds files: 
      * `foot_before_weight:` footprint BEFORE AK and PWF weighting;
      * `foot:` footprint AFTER AK and PWF weighting; 
      * `pwf, ak.norm, ap, ak.pwf:` computed pressure weighting and interpolated apriori and averaging kernel profiles from OCO-2/3; 
      * `wgt:` final weighting factors used for scaling trajectory-level footprint based on flags of `ak.wgt` and `pwf.wgt`. 


<!-- toc -->
Prerequisites
============
1. Download [OCO-2 Level 2 Lite files](https://disc.gsfc.nasa.gov/datasets/OCO2_L2_Lite_FP_9r/summary) and modify `oco2.path`; X-STILT will read in averaging kernels and pressure weighting functions from OCO-2 lite files to perform vertical weighting to footprint values for air parcels that releases from different altitudes and then provide the vertically compressed column footprints. See [this section](#generalization) if you desired to use other column sensors. 

2. Download [1 km ODIAC files in tif format](http://db.cger.nies.go.jp/dataset/ODIAC/DL_odiac2019.html) and modify `tiff.path`. ODIAC is the main emission product we adopted to estimate anthropogenic XCO<sub>2</sub> concentrations in ppm. 

3. Meteorological fields in ARL format, e.g., ones provided by [NOAA ARL](ftp://arlftp.arlhq.noaa.gov/archives/), and modify the relevant lines in [STEP 4](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L228-L234). 

4. Additional input data streams for performing transport error analyses:
   * [NOAA radiosonde data](https://ruc.noaa.gov/raobs/) for computing model-data wind errors; please choose wind speed unit of tenths of m s-1 and FSL format; Users could download all RAOB stations within a spatial area around their target city;
   * [Carbon-Tracker mole fraction data, e.g., CT-NRT](https://www.esrl.noaa.gov/gmd/ccgg/carbontracker/CT-NRT/) for getting the total XCO<sub>2</sub> in addition to FF XCO<sub>2</sub>. 

5. Additional input data streams for performing emission error analyses: 
   * Bottom-up emission inventories ensemble [FFDAS](http://ffdas.rc.nau.edu/index.html) and [EDGAR](https://edgar.jrc.ec.europa.eu/). 


Obtain column footprint
============
*Column Footprints* [ppm / (umol m<sup>-2</sup> s<sup>-1</sup>)] are the source-receptor sensivities or essentially the Jacobian Matrics between concentration (enhancements) and fluxes (for a given source/sink). Users can start with `run_xstilt.r` for model and parameter initializations.

1. Select a satellite overpass. By default, [STEP 1](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L106-L111) searches for all overpasses that have soundings falling into a larger spatial domain (i.e., [2 deg x 3 deg](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L78)) around your city as well as a smaller urban domain (i.e., [1 deg x 1 deg](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L102-L104)). Thus, users can modify the spatial domain based on their city sizes. [This part of the code](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L123-L126) will generate maps of observed XCO<sub>2</sub> and SIF if ```plotTF == T```.

2. Indicate the kind of simulation users would like to carry out in [STEP 2](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L135-L166). One needs to modify logical flags which included:

   >* `run_trajec = T or run_foot = T`: Backward trajectories + vertically weighted column footprints; 

   >* `run_hor_err = T`: Horizontal transport error analysis (one needs to calculate wind error statistics first, e.g., by comparing modeled winds with observed winds from NOAA RAOB, see [STEP 4](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L237-L241) for more info); 

   >* `run_ver_err = T`: Vertical transport error analysis (via scaling mixed layer height up and down with scaling factors stated in [STEP 4](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L253-L256));

   >* `run_sim = T`: Simulate FFCO2 XCO<sub>2</sub> enhancements (requires footprint), remember to turn off `run_trajec` and `run_foot`;

   >* `run_emiss_err = T`: Prior emission error analysis (requires footprint).

3. Modify parameters for placing column receptors (e.g., max height and vertical spacing between two vertical levels in meters, # of particles per level, receptor locations). By default, X-STILT will choose more receptors within the urban enhanced latitude range, see [STEP 3](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L201-L212).

4. Indicate the meteorological fields in [STEP 4](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L228-L234). 

5. Modify footprint information in [STEP 5](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L272-L320), e.g., horizontal resolution and spatial extent. 

6. Start running trajec or footprints via parallel computing in [STEP 6](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L322-L360).



Determine background XCO<sub>2</sub>
============
Start with main script of `compute_bg.r`. [M3 is the overpass specific background [ppm]](https://github.com/uataq/X-STILT/blob/master/compute_bg_XCO<sub>2</sub>.r#L90-L181) by releasing air parcels in a forward fashion from a city and determining urban plume and background region using 2D kernel density. Please refer to details described in [Sect. 2.3 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section2).



Estimate horizontal and vertical transport errors
============
Instead of using model ensembles for estimating errors, X-STILT propagates real-world random u-v- wind errors (with correlation length scale and timescale) into errors in XCO<sub>2</sub>. The fundamental approach is proposed and documented in [Lin and Gerbig, 2005](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2004GL021127), which is now extended to column CO2. We calculate and propagate XCO<sub>2</sub> errors from release levels, to receptor locations, and finally to overpasses. 

Briefly speaking, users need to first generate another set of trajectories with wind error statistics by turning on `run_hor_err`. Once trajectories with wind errors have been generated, one need to turn on `run_sim` to let X-STILT calculate the errors in XCO<sub>2</sub> based on two sets of trajectories (with/without wind errors). Hard to explain all the steps here, but feel free to contact Dien (dien.wu@utah.edu) if you’re interested in carrying out error estimates. Please also refer to details described in [Sect. 2.6 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section2).


Atmospheric inversion on XCO<sub>2</sub>
============
As discussed in [Sect. 4.2 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section4), we could come up the posterior scaling factor for anthropogenic emissions based on 5 overpasses over Riyadh, via a simple Bayesian Inversion. We treated the entire city as a whole and solve for one scaling factor, given biases in near-field wind direction. So, we did not solve for posterior emissions for every grid cell within a city Codes are not included in this repo. 


Generalization
============
As discussed in [Sect. 4.3 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section4), X-STILT can be applied to other column measurements and other species. One could simply change `site` and `obs_sensor` in the main scripts for modeling other cities besides Riyadh. Currently, the model can work with OCO-2 and 3, as well as TROPOMI CO, CH4, and NO2. Let us know if you are interested in other sensors/species. 


Example figures of column footprints and XCO<sub>2.ff</sub>
============
![](wgt_sum_xfoot_Riyadh_2015121610_gdas0p5_STILTv2_zoom8_-72hrs_100dpar.png)
Figure 1: Latitude integrated map of weighted column footprints [umol/m2/s] on 12/29/2014 from 70+ selected sounding/receptor over Riyadh.

![](wgt_sum_xco2_Riyadh_2015121610_gdas0p5_STILTv2_zoom8_-72hrs_100dpar.png)
Figure 2: Latitude integrated XCO<sub>2</sub>.ff contribution maps [ppm] on 12/29/2014 from 70+ selected sounding/receptor over Riyadh.


Reference
============
Wu, D., Lin, J. C., Fasoli, B., Oda, T., Ye, X., Lauvaux, T., Yang, E. G., and Kort, E. A.: A Lagrangian approach towards extracting signals of urban CO2 emissions from satellite observations of atmospheric column CO2 (XCO<sub>2</sub>): X-Stochastic Time-Inverted Lagrangian Transport model (“X-STILT v1”), *Geosci. Model Dev.*, 11, 4843-4871, https://doi.org/10.5194/gmd-11-4843-2018, 2018. 

Fasoli, B., Lin, J. C., Bowling, D. R., Mitchell, L., and Mendoza, D.: Simulating atmospheric tracer concentrations for spatially distributed receptors: updates to the Stochastic Time-Inverted Lagrangian Transport model's R interface (STILT-R version 2), *Geosci. Model Dev.*, 11, 2813-2824, https://doi.org/10.5194/gmd-11-2813-2018, 2018.

Lin, J.C., Gerbig, C., Wofsy, S.C., Andrews, A.E., Daube, B.C., Davis, K.J. and Grainger, C.A.: A near‐field tool for simulating the upstream influence of atmospheric observations: The Stochastic Time‐Inverted Lagrangian Transport (STILT) model. *Journal of Geophysical Research: Atmospheres*, 108(D16), https://doi.org/10.1029/2002JD003161. 2003. 

Lin, J. C., and Gerbig, C., Accounting for the effect of transport errors on tracer inversions, *Geophys. Res. Lett.*, 32, L01802, https://doi.org/10.1029/2004GL021127, 2005.