# X-STILT: Column-Stochastic Time-Inverted Lagrangian Transport model
## Descriptions:
X-STILT is an atmospheric transport model that deals with vertically integrated column CO2 and potentially other trace gases. Scripts and subroutines for X-STILT, modified from the Stochastic Time-Inverted Lagrangian Transport (STILT) model [*Lin et al*., 2003](https://doi.org/10.1029/2002JD003161) and STILT-R version 2 [*Fasoli et al*., 2018](https://doi.org/10.5194/gmd-11-2813-2018). Required datasets and methodologies are described in [*Wu et al*., 2018](https://doi.org/10.5194/gmd-11-4843-2018). 


We have merged column features with STILT-R version 2 (see github repo on https://uataq.github.io/stilt/) with improved model efficiency and more realistic "footprint" values. X-STILT model developments are ongoing towards a more flexible model framework that can be more easily applied to other column measurements. Contributions are welcomed and highly appreciated. Please contact Dien Wu (Dien.Wu@utah.edu) if you have any issues with the code. Thank you.

## Details for X-STILT based on STILT-R version 2 [*Fasoli et al*., 2018]:

Table of Contents:
- [**Prerequisites including input datasets**](#h1)
- [**Obtaining Column Footprint [ppm / (umol m-2 s-1)]**](#h2)
- [**Determining background XCO<sub>2</sub> value [ppm]**](#h3)
- [**Estimating horizontal and vertical transport errors**](#h4)
- [**Potential atmospheric inversion on XCO<sub>2</sub>**](#h5)
- [**Specific for your desired cities or space-based sensors**](#h6)
- [**Example Figures for an overpass over Riyadh**](#7)
- [**Reference**](#8)

<!-- toc -->
> ### **Prerequisites including input data** 
0. [Automatically install and load R packages](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L65-L68) required for X-STILT subroutines as stated in [dependencies](https://github.com/uataq/X-STILT/blob/master/r/dependencies.r). 

1. Download [OCO-2 Level 2 Lite files](https://disc.gsfc.nasa.gov/datasets/OCO2_L2_Lite_FP_9r/summary) and modify `oco2.path`; X-STILT will read in averaging kernels and pressure weighting functions from OCO-2 lite files to perform vertical weighting to footprint values for air parcels that releases from different altitudes and then provide the vertically compressed column footprints. See changes listed in ++ if using other column sensors. 

2. Download [1 km ODIAC files in tif format](http://db.cger.nies.go.jp/dataset/ODIAC/DL_odiac2019.html) and modify `tiff.path`. ODIAC is the main emission product we adopted to estimate anthropogenic XCO<sub>2</sub> concentrations in ppm. 

3. Meteorological fields in ARL format, e.g., the ones downloaded from ftp://arlftp.arlhq.noaa.gov/archives/, and modify [relevant lines in STEP 4](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L228-L234). 

4. Additional input data streams for performing transport error analyses:
   * [NOAA radiosonde data](https://ruc.noaa.gov/raobs/) for computing model-data wind errors; please choose wind speed unit of tenths of m s-1 and FSL format; Users could download all RAOB stations within a spatial area around their target city;
   * [Carbon-Tracker mole fraction data, e.g., CT-NRT](https://www.esrl.noaa.gov/gmd/ccgg/carbontracker/CT-NRT/) for getting the total XCO<sub>2</sub> in addition to FF XCO<sub>2</sub>. 

5. Additional input data streams for performing emission error analyses: 
   * Bottom-up emission inventories ensemble [FFDAS](http://ffdas.rc.nau.edu/index.html) and [EDGAR](https://edgar.jrc.ec.europa.eu/). 


> ### **Obtaining Column Footprint [ppm / (umol m<sup>-2</sup> s<sup>-1</sup>)]**

*Column Footprint* are the source-receptor sensivities or essentially the Jacobian Matrics between concentration (enhancements) and fluxes (for a given source/sink). Users can start with `run_xstilt.r` for model and parameter initializations.

1. Select one OCO-2 overpass. By defaul, [STEP 1](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L106-L111) searches for all overpasses that have soundings falling into a larger spatial domain (i.e., [2 deg x 3 deg](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L78)) around your city as well as a smaller urban domain (i.e., [1 deg x 1 deg](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L102-L104)). Thus, users can modify the spatial domain based on their city sizes. [This part of the code](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L123-L126) will generate maps of observed XCO<sub>2</sub> and SIF if ```plotTF == T```.

2. Indicate the kind of simulation users would like to start in [STEP 2](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L135-L166). One needs to modify logical flags which included:
   * `run_trajec = T or run_foot = T`: Backward trajectories + vertically weighted column footprints; 
   * `run_hor_err = T`: Horizontal transport error (one needs to calculate wind error statistics first); 
   * `run_ver_err = T`: Vertical transport error via scaling mixed layer height;
   * `run_sim = T`: Simulate FFCO2 XCO<sub>2</sub> enhancements (requires footprint), remember to turn off `run_trajec` and `run_foot`;
   * `run_emiss_err = T`: Prior emission error (requires footprint).

3. Modify parameters for placing column receptors (e.g., max height and vertical spacing between two vertical levels in meters, # of particles per level, receptor locations). By default, X-STILT will choose more receptors within the urban enhanced latitude range, see [STEP 3](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L201-L212).

4. Indicate the meteorological fields in [STEP 4](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L228-L234). 

5. Modify footprint information in [STEP 5](https://github.com/uataq/X-STILT/blob/master/run_xstilt.r#L272-L320), e.g., horizontal resolution and spatial extent. 

6. STEP 6 will finally run trajec or footprints via parallel computing.

> ### **Determining background XCO<sub>2</sub> value [ppm]**
Start with main script of `compute_bg_XCO<sub>2</sub>.r`. [M3 is the overpass specific background](https://github.com/uataq/X-STILT/blob/master/compute_bg_XCO<sub>2</sub>.r#L90-L181) by releasing air parcels in a forward fashion from a city and determining urban plume and background region using 2D kernel density. Please refer to details described in [Sect. 2.3 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section2).


> ### **Estimating horizontal and vertical transport errors**
Instead of using model ensembles for estimating errors, X-STILT propagates real-world random u-v- wind errors (with correlation length scale and timescale) into errors in XCO<sub>2</sub>. The fundamental approach is proposed and documented in [Lin and Gerbig, 2005](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2004GL021127), which is now extended to column CO2. We calculate and propagate XCO<sub>2</sub> errors from release levels, to receptor locations, and finally to overpasses. 

Briefly speaking, users need to first generate another set of trajectories with wind error statistics by turning on `run_hor_err`. Once trajectories with wind errors have been generated, one need to turn on `run_sim` to let X-STILT calculate the errors in XCO<sub>2</sub> based on two sets of trajectories (with/without wind errors). Hard to explain all the steps here, but feel free to contact Dien (dien.wu@utah.edu) if you’re interested in carrying out error estimates. Please also refer to details described in [Sect. 2.6 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section2).


> ### **Potential atmospheric inversion on XCO<sub>2</sub>**
As discussed in [Sect. 4.2 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section4), we could come up the posterior scaling factor for anthropogenic emissions based on 5 overpasses over Riyadh, via a simple Bayesian Inversion. We treated the entire city as a whole and solve for one scaling factor, given biases in near-field wind direction. So, we did not solve for posterior emissions for every grid cell within a city. Codes are not included in this repo. 


> ### **Specific for your desired cities or space-based sensors**
As discussed in [Sect. 4.3 in Wu et al. (2018)](https://www.geosci-model-dev.net/11/4843/2018/#section4), X-STILT can be applied to other column measurements and other species. One could simply change `site` in the main scripts for modeling other cities besides Riyadh. 

If the user would like to use X-STILT for other column measurements besides OCO-2 (e.g., TCCON, GOSAT, or TROPOMI), one needs to create a new subroutine for grabbing the averaging kernel profile (`ap`) and pressure weighting functions (`pwf`) for every pressure level (`pres`) as well as surface pressure (`psfc`) for each receptor location (indicated by `find.lat` and `find.lon`). The current subroutine for grabbing OCO-2 info is [get.oco2.info.r](https://github.com/uataq/X-STILT/blob/master/r/src/get.oco2.info.r). Please form your results in the same fashion as `all.info` and keep the list names as before to not break any subroutines downstream:
   ```
   all.info <- list(oco2.id = find.id, oco2.lat = find.lat,
                    oco2.lon = find.lon, ak.norm = ak.norm, 
                    pwf = pwf, pres = pres, ap = ap, 
                    oco2.grdhgt = grdhgt, oco2.psfc = psfc, 
                    oco2.foot = footprint, oco2.xco2 = XCO2, 
                    oco2.xco2.uncert = xco2.uncert)
   ```
   We are also working towards a more flexible version of X-STILT that could potentially work with other sensors. 


> ### **Example Figures of column footprints and XCO<sub>2.ff</sub> for an overpass over Riyadh**
![](wgt_sum_xfoot_Riyadh_2015121610_gdas0p5_STILTv2_zoom8_-72hrs_100dpar.png)
Figure 1: Latitude integrated map of weighted column footprints [umol/m2/s] on 12/29/2014 from 70+ selected sounding/receptor over Riyadh.

![](wgt_sum_xco2_Riyadh_2015121610_gdas0p5_STILTv2_zoom8_-72hrs_100dpar.png)
Figure 2: Latitude integrated XCO<sub>2</sub>.ff contribution maps [ppm] on 12/29/2014 from 70+ selected sounding/receptor over Riyadh.


>### **Reference**
Wu, D., Lin, J. C., Fasoli, B., Oda, T., Ye, X., Lauvaux, T., Yang, E. G., and Kort, E. A.: A Lagrangian approach towards extracting signals of urban CO2 emissions from satellite observations of atmospheric column CO2 (XCO<sub>2</sub>): X-Stochastic Time-Inverted Lagrangian Transport model (“X-STILT v1”), *Geosci. Model Dev.*, 11, 4843-4871, https://doi.org/10.5194/gmd-11-4843-2018, 2018. 

Fasoli, B., Lin, J. C., Bowling, D. R., Mitchell, L., and Mendoza, D.: Simulating atmospheric tracer concentrations for spatially distributed receptors: updates to the Stochastic Time-Inverted Lagrangian Transport model's R interface (STILT-R version 2), *Geosci. Model Dev.*, 11, 2813-2824, https://doi.org/10.5194/gmd-11-2813-2018, 2018.

Lin, J.C., Gerbig, C., Wofsy, S.C., Andrews, A.E., Daube, B.C., Davis, K.J. and Grainger, C.A.: A near‐field tool for simulating the upstream influence of atmospheric observations: The Stochastic Time‐Inverted Lagrangian Transport (STILT) model. *Journal of Geophysical Research: Atmospheres*, 108(D16), https://doi.org/10.1029/2002JD003161. 2005. 

Lin, J. C., and Gerbig, C., Accounting for the effect of transport errors on tracer inversions, *Geophys. Res. Lett.*, 32, L01802, https://doi.org/10.1029/2004GL021127, 2015.