# Predicting peak aphid densities using GAMs

## Training
In this repository we use Bayesian GAMs to determine when Aphid densities peak in field data collected using pan traps and net sweeps.

#### Data
The [original paper](https://resjournals.onlinelibrary.wiley.com/doi/full/10.1111/afe.12564) and [raw data](https://osf.io/j5vwy/) are open access. Idaho pan trap data were used for the pan trap model, and Vetch 2019 & 2020 for sweeps, both found in `./aphids_2023/osf/Aggregated Aphid Data`.

#### Scripts
Aphid densities were modeled against Aphid degree days with smooths, pooling all sites and years. More complex models were considered (including site and year random effects, and elevation as a fixed covariate), but either degraded model fit or required longer computation for minimal change in predicted peaks.

Code used for final GAM fitting is in `./aphids_2023/R/Pan Trap.R`, with some messier exploratory work in `./aphids_2023/R/ETL.R`.

## Extrapolating
Median cumulative degree days at the aphid peak were extracted from the fitted models. These CDDs were back-transformed to dates at new field sites by scraping climate data for those sites and modeling CDD as a function of Julian days.

#### Data
Predictions were made for each site in `./aphids_2023/osf/Aggregated Aphid Data/Vetch 2019 and 2020 with degree days.csv`. Climate data from 2000 to 2020 were scraped for each site from daymet.

#### Scripts
Code can be found in `./aphids_2023/R/Fetch Site DD.R`, which sources code from `./aphids_2023/R/cdd_density_functions.R` (originally from the phenogamML repo) to pull data with daymetr.

#### Outputs
Predictions generated by pan trap data are found in `./aphids_2023/outputs/peak_predictions_2023_pan_traps.csv`, and predictions generated by sweeps data are found in `./aphids_2023/outputs/peak_predictions_2023_pan_traps.csv`.

The point estimates for Julian day at peak aphid density as a function of CDD (median, lower, and upper estimates) are given by julian_peak_lower, _median, and _upper. These values are expressed at dates in 2023 in date_lower, _median, and _upper. For each point estimate, a 99% confidence interval is specified in days (ci_lower, _median, _upper, respectively).

# Phenology Models by Ecodata
We've done work spanning on arthropod phenology spanning multiple repos, typically with GAMs using both Bayesian (with the brms library) and frequentist (with the mgcv library) approaches. These repos contain useful functions including: calculating degree days from climate data, scraping climate data, fitting and visualising phenology models, and displaying interactive date/location based phenology curves in Shiny.

## 1. tick_modeling
#### Active May - October 2022
#### By: Michael
#### Repo: https://github.com/ecodata-technology/tick_modeling
#### Highlights:
#### Overview: 

(pic)

## 2. aphid_modeling
#### Active: May - October 2022
#### By: Michael
#### Repo: https://github.com/ecodata-technology/aphid_modeling
#### Highlights:
#### Overview: 
(pic)

## 3. phenogamShiny
#### Active: May 2022, March 2023
#### By: Rob, Michael, Alex
#### Repo: https://github.com/ecodata-technology/phenogamShiny
#### Highlights:
#### Overview: 
Alex built on the code for hosting on AWS and embedding in a marketing one-pager, mainly with added CSS for presentation. (Need to merge on the original still.)

## 4. phenogamML
#### Active: April 2022
#### By: Rob, Michael
#### Repo: https://github.com/ecodata-technology/phenogamML
#### Highlights:
#### Overview: 
`./src/cdd_density_functions`
Contains a daymet scrape using daymetr.
contains some python.

(pic)

## 5. ecodatamisc
#### Active: March 2022 - Jan 2023
#### By: Rob, Tim, Michael
#### Repo: https://github.com/ecodata-technology/ecodatamisc
#### Highlights: 
- plot_gam: plot mgcv (frequentist) gam response with a prediction ribbon. (get example for image).
- add_siteday_cdd: fetch CDD for each site, year, and julian day (using get_daymet rather than daymetr).

#### Overview: 
Library of various functions to streamline and automate common tasks in Ecodata work. Also has a get_daymet helper but superseded by the daymetr library.

(pic)

## 6. iNatML
#### Active: Dec 2021 - Mar 2022
#### By: Rob, Michael
#### Repo: https://github.com/robclark19/iNatML
#### Highlights:
#### Overview: 
(Pic)
