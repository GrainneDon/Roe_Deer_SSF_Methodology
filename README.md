Thermoregulatory Habitat Selection in Roe Deer (Capreolus capreolus) at Grimsö, Sweden
Author: Gráinne Margaret Donohue
Programme: MSc Geomatics with Remote Sensing and GIS
Institution: Stockholm University
Year: Spring 2026

Project Overview
This repository contains the R analysis script for a master's thesis investigating whether roe deer exhibit thermoregulatory habitat selection during an acute summer heat event at Grimsö, Örebro County, Sweden. The study uses a step selection function (SSF) framework adapted from Alston et al. (2020), comparing roe deer habitat selection and movement behaviour between a three-day heat event and a cool comparison period in July 2024.

Repository Contents

SSF_Clean.R — the full analysis script, including raw data processing, model fitting, RSS visualisation, step length analysis, habitat suitability mapping, KDE home ranges, and zonal statistics


R Version and Required Packages
This script was written and tested in R version 4.5.2.
The following packages are required:
rlibrary(amt)
library(terra)
library(tidyverse)
library(survival)
library(lubridate)
library(sf)
library(data.table)
library(spatialEco)
library(Hmisc)
library(corrplot)
library(car)
library(patchwork)
library(adehabitatHR)

Data Requirements
The following input files are required to run the script. These are not included in this repository.
Here is the updated data requirements section without file paths:

Data Requirements
The following input files are required to run the script. These are not included in this repository due to data sharing restrictions.

GPS collar data:

Cleaned GPS fixes for the heat event period
Cleaned GPS fixes for the cool comparison period

Temperature data:

Hourly air temperature from Kloten A weather station, downloaded from SMHI
Microclimate logger network data with coordinates, provided by Stockholm University

Spatial covariate rasters (SWEREF99 TM, 10m resolution):

Copernicus High Resolution Layer Tree Cover Density (2023)
Heat Load Index derived from Lantmäteriet Nationell Höjdmodell DEM
Euclidean distance to water features derived from Geofabrik OSM hydrographic data
SLU soil moisture dataset


Data Availability
The GPS collar data used in this study was provided by Professor Petter Kjellander, SLU Department of Ecology at Grimsö Research Station, and is not publicly available. The microclimate logger data was provided by Caroline Greiser, Stockholm University. Requests for access to the raw data should be directed to the data owners.
Weather station data is freely available from SMHI at: https://www.smhi.se
Copernicus Tree Cover Density data is freely available at: https://land.copernicus.eu
Lantmäteriet DEM data is available at: https://www.lantmateriet.se

Study Periods

Heat event: 20th–22nd July 2024 (daily maxima 24.9°C, 23.2°C, 19.6°C)
Cool comparison period: 4th–6th July 2024 (daily maxima 14.8°C, 14.5°C, 15.2°C)


Reference
Alston, J. M., Joyce, M. J., Merkle, J. A., & Moen, R. A. (2020). Temperature shapes movement and habitat selection by a heat-sensitive ungulate. Landscape Ecology, 35(9), 1961–1973. https://doi.org/10.1007/s10980-020-01072-y
