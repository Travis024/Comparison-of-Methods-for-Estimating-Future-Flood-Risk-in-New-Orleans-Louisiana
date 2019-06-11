# Comparison of Methods for Estimating Future Flood Risk in New Orleans Louisiana
A summer research project conducted by Anthony Wong and Travis Torline at the University of Colorado Boulder

This project will use sea level data from New Orleans, Louisiana to make improved projections of flood risk throughout the 21st century. The currently-used calculation of flood hazard relies on some “magic numbers”. We will estimate appropriate values for those parameters using an independent and statistically rigorous approach.

As sea levels and storm surges rise due to climate change, this project will highlight important modeling uncertainties and provide recommendations for augmenting the flood protection system in New Orleans. These results will benefit the flood risk modeling community as well as the local population of New Orleans.

## Week 1

Work through starter problems to get a better feel for the R programming language (including statistical functions like quantile, graphig function like polygon, and list objects), the extRemes packages, and NetCDF files (a common file type used in meteorology, climatology, and oceanography). This will begin to give me an idea of the types of data we'll be using in our analysis and hot to clean and work with that data.

View answers to the Starter Problems in the "StarterProblems.Rmd" file.

Learn more about the exTremes package: https://cran.r-project.org/web/packages/extRemes/extRemes.pdf
Learn more about NetCDF files: https://www.unidata.ucar.edu/software/netcdf/

## Week 2

Begin analysing prior reserach conducted by Mingxuan Zang concerning 8 different GEV models and how well they fit the data for the past 40 years of sea levels in New Orleans, Louisiana. These models contain mu (location), sigma (scale), and xi (shape) parameters as either stationary (denoted with a '0') or non-stationary (denoted with a '1'). Mingxuan worked on this research with Tony previously.

This data is uplaoded in the folder parameters_data.
