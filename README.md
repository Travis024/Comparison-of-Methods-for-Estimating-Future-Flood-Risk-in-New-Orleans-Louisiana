# Comparison of Methods for Estimating Future Flood Risk in New Orleans Louisiana
A summer research project conducted by Anthony Wong and Travis Torline at the University of Colorado Boulder

This project will use sea level data from New Orleans, Louisiana to make improved projections of flood risk throughout the 21st century. The currently-used calculation of flood hazard relies on some “magic numbers”. We will estimate appropriate values for those parameters using an independent and statistically rigorous approach.

As sea levels and storm surges rise due to climate change, this project will highlight important modeling uncertainties and provide recommendations for augmenting the flood protection system in New Orleans. These results will benefit the flood risk modeling community as well as the local population of New Orleans.

## Week 1

Work through starter problems to get a better feel for the R programming language (including statistical functions like quantile, graphig function like polygon, and list objects), the extRemes packages, and NetCDF files (a common file type used in meteorology, climatology, and oceanography). This will begin to give me an idea of the types of data we'll be using in our analysis and hot to clean and work with that data.

View answers to the Starter Problems in the "StarterProblems" folder.

Learn more about the exTremes package: https://cran.r-project.org/web/packages/extRemes/extRemes.pdf
Learn more about NetCDF files: https://www.unidata.ucar.edu/software/netcdf/

## Week 2

Begin analysing prior reserach conducted by Mingxuan Zang concerning 8 different GEV models and how flood risk will increase alongside global mean annual temperature for New Orleans, Louisiana. These models contain mu (location), sigma (scale), and xi (shape) parameters as either stationary (denoted with a '0') or non-stationary (denoted with a '1'). Mingxuan worked on this research with Tony previously.

The goal of my work is first to read in 8 files (one for each model) each containing 2064 sets of parameters. For each set of paramteters, I will find the '1 in 100' flood risk point, or the point where we risk seeing a flood of a height X every 1 in 100 years. I will do this for the years 2017 and 2065 to see how the '1 in 100' flood risk point changes over time.

This portion of the project can be found in the 'x100Calculations' folder.
