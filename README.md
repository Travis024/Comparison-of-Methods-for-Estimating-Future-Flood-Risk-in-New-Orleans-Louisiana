# Comparison of Methods for Estimating Future Flood Risk in New Orleans Louisiana
A summer research project conducted by Anthony Wong and Travis Torline at the University of Colorado Boulder

This project will use sea level data from New Orleans, Louisiana to make improved projections of flood risk throughout the 21st century. The currently-used calculation of flood hazard relies on some “magic numbers”. We will estimate appropriate values for those parameters using an independent and statistically rigorous approach.

As sea levels and storm surges rise due to climate change, this project will highlight important modeling uncertainties and provide recommendations for augmenting the flood protection system in New Orleans. These results will benefit the flood risk modeling community as well as the local population of New Orleans.

## Week 1

Work through starter problems to get a better feel for the R programming language (including statistical functions like quantile, graphig function like polygon, and list objects), the extRemes packages, and NetCDF files (a common file type used in meteorology, climatology, and oceanography). This will begin to give me an idea of the types of data we'll be using in our analysis and hot to clean and work with that data.

View answers to the Starter Problems in the "StarterProblems" folder.

Learn more about the exTremes package: https://cran.r-project.org/web/packages/extRemes/extRemes.pdf

Learn more about NetCDF files: https://www.unidata.ucar.edu/software/netcdf/

## Week 2 - 3

A GEV distribution has three parameters - mu (location), sigma (scale), and xi (shape). The values of these parameters can either be stationary (meaning that they don't change at all) or non-stationary (meaning that the value changes as a function of some other variable).

For this research, non-stationary values change as a function of global average temperature.

In this portion of the research, we analyze 8 .csv files of stationary and non-stationary parameter sets created by one of Tony's former students, Mingxuan Zang. Mingxuan whas uplaoded these files in the "parameters_data" folder along with the best fitting parameter sets. Each .csv file has a different combination of stationary and non-stationary parameters and 2064 values for those parameters.

It is our task to create a density distribution for each of the 8 .csv files showing 100-year storm levels for all of the 2064 parameter sets. A 100-year storm level is a storm with a sea level height that is only seen once every 100 years, and it is the minimum that cities should be preparing for.

Furthermore, we wish to show how these density distribtuions change over time as global mean temperature increases. By comapring the year 2016 to the year 2065, we can see hwo the ditribution of 100-year storms changes and whether we can expect to see higher or lower flood levels for these storms.

This portion of the project can be found in the 'x100Calculations' folder. Specifically, look for the file titled "x100Calculations".

## Week 4

The United States Army Corps of Engineers (USACE for short) has their own model for determining what storm levels and storm surges will look like in the future. This model involves calcualting a change in sea level, pulling stationary GEV parameters, and sampling a surge factor. What isn't clear about this method, however, is where the calulcation for the Surge Factor comes from. It appears to be pulled from a uniform distribution between 1.5 and 2, but this is not clearly stated anywhere.

This week, we will begin to model the USACE against our 8 parameter sets to see how they differ. This will be done by sampling changes in sea level from 2016 to 2065, sampling stationary GEV parameters, and samplng a surge factor from a unfiform distributoiion to calcuate 100-year storm level values using the USACE's method. We will then compare the density distribution of values from this method to our own method.

Our hope is to get the distribtuons to somewhat match up by determing what the distribtuon of surge factors should be. Should it be a uniform ditribution from 1.5 to 2, or maybe from 0 to 10? Is it possible that it's not even a uniform distribtuion at all?

This portion of the project can be found in the 'x100Calculations' folder. Specifically, look for the file titled "USACEFunction" for the function that generates values for the USACE's method and in the file titled "Driver File" for a slick and fast way to plot non-stationary distirbutions and the USACE method against one another.
