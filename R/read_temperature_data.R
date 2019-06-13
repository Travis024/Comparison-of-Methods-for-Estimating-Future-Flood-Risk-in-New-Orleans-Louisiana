##=============================================================================
## read_temperature_data.R
##
## read temperature data
## historical from NOAA:
##   NOAA National Centers for Environmental information, Climate at a Glance:
##   Global Time Series, published May 2017, retrieved on June 7, 2017 from
##   http://www.ncdc.noaa.gov/cag/
## future projections from CRNM under RCP8.5, part of CMIP5
##
## Questions? Tony Wong (anthony.e.wong@colorado.edu)
##=============================================================================

# read the projections temperature forcing from CRNM (CMIP5)
# note: these are in K, but going to normalize, so will take a difference and
# it's same as celsius
ncdata <- nc_open('../parameters_data/global.tas.aann.CNRM-CM5.historical+rcp85.r1i1p1.18500101-21001231.nc')
temperature_proj <- ncvar_get(ncdata, 'tas')
time_proj <- ncvar_get(ncdata, 'time')
nc_close(ncdata)

# 'time' on the netcdf ile is YYYYMMDD, where MMDD is July 2 each year (becuase
# of the averaging). so pluck off just the years
time_proj <- floor(time_proj/10000)

# read the historical forcing from NOAA
data.tmp <- read.table('../parameters_data/noaa_temperature_1880-2017.csv', header = TRUE, sep=',')
time_hist <- data.tmp$Year
temperature_hist <- data.tmp$Value

# extend historical back to 1850 with HadCRUT4
data.tmp = read.table('../parameters_data/HadCRUT.4.4.0.0.annual_ns_avg.txt')
time_hadcrut = data.tmp[,1]
temperature_hadcrut = data.tmp[,2]

# normalize all to 1901-2000
ind_norm <- which(time_hist==1901):which(time_hist==2000)
temperature_hist <- temperature_hist - mean(temperature_hist[ind_norm])

ind_norm <- which(time_hadcrut==1901):which(time_hadcrut==2000)
temperature_hadcrut <- temperature_hadcrut - mean(temperature_hadcrut[ind_norm])

ind_norm <- which(time_proj==1901):which(time_proj==2000)
temperature_proj <- temperature_proj - mean(temperature_proj[ind_norm])

# set up the forcing, wherein the historical starts, then projections finish
# 1850-1880 -- hadcrut4
# 1880-2016 -- NOAA historical
# 2017-2100 -- CRNM projection
time_forc <- min(time_hadcrut):max(time_proj)
temperature_forc <- rep(NA, length(time_forc))

ind_hadcrut <- which(time_hadcrut==time_forc[1]):(which(time_hadcrut==time_hist[1])-1)
ind_hist    <- 1:length(time_hist)
ind_proj    <- which(time_proj==(max(time_hist)+1)):which(time_proj==max(time_forc))

temperature_forc <- c(temperature_hadcrut[ind_hadcrut],
                      temperature_hist[ind_hist]      ,
                      temperature_proj[ind_proj]      )

# maximum temperature serves as an additinoal prior constraint on kappa0, kappa1
# that is, kappa1 > -kappa0/Tmax (otherwise, kappa = kappa0 + kappa1*T might be
# negative)  (not necessary for GEV models)
Tmax <- max(temperature_forc)

# Useful:

# function to trim temperature forcing to fit TG record unique years
# there might be missing years in TG record, so need to match each year and not
# just plop down an evenly spaced sequence
trimmed_forcing <- function(years_tidegauge, years_temperature, temperature) {
  output <- vector('list', 2); names(output) <- c('time','temperature')
  # check the beginning
  if(years_temperature[1] > years_tidegauge[1]) {print('ERROR - tide gauge record starts before temperature; add support for this situation')}
  # check the end
  if(max(years_temperature) < max(years_tidegauge)) {print('ERROR - tide gauge record ends after temperature; add support for this situation')}
  # match the indices of years_tidegauge within years_temperature
  imatch <- match(years_tidegauge, years_temperature)
  output$time <- years_temperature[imatch]
  output$temperature <- temperature[imatch]
  return(output)
}

#===============================================================================
# End
#===============================================================================
