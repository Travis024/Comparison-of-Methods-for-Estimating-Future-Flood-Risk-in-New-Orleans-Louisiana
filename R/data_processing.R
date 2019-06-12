## process_data.R

process_tg_data <- function(filename_in) {
  
  data.tg <- read.csv(filename_in)
  names(data.tg) <- c('year','month','day','hour','sl')
  
  # any years with less than 90% of the data?
  all_years <- sort(unique(data.tg$year))
  n_data_per_year <- 24*365 # 24 hours/day * 365 days/year (minimum)
  missing_data <- NULL
  for (year in all_years) {
    idx_this_year <- which(data.tg$year==year)
    n_data_this_year <- length(which(data.tg$sl[idx_this_year] > fillvalue))
    if (n_data_this_year/n_data_per_year < threshold_missing_data) {
      missing_data <- c(missing_data, year)
    }
  }
  # for NOLA, missing_data should be NULL
  if (!is.null(missing_data)) {print("WARNING: dropping some years for missing data")}
  
  # which years are good?
  good_years <- setdiff(all_years, missing_data)
  nyear <- length(good_years)
  
  # detrend by subtracting off annual means and get block maxima
  data.tg$sl.detrended <- data.tg$sl
  lsl_max <- rep(NA,nyear)  # annual block maxima
  
  for (year in good_years) {
    idx_this_year <- which(data.tg$year==year & data.tg$sl > fillvalue)
    mean_sl_this_year <- mean(data.tg$sl[idx_this_year])
    data.tg$sl.detrended[idx_this_year] <- data.tg$sl[idx_this_year] - mean_sl_this_year
    lsl_max[which(good_years==year)] <- max(data.tg$sl.detrended[idx_this_year])
  }
  
  output_matrix <- matrix(NA, nrow=nyear, ncol=2)
  colnames(output_matrix) <- c("year", "lsl_max")
  output_matrix[,"year"] <- good_years
  output_matrix[,"lsl_max"] <- lsl_max
  
  return(output_matrix)
}
