rdate <- function(x,
                  min_date = '2020-01-01',
                  max_date = '2021-12-31',
                  sort = TRUE) {

  dates <- sample(seq(as.Date(min_date), as.Date(max_date), by = "day"), x, replace = TRUE)
  if (sort == TRUE) {
    sort(dates)
  } else {
    dates
  }

}

rdate_increase <- function(dsdate, start=0, end=90) {
  ### Generates a fictional end date for an event of interest
  random_days <- sample(0:90, 1)
  gdate <- as.Date(dsdate) + days(random_days)
  gdate
}

rtime <- function(n=200000, start=7, end=19, lld=0.9) {
  # Specify the time range in minutes
  start_time <- 7 * 60  # 7 AM in minutes
  end_time <- 19 * 60   # 7 PM in minutes
  # Generate random minutes and seconds with higher likelihood for sampling between 7AM and 7PM
  random_minutes <- c(sample(0:(start_time - 1), round(n * 0.1), replace = TRUE),
                      sample(start_time:end_time, round(n * 0.9), replace = TRUE))
  random_seconds <- sample(0:59, n, replace = TRUE)
  times_list <- lapply(1:n, function(i) {
    paste(sprintf("%02d", random_minutes[i] %/% 60),
          sprintf("%02d", random_minutes[i] %% 60), sep = ":")
  })
  unlist(times_list)
}
