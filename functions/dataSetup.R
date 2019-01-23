dataSetup <- function(data, start.date, end.date = NULL, nowcast.date, days.back = NULL, f.priordelay) {
  #
  # Data setup
  #
  # Description
  # Function that sets up the data for nowcasting
  #
  # Arguments
  # data          Dateframe with two Date columns: onset.date and report.date
  # start.date    Starting date of outbreak
  # end.date      Ending date of outbreak
  #               In real-time, leave NULL so end.date = nowcast.date
  # nowcast.date  Nowcast date
  # days.back     Number of days back from nowcast.date to include in estimation procedure
  #               If NULL, it is set to two times the number of days in f.priordelay
  # f.priordelay  Prior delay PMF, from genPriorDelayDist
  #
  # Value
  # Dataframe with:
  # Date      date of disease onset
  # Delay     delay (days)
  # Reported  factor with levels: "Reported", "Not yet reported" and, retrospectively, "Future"
  # Day       factor with day of the week
  # Cases     number of cases
  # Est       include record in estimation procedure (1 = yes, 0 = no)
  # t         number of days since start.date
  # d         delay (days)
  # g         boundary constraint, log(reporting intensity)
  # b         boundary constraint indicator (1 = active, 0 = not active)
  
  #
  # Initial stuff
  #
  
  # If there is no end.date, set end.date equal to nowcast.date
  if (is.null(end.date)) end.date <- nowcast.date
  
  # Get maximum delay
  max.delay <- length(f.priordelay) - 1
  
  # If there is no days.back, set days.back to two times max.delay
  if (is.null(days.back)) days.back <- 2*max.delay
  
  # Get the dimensions of the reporting trapezium (= T x D1) and the entire outbreak (= T.true x D1)
  T      <- as.numeric(nowcast.date - start.date) + 1 # Number of days from start.date to nowcast.date
  T.true <- as.numeric(    end.date - start.date) + 1 # Number of days from start.date to end.date (truth, retrospectively)
  D      <- max.delay                                 # Maximum delay in days
  D1     <- D + 1                                     # Number of days from 0 to max.delay
  
  # Set vectors t, t.true and d
  t      <- 1:T      # Days since start.date: 1, 2, ..., T
  t.true <- 1:T.true # Days since start.date: 1, 2, ..., T.true
  d      <- 0:D      # Delays 0, 1, ..., max.delay
  
  #
  # Data operations
  #
  
  data <- data %>%
    # Filter records with start.date <= onset.data <= end.date
    filter(onset.date >= start.date & onset.date <= end.date) %>% 
    
    # Compute delay in days
    mutate(delay = (report.date - onset.date) %>% as.numeric) %>% 
    
    # Filter records with 0 <= delay <= max.delay
    filter(delay >= 0 & delay <= max.delay) %>% 
    
    # Categorize onset.date and delay
    # We need this to tabulate the cases by onset.date and delay
    mutate(
      onset.date.cat = onset.date  %>% cut(breaks = seq(from = start.date, to = end.date + 1, by = "day")),
      delay.cat      = delay %>% factor(levels = 0:max.delay)) %>% 
    
    # Remove (numeric) delay
    select(-delay)
  
  #
  # Construct reporting data
  #
  
  # Setup the reporting trapezium data as a grid by Date and Delay
  rep.data <- expand.grid(
    Date = data$onset.date.cat %>% levels %>% as.Date,
    Delay = d) %>% 
    
    mutate(
      # Add t (t.true, actually) and d, to assist in the calculations
      t = (Date - start.date + 1) %>% as.integer,
      d = Delay,
      
      # Add reporting category: Reported, Not yet reported, Future
      # Cases with t + d <= T have been reported
      Reported = ifelse(t + d <= T, yes = "Reported",
        # Cases with t > T are in the future
        no = ifelse(t > T, yes = "Future",
          # The rest has not been reported yet
          no = "Not yet reported")) %>% 
        # As factor
        factor(levels = c("Reported", "Not yet reported", "Future")),
      
      # Add day of the week
      # [t = 3, d = 0], [t = 2, d = 1], [t = 1, d = 2] = constant, etc., so we have
      # Monday is reference (trick: 2007-01-01 is Monday)
      Day = weekdays(
        x = Date + Delay,
        abbreviate = TRUE) %>%
        factor(
          levels = weekdays(
            x = seq(
              from = as.Date("2007-01-01"),
              to = as.Date("2007-01-07"),
              by = "1 day"),
            abbreviate = TRUE)),
      
      # Add tabulated cases by date and delay
      Cases = with(data, table(onset.date.cat, delay.cat)) %>% as.vector,
      
      # Include record in estimation procedure? (1 = yes, 0 = no)
      # This is where nowcast.date - days.back + 1 <= Onset.Date <= nowcast.date
      Est = (Date >= (nowcast.date - days.back + 1) & Date <= nowcast.date) %>% as.integer)

  #
  # Add boundary constraints to reporting data
  #

  # g is value to keep surface below eta <= g
  # b is where the constraint is active
  g <- matrix(0, nrow = T.true, ncol = D1)
  b <- matrix(0, nrow = T.true, ncol = D1)
  # Set g and b at t = 1
  g[1, ] <- log(f.priordelay)
  b[1, ] <- 1
  # Set g and b at max.delay
  g[, D1] <- log(f.priordelay[D1])
  b[, D1] <- 1
  # Add these to rep.data
  rep.data <- rep.data %>%
    mutate(
      b = b %>% as.vector,
      g = g %>% as.vector)

  #
  # Return output
  #
  
  return(rep.data)
}
