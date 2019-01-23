pitHist <- function(data, start.date, nowcast.dates, days.back, f.priordelay, 
  ord = 2, kappa = list(u = 1e6, b = 1e6, w = 1e-2, s = 1e-6),
  window = 7, J = 10) {
  #
  # Description
  # Function to create a PIT histogram for a sequence of nowcast dates and a number of days back
  #
  # Arguments
  # data           Dateframe with two Date columns: onset.date and report.date
  # nowcast.dates  Sequence of nowcast dates
  # days.back      Number of days back from nowcast.date to include in estimation procedure
  # f.priordelay   Prior delay PMF, from genPriorDelayDist
  # window         Number ofs day to make a nowcast for, counting back from nowcast date
  # J              Number of PIT histogram bins
  #
  # Value
  # Object of class "histogram". See help(hist)
  
  # Require packages
  require(surveillance)
  require(parallel)
  
  # Make list with observed counts N and predictive ecdf's F
  # Do this in parallel
  # NF.list is a list of length(nowcast.dates)
  NF.list <- mclapply(X = nowcast.dates, mc.cores = 4, FUN = function(nowcast.date) {
  #NF.list <- lapply(X = nowcast.dates, FUN = function(nowcast.date) {
    print(nowcast.date)
    # Data setup
    rep.data <- dataSetup(
      data         = data,
      start.date   = start.date,
      end.date     = nowcast.date,
      nowcast.date = nowcast.date,
      days.back    = days.back,
      f.priordelay = f.priordelay)
    # Model setup
    model.setup <- modelSetup(data = rep.data, ord = ord, kappa = kappa)
    # Run the nowcast
    nowcast.list <- nowcast(data = rep.data, model = model.setup, conf.level = 0.90)
    # Add observations and predictive ecdf
    return(list(
      N.pit = tail(aggregate(Cases ~ Date, FUN = sum, data = rep.data), n = window)$Cases,
      F.pit = tail(nowcast.list$F.nowcast, n = window)))
  })
  
  # Glue the list elements together in:
  # - vector with counts N.pit
  # - list of ecdf's F.pit
  N.pit <- do.call(what = c, args = lapply(NF.list, function(x) x$N.pit))
  F.pit <- do.call(what = c, args = lapply(NF.list, function(x) x$F.pit))
  
  # Return PIT histogram
  pit(x = N.pit, pdistr = F.pit, J = J)
}
