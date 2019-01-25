#
# Nowcasting the number of new symptomatic cases during infectious
# disease outbreaks using constrained P-spline smoothing
#
# Jan van de Kassteele - RIVM
#

#
# Source initial settings ----
#

source(file = "scripts/02 - initial settings.R")

#
# Read epidata ----
#

# Read pre-processed epidata
epi.data <- read.delim(file = "data/measles_NL_2013_2014.dat") %>% 
  # Make sure onset.date and report.date a Date class (was factor)
  mutate(
    onset.date  = onset.date  %>% as.Date,
    report.date = report.date %>% as.Date)

#
# Nowcasting loop ----
#

# Set prior delay distribution
f.priordelay <- genPriorDelayDist(mean.delay = 12, max.delay = 42, p = 0.99)

# Initial nowcast date
nowcast.date <- as.Date("2013-05-27")

# Open multipage pdf
pdf(file = "nowcasts_d2.pdf", width = 29.7/2.54, height = 21/2.54, paper = "a4r", onefile = TRUE)

# Loop over nowcast dates until specified date
while (nowcast.date <= as.Date("2013-10-01")) {
  
  # Print progress
  print(nowcast.date)
  
  # Data setup
  rep.data <- dataSetup(
    data         = epi.data,
    start.date   = as.Date("2013-05-01"), # Starting date of outbreak
    end.date     = as.Date("2013-11-01"), # Ending date of outbreak (in real-time, leave NULL so end.date = nowcast.date)
    nowcast.date = nowcast.date,          # Nowcast date
    days.back    = 2*42,                  # Number of days back from nowcast.date to include in estimation procedure
    f.priordelay = f.priordelay)          # Prior delay PMF
  
  # Model setup
  model.setup <- modelSetup(
    data = rep.data,
    ord = 2,
    kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6))
  
  # Nowcast
  nowcast.list <- nowcast(data = rep.data, model = model.setup, conf.level = 0.90)
  
  # Make nowcast plot
  nowcast.plot <- plotNowcast(data = rep.data, nowcast = nowcast.list, title = "Nowcast")
  plot(nowcast.plot)

  # Add 1 day
  nowcast.date <- nowcast.date + 1
  
  # End loop over nowcast dates
}

# Close pdf
dev.off()
