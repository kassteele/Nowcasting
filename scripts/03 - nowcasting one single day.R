#
# Nowcasting the number of new symptomatic cases during infectious
# disease outbreaks using constrained P-spline smoothing
#
# Jan van de Kassteele - RIVM
#

#
# Source initial settings ----
#

# Load packages
# Set locale to English
# Source functions
# Set colors and palettes
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
# Model setup ----
#

# Set prior reporting delay distribution
#
# Create a vector with a PMF on the starting date of the outbreak
# The assumption is that only 1 case is expected on the starting date of the outbreak
#
# For the measles outbreak:
# Prior mean reporting delay is 12 days
# Maximum reporting delay is 6 weeks (42 days), where 99% of all cases have been reported
f.priordelay <- genPriorDelayDist(mean.delay = 12, max.delay = 42, p = 0.99)

# Data setup
#
# Create a dataframe with the number of cross-tabulated cases and other things:
# Date      symptoms onset time
# Delay     reporting delay (days)
# t         integer showing the symptoms onset time since the lowest Date
# d         integer showing reporting delay (days)
# Reported  factor with levels: "Reported", "Not yet reported" and, retrospectively, "Future"
# Day       factor with day of the week
# Cases     number of cases
# Est       include record in estimation procedure (1 = yes, 0 = no)
# g         boundary constraint, log(reporting intensity)
# b         boundary constraint indicator (1 = active, 0 = not active)

# Data setup
rep.data <- dataSetup(
  data         = epi.data,
  start.date   = as.Date("2013-05-01"), # Starting date of outbreak
  end.date     = as.Date("2013-09-15"), # Ending date of outbreak (in real-time, leave NULL so end.date = nowcast.date)
  nowcast.date = as.Date("2013-09-01"), # Nowcast date
  days.back    = 2*42,                  # Number of days back from nowcast.date to include in estimation procedure
  f.priordelay = f.priordelay)          # Prior reporting delay PMF

# Model setup
# Default settings
model.setup <- modelSetup(
  data = rep.data,
  ord = 2,
  kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6))

#
# Nowcasting ----
#

# Doen voor 10 juni, 10 juli en 10 augustus
# ord 1 en 2

# Nowcast
nowcast.list <- nowcast(
  data = rep.data,
  model = model.setup,
  conf.level = 0.90)

# Show table
nowcast.list$nowcast

#
# Plots ----
#

# # Prior reporting delay PMF
# ggplot(mapping = aes(x = seq_along(f.priordelay) - 1, y = f.priordelay)) +
#   geom_col(fill = grey) +
#   labs(x = "Reporting delay", y = "Probability mass") +
#   theme_bw()

# Epicurve and reporting trapezium
(epicurve.plot <- plotEpicurve(data = rep.data, title = "Measles outbreak NL 2013-2014"))
(trapezium.plot <- plotTrapezium(data = rep.data, title = "Reporting trapezium"))

# Nowcast and time-varying reporting delay distribution
(nowcast.plot <- plotNowcast(data = rep.data, nowcast = nowcast.list, title = "Nowcast"))
(delaydist.plot <- plotDelayDist(data = rep.data, nowcast = nowcast.list, title = "Reporting delay distribution"))

# Save plots
ggsave(file = "documents/paper/fig_epicurve.png", plot = epicurve.plot, width = 16, height = 9, units = "cm", dpi = 300)
ggsave(file = "documents/paper/fig_trapezium.png", plot = trapezium.plot, width = 16, height = 6.5, units = "cm", dpi = 300)
ggsave(file = "documents/paper/fig_nowcast.png", plot = nowcast.plot, width = 16, height = 9, units = "cm", dpi = 300)
ggsave(file = "documents/paper/fig_delaydist.png", plot = delaydist.plot, width = 16, height = 6.5, units = "cm", dpi = 300)
