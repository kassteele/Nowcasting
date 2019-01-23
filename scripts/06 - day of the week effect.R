#
# Nowcasting the number of new symptomatic cases during infectious
# disease outbreaks using constrained P-spline smoothing
#
# Jan van de Kassteele - RIVM
#

#
# Source initial settings ----
#

# Libraries, settings, etc.
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

# Set prior delay distribution
f.priordelay <- genPriorDelayDist(mean.delay = 12, max.delay = 49, p = 0.99)

# Data setup
rep.data <- dataSetup(
  data         = epi.data,
  start.date   = as.Date("2013-05-01"), # Starting date of outbreak
  end.date     = as.Date("2013-11-01"), # Ending date of outbreak (in real-time, leave NULL so end.date = nowcast.date)
  nowcast.date = as.Date("2013-08-10"), # Nowcast date
  days.back    = 2*49,                  # Number of days back from nowcast.date to include in estimation procedure
  f.priordelay = f.priordelay)          # Prior delay PMF

# Model setup
model.setup <- modelSetup(
  data = rep.data,
  ord = 2,
  kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6))

# Now manually:
data <- rep.data
model <- model.setup
conf.level <- 0.90
# Open nowcast.R and run from 23 - 93

#
# Day-of-the-week effects
#

# Set days
days <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")

# Rate ratios
RR <- rbind(
  rep(1, 3), 
  alpha.beta.sim[(Ks + 1):(Ks + Kw), ] %>% apply(MARGIN = 1, quantile, prob = c(0.5, 0.025, 0.975)) %>% exp %>% t)
RR <- cbind(
  data.frame(factor(days, levels = days)),
  RR)
names(RR) <- c("Day", "RR", "Lower", "Upper")

# Plot
ggplot(data = RR, mapping = aes(x = Day, y = RR, ymin = Lower, ymax = Upper)) +
  geom_col(fill = grey) +
  geom_errorbar(width = 0.2) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.1, by = 0.2)) +
  labs(y = "Rate ratio") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank())

# Save
ggsave(file = "documents/paper/fig_dayeffect.png", width = 16, height = 8, units = "cm", dpi = 300)
