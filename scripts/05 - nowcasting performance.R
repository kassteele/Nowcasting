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

#
# Create PIT histograms ----
#

# Set Prior reporting delay distribution
f.priordelay <- genPriorDelayDist(mean.delay = 12, max.delay = 42, p = 0.99)

# Create PIT histograms for three phases
# _d1 = 1st order differences in time direction
# _d2 = 2nd order differences
phase1_d1 <- pitHist(
  data = epi.data,
  start.date = as.Date("2013-05-01"),
  nowcast.dates = seq(
    from = as.Date("2013-06-03"),
    to =   as.Date("2013-06-17"),
    by = "1 day"),
  days.back = 2*42,
  f.priordelay = f.priordelay,
  ord = 1, kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6),
  window = 7)

phase2_d1 <- pitHist(
  data = epi.data,
  start.date = as.Date("2013-05-01"),
  nowcast.dates = seq(
    from = as.Date("2013-07-03"),
    to =   as.Date("2013-07-17"),
    by = "1 day"),
  days.back = 2*42,
  f.priordelay = f.priordelay,
  ord = 1, kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6),
  window = 7)

phase3_d1 <- pitHist(
  data = epi.data,
  start.date = as.Date("2013-05-01"),
  nowcast.dates = seq(
    from = as.Date("2013-08-03"),
    to =   as.Date("2013-08-17"),
    by = "1 day"),
  days.back = 2*42,
  f.priordelay = f.priordelay,
  ord = 1, kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6),
  window = 7)

phase1_d2 <- pitHist(
  data = epi.data,
  start.date = as.Date("2013-05-01"),
  nowcast.dates = seq(
    from = as.Date("2013-06-03"),
    to =   as.Date("2013-06-17"),
    by = "1 day"),
  days.back = 2*42,
  f.priordelay = f.priordelay,
  ord = 2, kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6),
  window = 7)

phase2_d2 <- pitHist(
  data = epi.data,
  start.date = as.Date("2013-05-01"),
  nowcast.dates = seq(
    from = as.Date("2013-07-03"),
    to =   as.Date("2013-07-17"),
    by = "1 day"),
  days.back = 2*42,
  f.priordelay = f.priordelay,
  ord = 2, kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6),
  window = 7)

phase3_d2 <- pitHist(
  data = epi.data,
  start.date = as.Date("2013-05-01"),
  nowcast.dates = seq(
    from = as.Date("2013-08-03"),
    to =   as.Date("2013-08-17"),
    by = "1 day"),
  days.back = 2*42,
  f.priordelay = f.priordelay,
  ord = 2, kappa = list(u = 1e6, b = 1e6, w = 0.01, s = 1e-6),
  window = 7)

#
# PIT histograms ----
#

# Create pit.data for plotting
labs <- expand.grid(
  phase = c(
    "Growth phase, June 2013",
    "Peak phase, July 2013",
    "Decline phase, August 2013"),
  ord = c(
    "1st order differences",
    "2nd order differences"))
pit.data <- cbind(
  data.frame(
    mids = c(
      phase1_d1$mids,   phase2_d1$mids,   phase3_d1$mids,
      phase1_d2$mids,   phase2_d2$mids,   phase3_d2$mids),
    counts = c(
      phase1_d1$counts, phase2_d1$counts, phase3_d1$counts,
      phase1_d2$counts, phase2_d2$counts, phase3_d2$counts)),
  labs[rep(1:6, each = 10), ])

# Show mean absolute difference with 1
pit.data %>%
  mutate(mad = abs(counts - 1)) %>% 
  group_by(ord, phase) %>% 
  summarize(mad = mean(mad) %>% round(digits = 2))

# 1 1st order differences Growth phase, June 2013     0.84
# 2 1st order differences Peak phase, July 2013       0.87
# 3 1st order differences Decline phase, August 2013  0.72
# 4 2nd order differences Growth phase, June 2013     0.7 
# 5 2nd order differences Peak phase, July 2013       0.23
# 6 2nd order differences Decline phase, August 2013  0.62

#
# Plot ----
#

# Plot PIT histograms
ggplot(data = pit.data, mapping = aes(x = mids, y = counts)) +
  geom_hline(yintercept = 1, lty = 2, size = 0.2) +
  geom_col(fill = grey) +
  facet_grid(rows = vars(ord), cols = vars(phase)) +
  labs(x = "PIT", y = "Relative frequency") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())

# Save plot
ggsave(
  filename = "documents/paper/fig_pithist.png",
  width = 16, height = 11, units = "cm", dpi = 300)


