#
# Init ----
#

# Load packages
library(tidyverse)
library(openxlsx)

#
# Read epidata ----
#

# We should end up with a dataframe with two columns:
# onset.date  = symptoms onset time
# report.date = date of report

# Read epi data from file
epi.data <- read.xlsx(xlsxFile = "data/Mazelen_rapport_140930.xlsx", sheet = 2)

# Some filtering and modifications
epi.data <- epi.data %>%
  # EPILANDBESMETTING: Nederland
  filter(EPILANDBESMETTING == "Nederland") %>%
  mutate(
    # ZIE1EZIEKTEDT     = disease onset date before 2011-06-24
    # MORBZIE1EZIEKTEDT = disease onset date after 2011-06-24
    # CREATIEDT         = report date
    # Make dates a Date class (was POSIXct class)
    onset.date1 = ZIE1EZIEKTEDT     %>% strptime(format = "%d-%m-%Y") %>% as.Date,
    onset.date2 = MORBZIE1EZIEKTEDT %>% strptime(format = "%d-%m-%Y") %>% as.Date,
    onset.date  = ifelse(test = !is.na(onset.date1), yes = onset.date1, no = onset.date2) %>% as.Date(origin = "1970-01-01"),
    report.date = CREATIEDT         %>% strptime(format = "%d-%m-%Y") %>% as.Date) %>%
  # Drop all other variables
  select(onset.date, report.date) %>%
  # NL 2013-2014 measles outbreak
  # About the outbreak: http://www.rivm.nl/Onderwerpen/M/Mazelen/Mazelenepidemie_2013_2014
  # Start date: 2013-05-01
  # End date:   2014-02-28
  filter(onset.date >= as.Date("2013-05-01") & onset.date <= as.Date("2014-02-28")) %>%
  # Optional: re-order rows by onset.date, and within onset.date by report.date
  arrange(onset.date, report.date)

#
# Save epidata ----
#

# Save dataset for further use
epi.data %>% 
  write.table(file = "data/measles_NL_2013_2014.dat", sep = "\t", quote = FALSE, row.names = FALSE)
