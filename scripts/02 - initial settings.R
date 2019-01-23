#
# Init ----
#

# Load packages
library(Matrix)
library(scales)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

# Set time locale to English for English day and month names
if (Sys.info()["sysname"] == "Windows") {
  # Windows
  Sys.setlocale(category = "LC_TIME", locale = "English")
} else {
  # Linux
  Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8")
}

# Source functions
list.files(path = "functions", full.names = TRUE) %>% 
  walk(source)

# Set colors
blue <- brewer.pal(n = 9, name = "Set1")[2]
oran <- brewer.pal(n = 9, name = "Set1")[5]
grey <- brewer.pal(n = 9, name = "Set1")[9]

# Set color palettes
blue.pal <- brewer.pal(n = 9, name = "Blues")   %>% colorRampPalette
oran.pal <- brewer.pal(n = 9, name = "Oranges") %>% colorRampPalette
grey.pal <- brewer.pal(n = 9, name = "Greys")   %>% colorRampPalette
