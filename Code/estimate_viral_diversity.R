
library(tidyverse)
library(patchwork)
library(MetBrewer)

setwd("~/Documents/Github/straightlinewasalie/Data")
# Access virus data

library(virionData)
get_versioned_data(version = "17636397", dir_path = "./Temp")
