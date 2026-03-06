
library(tidyverse)
library(vroom)
library(iNEXT)

setwd("~/Github/straightlinewasalie/Data")

set.seed(42)

# Set up virus data
virion <- vroom("./Temp/17636397/virion.csv.gz")
virion %>%
  mutate(Host = str_to_lower(Host), Virus = str_to_lower(Virus)) %>%
  filter(DetectionMethod %in% c("PCR/Sequencing","Isolation/Observation")) %>%
  filter(HostClass=='mammalia') %>%
  # filter(!(Host=='homo sapiens')) %>%
  filter(!is.na(Virus), !is.na(Host)) %>%
  select(Host, Virus) %>% 
  distinct() -> clo

clo %>% 
  count(Virus) %>%
  count(n) -> counts
h = length(unique(clo$Host))
chao2 = sum(counts$nn) + ((h-1)/h)*(counts$nn[1]^2)/(2*counts$nn[2])

incidence <- clo %>% 
  mutate(Count = 1) %>% 
  pivot_wider(names_from = Host, values_from = Count) %>% 
  replace(is.na(.), 0) %>%
  select(-Virus) %>% as.matrix()
estimateD(incidence, q = 0, datatype = "incidence_raw", level = 6759, base = 'size')

clo %>%
  count(Host) %>% mutate(singleton = (n==1)) %>%
  summarize(singleton)
