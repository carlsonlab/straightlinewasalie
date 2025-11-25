
library(tidyverse)
library(patchwork)
library(MetBrewer)

library(tidyverse)
library(patchwork)
library(MetBrewer)
library(MoMAColors)
library(vroom)

setwd("~/Documents/Github/straightlinewasalie/Data")

# Access virus data
# library(virionData)
# get_versioned_data(version = "17636397", dir_path = "./Temp")

# Set up virus data
virion <- vroom("./Temp/17636397/virion.csv.gz")
virion %>%
  mutate(Host = str_to_lower(Host), Virus = str_to_lower(Virus)) %>%
  filter(DetectionMethod %in% c("PCR/Sequencing","Isolation/Observation")) %>%
  filter(HostClass=='mammalia') %>%
  filter(!(Host=='homo sapiens')) %>%
  filter(!is.na(Host)) %>% filter(!is.na(Virus)) %>%
  select(Host, Virus) %>% 
  distinct() -> clo

Hosts <- length(unique(clo$Host))
Viruses <- length(unique(clo$Virus))

df1 <- data.frame(i = 1:10000, hosts = NA, viruses = NA)
for (i in 1:10000) {
  nhosts <- round(runif(1, min = 1, max = Hosts))
  hosts <- sample(unique(clo$Host), nhosts)
  clo_i <- clo %>% filter(Host %in% hosts)
  nviruses <- length(unique(clo_i$Virus))
  df1[i,] <- c(i, nhosts, nviruses)
  print(i)
}

df2 <- data.frame(i = 1:10000, hosts = NA, viruses = NA, links = NA)
for (i in 1:10000) {
  nrows <- nrow(clo)
  nrow_i <- round(runif(1, min = 1, max = nrows))
  rows <- sample(c(1:nrows), nrow_i)
  clo_i <- clo[rows,]
  nhosts <-length(unique(clo_i$Host))
  nviruses <- length(unique(clo_i$Virus))
  df2[i,] <- c(i, nhosts, nviruses, nrow_i)
  print(i)
}

thomas <- met.brewer("Thomas", n=8)
df1 %>% 
  ggplot(aes(x = hosts, y = viruses)) + 
  theme_bw() + 
  geom_point(alpha = 0.02, col = thomas[8]) + 
  geom_abline(intercept = 0, slope = slope, col = 'lightgrey') +
  ggtitle('Complete sampling') + xlab("Hosts") + ylab("Viruses") -> g1
df2 %>% 
  ggplot(aes(x = hosts, y = viruses)) + 
  geom_abline(intercept = 0, slope = slope, col = 'lightgrey') + 
  theme_bw() + geom_point(alpha = 0.02, col = thomas[2]) + 
  ggtitle('Partial sampling') + xlab("Hosts") + ylab("Viruses") -> g2

g2+g1
