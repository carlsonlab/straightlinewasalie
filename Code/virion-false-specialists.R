
library(tidyverse)
library(patchwork)
library(MetBrewer)

library(tidyverse)
library(patchwork)
library(MetBrewer)
library(MoMAColors)
library(vroom)

setwd("~/Github/straightlinewasalie/Data")

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
  rowwise() %>%
  mutate(Year = min(PublicationYear,ReleaseYear,CollectionYear,na.rm=TRUE)) %>% 
  filter(!(Year==Inf)) -> clo

clo %>% select(Host, Virus, Year) %>%
  distinct() %>% 
  filter(!is.na(year)) %>%
  group_by(Host, Virus) %>%
  summarize(Year = min(Year)) -> clo

clo %>% pull(Host) %>% unique() %>% length() -> total_hosts
clo %>% pull(Virus) %>% unique() %>% length() -> total_viruses
clo %>%
  select(Host, Virus) %>% distinct() %>%
  group_by(Virus) %>% summarize(HostRange = n()) %>%
  filter(HostRange==1) %>% pull(Virus) -> true_specialists

##### Charting it through real-time

params <- data.frame(
  year = c(1960:2024),
  hosts = NA, viruses = NA, 
  prop_false = NA
)

for(year in 1960:2024) {
  
  clo %>% filter(Year <= year) -> clo_year
  
  clo_year %>% pull(Host) %>% unique() %>% length() -> hosts
  clo_year %>% pull(Virus) %>% unique() %>% length() -> viruses
  
  clo_year %>% 
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>%
    filter(HostRange==1) %>% pull(Virus) -> obs_specialists
  
  prop_false <- 1-(sum(as.numeric(obs_specialists %in% true_specialists))/length(obs_specialists))
  
  params[params$year == year,] <- c(year, hosts, viruses, prop_false)
  print(year)
}

##### Charting it through sampling

params2 <- data.frame(
  hosts = NA, viruses = NA,
  prop_false = NA
)

hosts <- unique(clo$Host)

for(i in 1:1000) {
  
  n_hosts_i <- sample(c(5:total_hosts),1)
  hosts_i <- sample(hosts, n_hosts_i)
  
  clo %>% filter(Host %in% hosts_i) -> clo_i
  n_viruses_i <- length(unique(clo_i$Virus))
  clo_i %>% 
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>%
    filter(HostRange==1) %>% pull(Virus) -> obs_specialists
  prop_false <- 1-(sum(as.numeric(obs_specialists %in% true_specialists))/length(obs_specialists))
  
  params2[i,] <- c(n_hosts_i, n_viruses_i, prop_false)
  print(i)
}

##### Make plot

params2 %>% 
  ggplot(aes(x = hosts, y = prop_false)) + 
  geom_point() + theme_bw() + 
  geom_point(data = params, aes(x = hosts, y = prop_false, color = year))
