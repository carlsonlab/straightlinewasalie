
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

params <- data.frame(
  year = c(1960:2024),
  viruses = NA, hosts = NA, 
  sigma = NA, eta = NA
)

for(year in 1960:2024) {
  
  clo %>% filter(Year <= year) -> clo_year
  
  clo_year %>%
    pull(Virus) %>% unique() %>% length() -> viruses
  clo_year %>%
    pull(Host) %>% unique() %>% length() -> hosts
  
  clo_year %>%
    group_by(Virus) %>%
    count() %>% pull(n) %>% mean() -> eta 
  
  clo_year %>%
    group_by(Host) %>%
    count() %>% pull(n) %>% mean() -> sigma  
  
  params[params$year == year,] <- c(year, viruses, hosts, sigma, eta)
  
}

# ['#feebe2','#fbb4b9','#f768a1','#c51b8a','#7a0177']


##########################################
##########################################
##########################################

nHosts <- length(unique(clo$Host))
nViruses <- length(unique(clo$Virus))

params2 <- data.frame(
  viruses = NA, hosts = NA, 
  sigma = NA, eta = NA
)

for (i in 5:nHosts){
  
  hosts.i <- sample(unique(clo$Host), i)
  clo.i <- clo %>% filter(Host %in% hosts.i)
  
  clo.i %>%
    pull(Virus) %>% unique() %>% length() -> viruses
  clo.i %>%
    pull(Host) %>% unique() %>% length() -> hosts
  clo.i %>%
    group_by(Virus) %>%
    count() %>% pull(n) %>% mean() -> eta 
  clo.i %>%
    group_by(Host) %>%
    count() %>% pull(n) %>% mean() -> sigma  
  
  params2[i-4,] <- c(viruses, hosts, sigma, eta)
  print(i)
  
}

beths <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")


params %>%
  ggplot(aes(x = hosts, y = viruses, color = year)) + 
  geom_point(size = 3) + 
  theme_bw() +
  xlab('\n Hosts') + ylab('Viruses \n') +
  scale_color_gradient(low=beths[3], high=beths[7], name="Year") -> g1; g1

params %>%
  ggplot(aes(x = sigma, y = eta, color = year)) + 
  geom_point(size = 3) + 
  theme_bw() +
  xlab('\n Average virus diversity per host') + ylab('Average host range \n')+
  scale_color_gradient(low=beths[3], high=beths[7], name="Year") -> g2; g2

g1 + g2 + plot_layout(guides = 'collect') 

ggsave("~/Documents/Github/straightlinewasalie/Figures/virion-change-through-time.pdf", width = 8.5, height = 4)
