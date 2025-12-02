
library(tidyverse)
library(patchwork)

clo <- read_csv("~/Documents/Github/clover/clover/clover_0.1_mammalviruses/CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv")

clo %>% select(Host, Virus, PublicationYear) %>%
  distinct() %>% 
  filter(!is.na(PublicationYear)) %>%
  filter(!Host=='homo sapiens') %>%
  group_by(Host, Virus) %>%
  summarize(PublicationYear = min(PublicationYear)) -> clo

clo %>% pull(Host) %>% unique() %>% length() -> total_hosts
clo %>% pull(Virus) %>% unique() %>% length() -> total_viruses

clo %>%
  group_by(Virus) %>%
  count() -> hostrange 
# hostrange %>% ggplot(aes(x = n)) + geom_histogram() + scale_x_log10()

hostrange %>%
  filter(n==1) %>% 
  pull(Virus) -> specialists

length(specialists) -> total_specialists
total_viruses - total_specialists -> total_generalists

params <- data.frame(
  year = c(1960:2017),
  viruses = NA, hosts = NA, sigma_s = NA, sigma_g = NA, eta_g = NA, b = NA, true_specialists = NA, obs_specialists = NA, obs_generalists = NA, prop_false_specialists = NA
)

for(year in 1960:2017) {

clo %>% filter(PublicationYear <= year) -> clo_year
  
clo_year %>%
    pull(Virus) %>% unique() %>% length() -> viruses
clo_year %>%
    pull(Host) %>% unique() %>% length() -> hosts
  
clo_year %>%
  group_by(Virus) %>%
  count() %>%
  filter(n==1) %>% 
  pull(Virus) -> specialists_year
length(specialists_year) -> obs_specialists

clo_year %>%
  filter(Virus %in% specialists_year) -> specialist_network
clo_year %>%
  filter(!(Virus %in% specialists_year)) -> generalist_network

generalist_network %>% 
  pull(Virus) %>% unique() %>% length() -> obs_generalists

specialist_network %>% 
  pull(Host) %>% unique() %>% length() -> hosts_of_specialists
generalist_network %>% 
  pull(Host) %>% unique() %>% length() -> hosts_of_generalists

specialist_network %>%
  group_by(Host) %>% 
  count() %>% pull(n) %>% mean() -> sigma_s
sigma_s*(hosts_of_specialists/hosts) -> sigma_s
generalist_network %>%
  group_by(Host) %>% 
  count() %>% pull(n) %>% mean() -> sigma_g
sigma_g*(hosts_of_generalists/hosts) -> sigma_g
generalist_network %>%
  group_by(Virus) %>% 
  count() %>% pull(n) %>% mean() -> eta_g
(eta_g - 1)/hosts -> b

specialist_network %>% 
  ungroup() %>% select(Virus) %>% distinct() %>%
  mutate(true_specialist = (Virus %in% specialists)) %>%
  pull(true_specialist) %>% as.numeric() %>% sum() -> true_specialists
1-(true_specialists/obs_specialists) -> prop_false_specialists

params[params$year == year,] <- c(year, viruses, hosts, sigma_s, sigma_g, eta_g, b, true_specialists, obs_specialists, obs_generalists, prop_false_specialists)

}

params %>%
  mutate(scaling_generalists = (sigma_g)/(1+b*hosts)) %>%
  mutate(scaling_specialists = sigma_s) %>% 
  mutate(est_specialists = scaling_specialists*total_hosts,
         est_generalists = scaling_generalists*total_hosts,
         est_total = est_specialists+est_generalists) -> params

params %>%
  select(year, hosts, viruses) %>% pivot_longer(cols = c(hosts, viruses)) %>%
  ggplot(aes(x = year, y = value, group = name, col = name)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank(), legend.position = c(0.18, 0.85)) +
  ggtitle('Total species richness') + ylim(0, max(params$hosts)*1.05) -> g1

params %>%
  select(year, obs_generalists, obs_specialists, true_specialists) %>% pivot_longer(cols = c(obs_generalists, obs_specialists, true_specialists)) %>%
  ggplot(aes(x = year, y = value, group = name, col = name)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank(), legend.position = c(0.18, 0.85)) +
  ggtitle('Virus species richness') + ylim(0, max(params$obs_generalists)*1.05) -> g2

params %>%
  select(year, sigma_s, sigma_g) %>% pivot_longer(cols = c(sigma_s, sigma_g)) %>%
  ggplot(aes(x = year, y = value, group = name, col = name)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank(), legend.position = c(0.18, 0.85)) +
  ggtitle('Viruses per individual host species') + ylim(0, max(sigma_g)*1.05)  -> g3

params %>%
  select(year, eta_g) %>% 
  ggplot(aes(x = year, y = eta_g)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank()) +
  ggtitle('Host range of observed generalist viruses') + ylim(0, max(params$eta_g)*1.05)  -> g4

params %>%
  select(year, scaling_specialists, scaling_generalists) %>% pivot_longer(cols = c(scaling_specialists, scaling_generalists)) %>%
  ggplot(aes(x = year, y = value, group = name, col = name)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank(), legend.position = c(0.18, 0.85)) +
  ggtitle('Viruses per total host species') + ylim(0, max(params$scaling_specialists)*1.05)  -> g5 

params %>%
  select(year, prop_false_specialists) %>% 
  ggplot(aes(x = year, y = prop_false_specialists)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank()) +
  ggtitle('Proportion of false specialists') + ylim(0, 1) -> g6

g1 + g2 + g3 + g4 + g5 + g6

params %>%
  ggplot(aes(x = year, y = est_specialists)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank()) +
  geom_hline(yintercept = total_specialists) + 
  ggtitle('Estimated total specialist viruses') + ylim(0, max(params$est_specialists)*1.05) -> g7

params %>%
  ggplot(aes(x = year, y = est_generalists)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank()) +
  geom_hline(yintercept = total_generalists) + 
  ggtitle('Estimated total generalist viruses') + ylim(0, max(params$est_generalists)*1.05) -> g8

params %>%
  ggplot(aes(x = year, y = est_total)) +
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + theme(legend.title = element_blank()) +
  geom_hline(yintercept = total_viruses) + 
  ggtitle('Estimated total viruses') + ylim(0, max(params$est_total)*1.05) -> g9

g1 + g2 + g6 + g3 + g4 + g5 + g7 + g8 + g9


params %>%
  ggplot(aes(x = hosts, y = viruses, color = year)) + 
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + xlim(0, total_hosts*1.05) + ylim(0, total_viruses*1.05)  -> g

params %>%
  ggplot(aes(x = hosts, y = viruses, color = year)) + 
  theme_bw() + geom_point() + ylab(NULL) + xlab(NULL) + xlim(0, total_hosts*1.05) + ylim(0, total_viruses*1.05) 
for (year in 1960:2017){
  slope <- params[params$year == year,]$viruses/params[params$year == year,]$hosts
  year_fake <- data.frame(year = year, hosts = c(1:(total_hosts*1.5)), viruses = slope*c(1:(total_hosts*1.5)))
  if(year == 1960) {all_years_fake <- year_fake} else {all_years_fake <- rbind(all_years_fake, year_fake)}
}
g + geom_line(data = all_years_fake, aes(x = hosts, y = viruses, group = year, color = year), alpha = 0.1)


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# 
# library(tidyverse)
# library(patchwork)
# 
# clo <- read_csv("~/Documents/Github/clover/clover/clover_0.1_mammalviruses/CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv")
# 
# clo %>% select(Host, Virus, PublicationYear) %>%
#   distinct() %>% 
#   filter(!is.na(PublicationYear)) %>%
#   filter(!Host=='homo sapiens') %>%
#   group_by(Host, Virus) %>%
#   summarize(PublicationYear = min(PublicationYear)) -> clo
# 
# clo %>% pull(Host) %>% unique() %>% length() -> total_hosts
# clo %>% pull(Virus) %>% unique() %>% length() -> total_viruses
# 
# params <- data.frame(
#   year = c(1960:2017),
#   viruses = NA, hosts = NA, 
#   sigma = NA, eta = NA
# )
# 
# for(year in 1960:2017) {
#   
#   clo %>% filter(PublicationYear <= year) -> clo_year
#   
#   clo_year %>%
#     pull(Virus) %>% unique() %>% length() -> viruses
#   clo_year %>%
#     pull(Host) %>% unique() %>% length() -> hosts
#   
#   clo_year %>%
#     group_by(Virus) %>%
#     count() %>% pull(n) %>% mean() -> eta 
#   
#   clo_year %>%
#     group_by(Host) %>%
#     count() %>% pull(n) %>% mean() -> sigma  
#   
#   params[params$year == year,] <- c(year, viruses, hosts, sigma, eta)
#   
# }
# 
# params %>%
#   ggplot(aes(x = hosts, y = viruses, color = year)) + 
#   geom_point(size = 3) + 
#   theme_bw() +
#   xlab('\n Hosts') + ylab('Viruses \n') -> g1
# 
# params %>%
#   ggplot(aes(x = sigma, y = eta, color = year)) + 
#   geom_point(size = 3) + 
#   theme_bw() +
#   xlab('\n Average virus diversity per host') + ylab('Average host range \n') -> g2
# 
# g1 + g2 + plot_layout(guides = 'collect')
library(tidyverse)
library(patchwork)

clo <- read_csv("~/Documents/Github/clover/clover/clover_0.1_mammalviruses/CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv")

clo %>% select(Host, Virus, PublicationYear) %>%
  distinct() %>% 
  filter(!is.na(PublicationYear)) %>%
  filter(!Host=='homo sapiens') %>%
  group_by(Host, Virus) %>%
  summarize(PublicationYear = min(PublicationYear)) -> clo

clo %>% pull(Host) %>% unique() %>% length() -> total_hosts
clo %>% pull(Virus) %>% unique() %>% length() -> total_viruses

params <- data.frame(
  year = c(1960:2017),
  viruses = NA, hosts = NA, 
  sigma = NA, eta = NA
)

for(year in 1960:2017) {
  
  clo %>% filter(PublicationYear <= year) -> clo_year
  
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

params %>%
  ggplot(aes(x = hosts, y = viruses, color = year)) + 
  geom_point(size = 3) + 
  theme_bw() +
  xlab('\n Hosts') + ylab('Viruses \n') -> g1

params %>%
  ggplot(aes(x = sigma, y = eta, color = year)) + 
  geom_point(size = 3) + 
  theme_bw() +
  xlab('\n Average virus diversity per host') + ylab('Average host range \n') -> g2

g1 + g2 + plot_layout(guides = 'collect')


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

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

params %>%
  ggplot(aes(x = hosts, y = viruses, color = year)) + 
  geom_point(size = 3) + 
  theme_bw() +
  xlab('\n Hosts') + ylab('Viruses \n') +
  scale_color_gradientn(colors=moma.colors("OKeeffe", direction = 1), name = 'Year') -> g1; g1

params %>%
  ggplot(aes(x = sigma, y = eta, color = year)) + 
  geom_point(size = 3) + 
  theme_bw() +
  xlab('\n Average virus diversity per host') + ylab('Average host range \n')+
  scale_color_gradientn(colors=moma.colors("OKeeffe", direction = 1), name = 'Year') -> g2

g1 + g2 + plot_layout(guides = 'collect') 
  

