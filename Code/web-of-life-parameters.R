
library(tidyverse)
library(patchwork)
library(igraph)
library(codependent)

setwd("~/Documents/Github/twothings/Data/web-of-life_2025-11-16_214145")
networks <- list.files()
networks <- setdiff(networks, c("README", "references.csv")) 
networks <- networks[!str_detect(networks, "FW")]

refs <- read_csv("references.csv") %>%
  select(ID, `Type of interactions`) %>%
  rename(names = ID, type = `Type of interactions`) -> refs

network_data <- data.frame(names = networks)
network_data %>%
  mutate(names = str_replace(names, ".csv", "")) %>%
  left_join(refs) -> network_data

network_data$hosts <- NA
network_data$symbionts <- NA
network_data$sigma_s <- NA
network_data$sigma_g <- NA 
network_data$eta_g <- NA
network_data$eta_bar <- NA
network_data$modularity <- NA
network_data$host_sharing <- NA
network_data$z <- NA

for (i in 1:length(networks)) {

net <- read_csv(networks[[i]])

if(colnames(net)[[2]]=="Num. of hosts sampled") { net <- net %>% select(-"Num. of hosts sampled") }
if(colnames(net)[[2]]=="Abundance") { net <- net %>% select(-"Abundance") }
if(colnames(net)[[2]]=="Number of flowers") { net <- net %>% select(-"Number of flowers") }

net_as_mat <- net

net %>%
  rename(Host = `...1`) %>%
  pivot_longer(-Host, names_to = "Symbiont", values_to = "Presence") %>%
  mutate(Presence = as.numeric(Presence > 0)) %>%
  filter(Presence == 1) %>% select(-Presence) %>% distinct() -> net

length(unique(net$Host)) -> hosts
length(unique(net$Symbiont)) -> symbionts

net %>%
  group_by(Symbiont) %>%
  count() %>%
  filter(n==1) %>% 
  pull(Symbiont) -> specialists
length(specialists) -> obs_specialists

net %>%
  filter(Symbiont %in% specialists) -> specialist_network
net %>%
  filter(!(Symbiont %in% specialists)) -> generalist_network

generalist_network %>% 
  pull(Symbiont) %>% unique() %>% length() -> obs_generalists

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
  group_by(Symbiont) %>% 
  count() %>% pull(n) %>% mean() -> eta_g
net %>%
  group_by(Symbiont) %>% 
  count() %>% pull(n) %>% mean() -> eta_bar

modularity(cluster_walktrap(graph_from_edgelist(as.matrix(net)))) -> modularity 

net_mat <- as.matrix(net_as_mat[,-1])
rownames(net_mat) <- net_as_mat[,1] %>% pull()
host_sharing <- bipartite::networklevel(net_mat, index = "mean number of shared partners")[[1]]

z <- tryCatch(codependent::copredict(as.data.frame(net), n.indep=100, iter=10, plot = FALSE)[[2]][2], error = function(e) {NA})

network_data[paste0(network_data$names,".csv") == networks[[i]],3:11] <- c(hosts, symbionts, sigma_s, sigma_g, eta_g, eta_bar, modularity, host_sharing, z)
}

network_data %>% head()

network_data %>%
  ggplot(aes(x = type, y = eta_g, group = type, color = type)) + geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = z, y = host_sharing, group = type, color = type)) + geom_point()
network_data %>%
  ggplot(aes(x = z, y = eta_bar, group = type, color = type)) + geom_point()

network_data %>%
  ggplot(aes(x = z, y = eta_bar, group = type, color = type)) + geom_point()
