
library(tidyverse)
library(patchwork)
library(igraph)
library(codependent)
library(bipartite)

setwd("~/Documents/Github/straightlinewasalie/Data/Edgelists/web-of-life_2025-11-16_214145")




# Organize network datasets

networks <- list.files()
networks <- setdiff(networks, c("README", "references.csv")) 
networks <- networks[!str_detect(networks, "FW")]

refs <- read_csv("references.csv") %>%
  select(ID, `Type of interactions`) %>%
  rename(names = ID, type = `Type of interactions`) -> refs


# Generate the data frame

network_data <- data.frame(names = networks)
network_data %>%
  mutate(names = str_replace(names, ".csv", "")) %>%
  left_join(refs) -> network_data

network_data$hosts <- NA
network_data$symbionts <- NA
network_data$specialists <- NA
network_data$generalists <- NA
network_data$sigma_s <- NA
network_data$sigma_g <- NA 
network_data$eta_g <- NA
network_data$eta_bar <- NA
network_data$modularity <- NA
network_data$modularity_generalists <- NA
network_data$nestedness <- NA
network_data$nestedness_generalists <- NA
network_data$connectance <- NA
network_data$connectance_generalists <- NA
network_data$host_sharing <- NA
network_data$z <- NA



# Start of for loop

for (i in 1:length(networks)) {

# Read and clean file
  
net <- read_csv(networks[[i]])
if(colnames(net)[[2]]=="Num. of hosts sampled") { net <- net %>% select(-"Num. of hosts sampled") }
if(colnames(net)[[2]]=="Abundance") { net <- net %>% select(-"Abundance") }
if(colnames(net)[[2]]=="Number of flowers") { net <- net %>% select(-"Number of flowers") }

# Network transforms

net_as_mat <- net
net_mat <- as.matrix(net_as_mat[,-1])
rownames(net_mat) <- net_as_mat[,1] %>% pull()
net_mat <- ifelse(net_mat > 0, 1, 0)

net %>%
  rename(Host = `...1`) %>%
  pivot_longer(-Host, names_to = "Symbiont", values_to = "Presence") %>%
  mutate(Presence = as.numeric(Presence > 0)) %>%
  filter(Presence == 1) %>% select(-Presence) %>% distinct() -> net

length(unique(net$Host)) -> hosts
length(unique(net$Symbiont)) -> symbionts

# Split specialists and generalists 

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

# Make a generalist network too for good measure

gen_mat <- net_mat[rownames(net_mat) %in% generalist_network$Host,
                   colnames(net_mat) %in% generalist_network$Symbiont,
                   drop = FALSE]# Summary statistics

# Start generating some statistics

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

# Network parameters
modularity <- modularity(cluster_walktrap(graph_from_edgelist(as.matrix(net))))
nestedness <- nested(net_mat, method = 'NODF2')
connectance <- networklevel(net_mat, index = "connectance")
host_sharing <- networklevel(net_mat, index = "mean number of shared partners")[[1]]
z <- tryCatch(codependent::copredict(as.data.frame(net), n.indep=100, iter=10, plot = FALSE)[[2]][2], error = function(e) {NA})

# Generalist network parameters, but, proofing this against the possibility a network is all specialists
if(length(gen_mat)==0) {
  modularity_generalists <- NA
  modularity_generalists <- NA
  modularity_generalists <- NA
} else {
  modularity_generalists <- modularity(cluster_walktrap(graph_from_edgelist(as.matrix(generalist_network))))
  nestedness_generalists <- nested(gen_mat, method = 'NODF2')
  connectance_generalists <- networklevel(gen_mat, index = "connectance")
}

# Write out the data
network_data[paste0(network_data$names,".csv") == networks[[i]],3:ncol(network_data)] <- c(hosts, symbionts, obs_specialists, obs_generalists, sigma_s, sigma_g, eta_g, eta_bar, modularity, modularity_generalists, nestedness, nestedness_generalists, connectance, connectance_generalists, host_sharing, z)
}

# Some final changes to the data frame

network_data$prop_s <- network_data$specialists/network_data$symbionts 
network_data$ratio_s <- network_data$specialists/network_data$generalists 

# End of data generation. Start of analysis





# Some parameter comparisons by ecological interaction type
network_data %>%
  ggplot(aes(x = type, y = sigma_s, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = type, y = sigma_g, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = type, y = eta_g, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = type, y = modularity, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = type, y = nestedness, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = type, y = connectance, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()

# Old patterns we wanted to reproduce and then explain to ourselves
network_data %>%
  ggplot(aes(x = z, y = modularity, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = z, y = host_sharing, group = type, color = type)) +
  geom_point()

# Codependency is just scaling stuff. Bam. Explained!
network_data %>%
  ggplot(aes(x = z, y = eta_bar, group = type, color = type)) + geom_point()
network_data %>%
  ggplot(aes(x = z, y = eta_g, group = type, color = type)) + geom_point()
network_data %>%
  ggplot(aes(x = z, y = sigma_s, group = type, color = type)) + geom_point()
network_data %>%
  ggplot(aes(x = z, y = sigma_g, group = type, color = type)) + geom_point
network_data %>%
  ggplot(aes(x = z, y = prop_s, group = type, color = type)) + geom_point()
network_data %>%
  ggplot(aes(x = z, y = ratio_s, group = type, color = type)) + geom_point()

# Are specialists spandrels? Maybe
network_data %>%
  ggplot(aes(x = prop_s, y = eta_g, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = modularity, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = nestedness, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = connectance, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = host_sharing, group = type, color = type)) +
  geom_point()

# What's going on in the generalist part of the network?
network_data %>%
  ggplot(aes(x = type, y = modularity_generalists, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = type, y = nestedness_generalists, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = modularity, y = modularity_generalists, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = nestedness, y = nestedness_generalists, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = connectance, y = connectance_generalists, group = type, color = type)) +
  geom_point()

# Are metrics on the generalists "better" somehow?
network_data %>%
  ggplot(aes(x = modularity, y = nestedness, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = modularity_generalists, y = nestedness_generalists, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = connectance, y = nestedness, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = connectance_generalists, y = nestedness_generalists, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = connectance, y = modularity, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = connectance_generalists, y = modularity_generalists, group = type, color = type)) +
  geom_point()






# A clean, nice figure about network constraints

network_subset <- network_data %>% filter(!(type %in% c("Plant-Ant")))

network_subset %>%
  ggplot(aes(x = modularity, y = modularity_generalists, group = type, color = type)) +
  geom_point() + 
  theme_bw() +
  xlim(0,1) + ylim(0,1)
network_subset %>%
  ggplot(aes(x = connectance, y = connectance_generalists, group = type, color = type)) +
  geom_point() + 
  theme_bw()+
  xlim(0,1) + ylim(0,1)
network_subset %>%
  ggplot(aes(x = nestedness, y = nestedness_generalists, group = type, color = type)) +
  geom_point() + 
  theme_bw()
