
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
  nestedness_generalists <- tryCatch(nested(gen_mat, method = 'NODF2'), error = function(e) {NA})
  connectance_generalists <- tryCatch(networklevel(gen_mat, index = "connectance"), warning = function(e) {NA})
}

# Write out the data
network_data[paste0(network_data$names,".csv") == networks[[i]],3:ncol(network_data)] <- c(hosts, symbionts, obs_specialists, obs_generalists, sigma_s, sigma_g, eta_g, eta_bar, modularity, modularity_generalists, nestedness, nestedness_generalists, connectance, connectance_generalists, host_sharing, z)
}

# Some final changes to the data frame:
# No tiny tiny networks
network_data <- network_data %>% filter(hosts >= 5, symbionts >= 5)
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
network_data %>%
  ggplot(aes(x = type, y = prop_s, group = type, color = type)) + 
  geom_boxplot() + geom_jitter()
network_data %>%
  ggplot(aes(x = type, y = ratio_s, group = type, color = type)) + 
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
  ggplot(aes(y = z, x = eta_bar, group = type, color = type)) + geom_point() + scale_x_log10()
network_data %>%
  ggplot(aes(y = z, x = eta_g, group = type, color = type)) + geom_point() + scale_x_log10()
network_data %>%
  ggplot(aes(y = z, x = sigma_s, group = type, color = type)) + geom_point() + scale_x_log10()
network_data %>%
  ggplot(aes(y = z, x = sigma_g, group = type, color = type)) + geom_point() + scale_x_log10()
network_data %>%
  ggplot(aes(y = z, x = prop_s, group = type, color = type)) + geom_point()
network_data %>%
  ggplot(aes(y = z, x = ratio_s, group = type, color = type)) + geom_point()

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

# Proportion specialists correlates with generalist network attributes but doesn't constrain them
network_data %>%
  ggplot(aes(x = prop_s, y = modularity_generalists, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = nestedness_generalists, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = connectance_generalists, group = type, color = type)) +
  geom_point()

# Proportion specialists correlates with generalist network attributes but doesn't constrain them
network_data %>%
  ggplot(aes(x = prop_s, y = eta_bar, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = eta_g, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = sigma_s, group = type, color = type)) +
  geom_point()
network_data %>%
  ggplot(aes(x = prop_s, y = sigma_g, group = type, color = type)) +
  geom_point()


# Quick check: are constraints real?
network_data %>%
  ggplot(aes(x = prop_s, y = eta_g, group = type, color = log(hosts))) +
  geom_point(size = 2) + theme_bw() +
  scale_color_gradient2(low='red', mid='yellow', high='blue')






# A clean, nice figure about network constraints

table(network_data$type)
network_subset <- network_data %>% 
  filter(!(type %in% c("Plant-Ant", "Plant-Herbivore","Anemone-Fish"))) %>%
  mutate(type = recode(type, !!!c("Pollination" = "Plant-Pollinator", "Seed Dispersal" = "Plant-Seed Disperser"))) %>%
  mutate(type = factor(type, levels = c("Host-Parasite","Plant-Seed Disperser","Plant-Pollinator"))) %>%
  mutate(nestedness = nestedness/100, nestedness_generalists = nestedness_generalists/100) %>%
  arrange(desc(type))

beths <- c("#DD7188","#6F8844","#90C2DB","#8B6654","#F4D337","#5AB8DA","#AC3F3E")
beth0 <- beths[c(1,5,6)]

network_subset %>%
  ggplot(aes(x = prop_s, y = modularity)) + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Proportion of specialists") + ylab("Modularity (full network)") +
  geom_smooth(method = 'lm', color = 'gray20', fill = 'black') +
  scale_fill_manual(values = beth0, name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) -> scatter_prop_m

network_subset %>%
  ggplot(aes(x = prop_s, y = connectance)) + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Proportion of specialists") + ylab("Connectance (full network)") +
  geom_smooth(method = 'lm', color = 'gray20', fill = 'black') +
  scale_fill_manual(values = beth0, name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) -> scatter_prop_c
  
network_subset %>%
  ggplot(aes(x = prop_s, y = nestedness)) + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Proportion of specialists") + ylab("Nestedness (full network)") +
  geom_smooth(method = 'lm', color = 'gray20', fill = 'black') +
  scale_fill_manual(values = beth0, name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) -> scatter_prop_n

scatter_prop_c + scatter_prop_m + scatter_prop_n -> scatter_prop



network_subset %>%
  ggplot(aes(x = modularity, y = modularity_generalists, group = type, color = type)) +
  geom_segment(aes(x = 0, xend = 0.7, y = 0, yend = 0.7), color = 'gray20', lwd = 1, lineend = "round", linetype = 'dashed') + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Modularity (full network)") + ylab("Modularity (generalists)") +
  xlim(0,0.7) + ylim(0,0.7) +
  scale_fill_manual(values = beth0, name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) -> scatter_m

network_subset %>%
  ggplot(aes(x = connectance, y = connectance_generalists, group = type, color = type)) +
  geom_segment(aes(x = 0, xend = 0.75, y = 0, yend = 0.75), color = 'gray20', lwd = 1, lineend = "round", linetype = 'dashed') + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Connectance (full network)") + ylab("Connectance (generalists)") +
  xlim(0,0.75) + ylim(0,0.75) +
  scale_fill_manual(values = beth0, name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) -> scatter_c

network_subset %>%
  ggplot(aes(x = nestedness, y = nestedness_generalists, group = type, color = type)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = 'gray20', lwd = 1, lineend = "round", linetype = 'dashed') + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() + 
  xlab("Nestedness (full network)") + ylab("Nestedness (generalists)") +
  xlim(0,1) + ylim(0,1) +
  scale_fill_manual(values = beth0, name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) -> scatter_n

scatter <- scatter_c + scatter_m + scatter_n



network_subset %>%
  select(type, modularity, modularity_generalists) %>% 
  pivot_longer(cols = c(modularity, modularity_generalists)) %>% 
  mutate(name = recode(name, !!!c("modularity" = "Full network", "modularity_generalists" = "Generalists"))) %>%

  ggplot(aes(x = type, y = value)) +
  geom_point(aes(group = name, fill = type), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.8, shape = 21, alpha = 0.25, color = 'black') + 
  geom_boxplot(aes(fill = type, alpha = name), outliers = FALSE, lwd = 0.4) + 
  theme_bw() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()) + 
  scale_alpha_manual(values = c(0.8, 0.4), guide = guide_legend(override.aes = list(fill = "black")), name = NULL) + 
  xlab("Type of symbiosis") + ylab("Modularity") + 
  guides(fill = 'none', color = 'none') + 
  scale_fill_manual(values = beth0, name = NULL) +
  
  ylim(0,1) -> bars_m

network_subset %>%
  select(type, connectance, connectance_generalists) %>% 
  pivot_longer(cols = c(connectance, connectance_generalists)) %>% 
  mutate(name = recode(name, !!!c("connectance" = "Full network", "connectance_generalists" = "Generalists"))) %>%

  ggplot(aes(x = type, y = value)) +
  geom_point(aes(group = name, fill = type), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.8, shape = 21, alpha = 0.25, color = 'black') + 
  geom_boxplot(aes(fill = type, alpha = name), outliers = FALSE, lwd = 0.4) + 
  theme_bw() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()) + 
  scale_alpha_manual(values = c(0.8, 0.4), guide = guide_legend(override.aes = list(fill = "black")), name = NULL) + 
  xlab("Type of symbiosis") + ylab("Connectance") + 
  guides(fill = 'none', color = 'none') + 
  scale_fill_manual(values = beth0, name = NULL) +
  
  ylim(0,1) -> bars_c

network_subset %>%
  select(type, nestedness, nestedness_generalists) %>% 
  pivot_longer(cols = c(nestedness, nestedness_generalists)) %>% 
  mutate(name = recode(name, !!!c("nestedness" = "Full network", "nestedness_generalists" = "Generalists"))) %>%
  
  ggplot(aes(x = type, y = value)) +
  geom_point(aes(group = name, fill = type), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.8, shape = 21, alpha = 0.25, color = 'black') + 
  geom_boxplot(aes(fill = type, alpha = name), outliers = FALSE, lwd = 0.4) + 
  theme_bw() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank()) + 
  scale_alpha_manual(values = c(0.8, 0.4), guide = guide_legend(override.aes = list(fill = "gray0")), name = NULL) + 
  xlab("Type of symbiosis") + ylab("Nestedness") +  
  guides(fill = 'none', color = 'none') +
  scale_fill_manual(values = beth0, name = NULL) + 
  
  ylim(0,1) -> bars_n

bars <- bars_c + bars_m + bars_n

((scatter_prop_c + scatter_prop_m + scatter_prop_n)/(scatter_c + scatter_m + scatter_n)/(bars_c + bars_m + bars_n)) + 
  plot_layout(guides = 'collect') & theme(legend.position = 'top', legend.text = element_text(size = 11), legend.key.height = unit(1, 'cm'))

ggsave(filename = "C:/Users/cjc277/OneDrive - Yale University/Documents/Github/straightlinewasalie/Figures/constraints.pdf", 
       width = 10, height = 10, units = "in")





# A clean, nice figure about power laws

network_subset_2 <- network_data %>% 
  mutate(type = recode(type, !!!c("Pollination" = "Plant-Pollinator", 
                                  "Seed Dispersal" = "Plant-Seed Disperser")))

network_subset_2 %>%
  ggplot(aes(y = z, x = prop_s, color = type)) + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Proportion of specialists") + ylab("Power law exponent (z)") +
  geom_smooth(method = 'lm', color = 'gray20', fill = 'black') +
  scale_fill_manual(values = beths[c(4,1,7,2,6,5)], name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) -> p1
network_subset_2 %>%
  ggplot(aes(y = z, x = sigma_s, color = type)) + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Specialist symbionts per host (log10)") + ylab("Power law exponent (z)") +
  geom_smooth(method = 'lm', color = 'gray20', fill = 'black') +
  scale_fill_manual(values = beths[c(4,1,7,2,6,5)], name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) +
  scale_x_log10() -> p2
network_subset_2 %>%
  ggplot(aes(y = z, x = sigma_g, color = type)) + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Generalist symbionts per host (log10)") + ylab("Power law exponent (z)") +
  geom_smooth(method = 'lm', color = 'gray20', fill = 'black') +
  scale_fill_manual(values = beths[c(4,1,7,2,6,5)], name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) +
  scale_x_log10() -> p3
network_subset_2 %>%
  ggplot(aes(y = z, x = eta_g, color = type)) + 
  geom_point(aes(fill = type), size = 2, stroke = 0.5, shape = 21, alpha = 0.7, color = 'gray20') + 
  theme_bw() +
  xlab("Average host range of generalists (log10)") + ylab("Power law exponent (z)") +
  geom_smooth(method = 'lm', color = 'gray20', fill = 'black') +
  scale_fill_manual(values = beths[c(4,1,7,2,6,5)], name = NULL) + 
  guides(fill = guide_legend(override.aes = list(size=5))) +
  scale_x_log10() -> p4

((p1+p2)/(p3+p4)) + plot_layout(guides='collect') & theme(legend.position = 'bottom', legend.text = element_text(size = 9.5))


ggsave(filename = "~/Documents/Github/straightlinewasalie/Figures/power-law-exponents.pdf", 
       width = 10, height = 11, units = "in")
