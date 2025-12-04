
#####################
### 0. Preamble #####
#####################

# Load packages

library(tidyverse)
library(codependent)
library(iNEXT)
library(patchwork)
library(MoMAColors)

# Load the data

library(patchwork)
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
  filter(ICTVRatified==TRUE) %>%
  rename(Parasite = Virus) %>%
  select(Host,Parasite) %>%
  distinct() %>%
  filter(!is.na(Host), !is.na(Parasite)) -> helm

# Summary variables and statistics

hosts <- unique(helm$Host)
pars <- unique(helm$Parasite)
H <- length(hosts)
P <- length(pars)

#####################
### 1.Sampling  #####
#####################

results <- data.frame(prop = NA,
                      est.1 = NA,
                      est.2 = NA,
                      est.3 = NA)
iter = 300

for (i in 1:iter){

prop.H = sample(c(5:95), 1)/100

H.i <- round(H*prop.H)
hosts.i <- sample(hosts, H.i)
helm.i <- helm %>% filter(Host %in% hosts.i)
pars.i <- unique(helm.i$Parasite)
P.i <- length(pars.i)
specs.i <- helm.i %>% count(Parasite) %>% filter(n==1) %>% pull(Parasite)
S.i <- length(specs.i)
gens.i <- helm.i %>% count(Parasite) %>% filter(n>1) %>% pull(Parasite)
G.i <- length(gens.i)

#####################
### 2. Methods  #####
#####################

###
### 2a. Linear 
###

est.1 <- (P.i/H.i)*H

###
### 2b. Chao
###

incidence <- helm.i %>% 
              mutate(Count = 1) %>% 
              pivot_wider(names_from = Host, values_from = Count) %>% 
              replace(is.na(.), 0) %>%
              select(-Parasite) %>% as.matrix()
est.2 <- estimateD(incidence, q = 0, datatype = "incidence_raw", level = H, base = 'size')$qD

###
### 2c. Linear plus Chao
###

est.spec <- (S.i/H.i)*H

incidence <- helm.i %>% 
  filter(Parasite %in% gens.i) %>%
  mutate(Count = 1) %>% 
  pivot_wider(names_from = Host, values_from = Count) %>% 
  replace(is.na(.), 0) %>%
  select(-Parasite) %>% as.matrix()
est.gen <- estimateD(incidence, q = 0, datatype = "incidence_raw", level = H, base = 'size')$qD

est.3 <- est.spec + est.gen

results[i,] <- c(prop.H, est.1, est.2, est.3)
print(results[i,])
}

###
### 3. Visualization
###

# write_csv(results, "~/Documents/Github/twothings/diversityestimates.csv")

beths <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")

results %>%
  pivot_longer(cols = c(est.1, est.2, est.3),
               names_to = 'method',
               values_to = 'estimate') %>%
  mutate(method = recode(method, !!!c('est.1' = 'Linear (all viruses)',
                                      'est.2' = 'iNEXT (all viruses)',
                                      'est.3' = 'Linear (specialists) + iNEXT (generalists)'))) %>%
  ggplot(aes(x = prop, y = estimate/P, group = method, color = method)) + 
  #geom_vline(xintercept = 0.16, linetype = 2) + 
  geom_hline(yintercept = 1, linetype = 'dashed') +   
  geom_point(alpha = 1, pch = 16, size = 2) + 
  ylim(0, 3) + 
  xlab('\n Proportion of hosts sampled') + ylab('\n Estimated / true virus diversity \n') + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.key.height = unit(-2, "cm")) + 
  guides(color = guide_legend(nrow = 4)) + labs(color = "Method") +
  scale_colour_manual(values = beths[c(10,6,3)])

ggsave("~/Documents/Github/straightlinewasalie/Figures/iNEXT.pdf", height = 7, width = 6) 
ggsave("~/Documents/Github/straightlinewasalie/Figures/iNEXT.jpg", height = 7, width = 6) 


 

