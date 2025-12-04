
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

results %>%
  pivot_longer(cols = c(est.1, est.2, est.3),
               names_to = 'method',
               values_to = 'estimate') %>%
  mutate(method = recode(method, !!!c('est.1' = 'Linear (all symbionts)',
                                      'est.2' = 'iNEXT (all symbionts)',
                                      'est.3' = 'Linear (specialists) + iNEXT (generalists)'))) %>%
  ggplot(aes(x = prop, y = estimate/P, group = method, color = method)) + 
  geom_vline(xintercept = 0.16, linetype = 2) + 
  geom_hline(yintercept = 1) +   
  geom_point(alpha = 0.35, stroke = NA, size = 2.5) + 
  scale_y_log10(limits = c(0.1, 10)) +   
  annotate("rect", xmin = 0.05, xmax = 0.15, 
           ymin = 0.3, ymax = 16, 
           alpha = .25) + 
  xlab('\n Proportion of hosts sampled') + ylab('\n Estimated parasite diversity / true parasite diversity \n') + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.key.height = unit(-2, "cm")) + 
  guides(color = guide_legend(nrow = 4)) + labs(color = "Method") +
  scale_colour_manual(values = moma.colors("Panton", 7)[c(3,5,4)])
 


 

