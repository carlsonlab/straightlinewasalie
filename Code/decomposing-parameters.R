
library(tidyverse)
library(patchwork)
library(MetBrewer)
library(MoMAColors)
library(vroom)

# setwd("~/Documents/Github/straightlinewasalie/Data")

# Access virus data
# library(virionData)
# get_versioned_data(version = "17636397", dir_path = "./Temp")

# helm <- read.csv("Data/helminths.csv")

# Set up virus data
virion <- vroom("./Data/Temp/17636397/virion.csv.gz")
virion %>%
  mutate(Host = str_to_lower(Host), Virus = str_to_lower(Virus)) %>%
  filter(DetectionMethod %in% c("PCR/Sequencing","Isolation/Observation")) %>%
  filter(HostClass=='mammalia') %>%
  filter(!(Host=='homo sapiens')) %>%
  filter(!is.na(Host)) %>% filter(!is.na(Virus)) %>%
  select(Host, Virus, VirusGenus, VirusFamily, ICTVRatified) %>%
  distinct() -> virion
virion %>%
  select(Host, Virus) %>%
  rename(Parasite = Virus) -> h

hosts <- unique(h$Host)
n.hosts <- length(hosts)
pars <- unique(h$Parasite)
n.pars <- length(pars)
specs <- table(h$Parasite) %>% data.frame() %>% filter(Freq == 1) %>% pull(Var1)
n.specs <- length(specs)
gens <- table(h$Parasite) %>% data.frame() %>% filter(Freq > 1) %>% pull(Var1)
n.gens <- length(gens)

### EXPERIMENTAL DATA

reps <- 10
for (j in 1:reps) {
  
  exp.j <- data.frame(iter = c(1:n.hosts),
                      h = NA,
                      p = NA,
                      phs = NA,
                      phg = NA,
                      hpg = NA)
  
  for (i in 1:n.hosts) {
    
    hosts.i <- sample(hosts, i)
    h.i <- h %>% filter(Host %in% hosts.i)
    h.spec.i <- h.i %>% filter(Parasite %in% specs)
    h.gen.i <- h.i %>% filter(Parasite %in% gens)
    
    exp.j$h[i] <- length(unique(h.i$Host))
    exp.j$p[i] <- length(unique(h.i$Parasite))
    exp.j$ps[i] <- length(unique(h.spec.i$Parasite))
    exp.j$pg[i] <- length(unique(h.gen.i$Parasite))
    exp.j$phs[i] <- exp.j$ps[i]/exp.j$h[i]
    exp.j$phg[i] <- mean(table(h.gen.i$Host))
    exp.j$hpg[i] <- mean(table(h.gen.i$Parasite))
    
    print(i)
    
  }
  
  if(j == 1) {
    exp <- exp.j
  } else {
    exp <- bind_rows(exp, exp.j)
  }
}



### GRAPHS

beths <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")

exp %>%
  ggplot(aes(x = h, y = p)) + 
  geom_point(size = 0.1, alpha = 0.5, pch = 16, color = beths[[4]]) + 
  geom_line(aes(x = h, y = p), col = 'black', data = est, size = 0.8, lty = 5) + 
  theme_bw() + xlab("Hosts") + ylab("Total viruses") -> g1; g1
exp %>%
  ggplot(aes(x = h, y = ps)) + 
  geom_point(size = 0.1, alpha = 0.5, pch = 16, color = beths[[7]]) +
  geom_line(aes(x = h, y = s), col = 'black', data = est, size = 0.8, lty = 5) + 
  theme_bw() + xlab("Hosts") + ylab("Specialist viruses")  -> g2; g2
exp %>%
  ggplot(aes(x = h, y = phs)) + 
  geom_point(size = 0.1, alpha = 0.5, pch = 16, color = beths[[7]]) + 
  geom_line(aes(x = h, y = phs), col = 'black', data = est, size = 0.8, lty = 5) + 
  theme_bw() + xlab("Hosts") + ylab("Specialist viruses per host")  -> g3; g3
exp %>%
  ggplot(aes(x = h, y = pg)) + 
  geom_point(size = 0.1, alpha = 0.5, pch = 16, color = beths[[6]]) +
  geom_line(aes(x = h, y = g), col = 'black', data = est, size = 0.8, lty = 5) + 
  theme_bw() + xlab("Hosts") + ylab("Generalist viruses")  -> g4; g4
exp %>%
  ggplot(aes(x = h, y = hpg)) + 
  geom_point(size = 0.1, alpha = 0.5, pch = 16, color = beths[[6]]) +
  geom_line(aes(x = h, y = hpg), col = 'black', data = est, size = 0.8, lty = 5) + 
  theme_bw() + xlab("Hosts") + ylab("Hosts per generalist virus")  -> g5; g5
exp %>%
  ggplot(aes(x = h, y = phg)) + 
  geom_point(size = 0.1, alpha = 0.5, pch = 16, color = beths[[6]]) + 
  geom_line(aes(x = h, y = phg), col = 'black', data = est, size = 0.8, lty = 5) + 
  theme_bw()+ xlab("Hosts") + ylab("Generalist viruses per host") -> g6; g6

(g1+g2+g3)/(g4+g5+g6) + plot_annotation(tag_levels = "A")

ggsave("./Figures/parameters.pdf", height = 6, width = 9.5)
ggsave("./Figures/parameters.jpg", height = 6, width = 9.5, dpi = 600)
