
library(tidyverse)
library(patchwork)
library(MetBrewer)

straightLine <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")

met <- c("#F4C40F", "#FE9B00", "#D8443C", "#DE597C", "#E87B89", "#9F5691", "#633372", "#1F6E9C", "#2B9B81", "#92C051")


clo1 <- read_csv("Data/CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv") %>%
  select(Host, Virus) %>%
  distinct() %>%
  filter(!Host=='homo sapiens')

clo2 <- read_csv("Data/helminths.csv") %>%
  rename(Virus = Parasite) %>%
  filter(!Host=='homo sapiens') %>% distinct()

for (cl in 1:2){
  if(cl==1){
    clo <- clo1
  } else {
    clo <- clo2
  }

Hosts <- length(unique(clo$Host))
Viruses <- length(unique(clo$Virus))

clo %>%
  select(Host, Virus) %>% distinct() %>%
  group_by(Virus) %>% summarize(HostRange = n()) %>%
  filter(HostRange==1) %>% pull(Virus) -> true_specialists

df1 <- data.frame(i = 1:10000, hosts = NA, false_specialists = NA)
for (i in 1:10000) {
  nhosts <- round(runif(1, min = 1, max = Hosts))
  hosts <- sample(unique(clo$Host), nhosts)
  clo_i <- clo %>% filter(Host %in% hosts)
  clo_i %>%
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>%
    filter(HostRange==1) %>% pull(Virus) -> obs_specialists
  obs_true_specialists <- obs_specialists[obs_specialists %in% true_specialists]
  obs_false_specialists <- obs_specialists[!(obs_specialists %in% true_specialists)]
  n_obs_true <- length(unique(obs_true_specialists))
  n_obs_false <- length(unique(obs_false_specialists))
  prop_false <- (n_obs_false)/(n_obs_false+n_obs_true)
  df1[i,] <- c(i, nhosts, prop_false)
  print(i)
}

df2 <- data.frame(i = 1:10000, hosts = NA, false_specialists = NA)
for (i in 1:10000) {
  nrows <- nrow(clo)
  nrow_i <- round(runif(1, min = 1, max = nrows))
  rows <- sample(c(1:nrows), nrow_i)
  clo_i <- clo[rows,]
  clo_i %>% pull(Host) %>% unique() %>% length() -> nhosts_i
  clo_i %>%
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>%
    filter(HostRange==1) %>% pull(Virus) -> obs_specialists
  obs_true_specialists <- obs_specialists[obs_specialists %in% true_specialists]
  obs_false_specialists <- obs_specialists[!(obs_specialists %in% true_specialists)]
  n_obs_true <- length(unique(obs_true_specialists))
  n_obs_false <- length(unique(obs_false_specialists))
  prop_false <- (n_obs_false)/(n_obs_false+n_obs_true)
  df2[i,] <- c(i, nhosts_i, prop_false)
  print(i)
}

thomas <- met.brewer("Thomas", n=8)

if(cl==1){
  df1$type = 'Idealized (complete)'
  df2$type = 'Realistic (incomplete)'
  df <- bind_rows(df1, df2)
  df %>%
    ggplot(aes(x = hosts, y = false_specialists, color = type, group = type)) +
    theme_bw() +
    geom_point(alpha = 0.05) +
    scale_color_manual(values = c(thomas[8], thomas[2])) +
    ggtitle('Mammalian viruses') + xlab("Hosts") + ylab("Proportion of false specialists") +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.position = c(0.25, 0.1), legend.title = element_blank()) -> g1
} else {
  df1$type = 'Idealized (complete)'
  df2$type = 'Realistic (incomplete)'
  df3 <- bind_rows(df1, df2)
  df3 %>%
    ggplot(aes(x = hosts, y = false_specialists, color = type, group = type)) +
    theme_bw() +
    geom_point(alpha = 0.05) +
    scale_color_manual(values = c(thomas[8], thomas[2])) +
    ggtitle('Vertebrate parasites') + xlab("Hosts") + ylab("Proportion of false specialists") +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.position = c(0.25, 0.1), legend.title = element_blank()) -> g2
}

}

