
setwd("~/Documents/Github/")

h <- read_csv("helminths/Data/associations cleaned v2.csv")
d <- read_csv("twothings/Data/Dobson.csv") %>%
  rename(PoulinPPH = DobsonPPH, PoulinHPP = DobsonHPP)

h %>% mutate(hostgroup = recode(hostgroup, !!!c("Chondrostei" = "Osteichthyes",
                                                "Cladistei" = "Osteichthyes",
                                                "Holostei" = "Osteichthyes",
                                                "Teleostei" = "Osteichthyes"))) -> h

d$NHMPPH <-NA
d$NHMHPP <-NA

for (i in 1:24) {
  h %>%
    filter(hostgroup == d$HostGroup[i]) %>%
    filter(group == d$ParGroup[i]) %>%
    select(Host, Parasite) %>% distinct() -> hi
  d$NHMHPP[i] <- hi %>% count(Host) %>% pull(n) %>% mean()
  d$NHMPPH[i] <- hi %>% count(Parasite) %>% pull(n) %>% mean()
}

d %>% ggplot(aes(x = PoulinPPH, y = NHMPPH)) + geom_point() + theme_bw(); cor(d$NHMPPH, d$PoulinPPH, use = "complete.obs")
d %>% ggplot(aes(x = PoulinHPP, y = NHMHPP)) + geom_point() + theme_bw(); cor(d$NHMHPP, d$PoulinHPP, use = "complete.obs")

d %>% 
  select(HostGroup, ParGroup, PoulinPPH, NHMPPH) %>%
  pivot_longer(cols = c("PoulinPPH", "NHMPPH"), names_to = "Study", values_to = "PPH") %>%
  mutate(Study = str_replace_all(Study, "PPH", "")) %>% 
  mutate(Study = factor(Study, levels = c("Poulin", "NHM"))) %>%
  unite(Group, c(HostGroup, ParGroup)) %>%
  ggplot(aes(x = Study, y = PPH, group = Group)) + geom_point() + geom_line() + theme_bw() + 
  ylab("Parasites per host") + 
  ggtitle(paste0("Correlation: ",  round(cor(d$NHMPPH, d$PoulinPPH, use = "complete.obs"), 3), ", p = ", cor.test(d$NHMPPH, d$PoulinPPH, use = "complete.obs")$p.value %>% signif(digits = 2))) -> d1 

d %>% 
  select(HostGroup, ParGroup, PoulinHPP, NHMHPP) %>%
  pivot_longer(cols = c("PoulinHPP", "NHMHPP"), names_to = "Study", values_to = "HPP") %>%
  mutate(Study = str_replace_all(Study, "HPP", "")) %>% 
  mutate(Study = factor(Study, levels = c("Poulin", "NHM"))) %>%
  unite(Group, c(HostGroup, ParGroup)) %>%
  ggplot(aes(x = Study, y = HPP, group = Group)) + geom_point() + geom_line() + theme_bw() + 
  ylab("Hosts per parasite") + 
  ggtitle(paste0("Correlation: ",  round(cor(d$NHMHPP, d$PoulinHPP, use = "complete.obs"), 3), ", p = ", cor.test(d$NHMHPP, d$PoulinHPP, use = "complete.obs")$p.value %>% signif(digits = 2))) -> d2

d %>% mutate(NHMRatio = NHMPPH/NHMHPP, PoulinRatio = PoulinPPH/PoulinHPP) -> d
d %>%
  select(HostGroup, ParGroup, PoulinRatio, NHMRatio) %>%
  pivot_longer(cols = c("PoulinRatio", "NHMRatio"), names_to = "Study", values_to = "Ratio") %>%
  mutate(Study = str_replace_all(Study, "Ratio", "")) %>% 
  mutate(Study = factor(Study, levels = c("Poulin", "NHM"))) %>%
  unite(Group, c(HostGroup, ParGroup)) %>%
  ggplot(aes(x = Study, y = Ratio, group = Group)) + geom_point() + geom_line() + theme_bw() + 
  ylab("Ratio") + 
  ggtitle(paste0("Correlation: ",  round(cor(d$NHMRatio, d$PoulinRatio, use = "complete.obs"), 3), ", p = ", cor.test(d$PoulinRatio, d$NHMRatio, use = "complete.obs")$p.value %>% signif(digits = 2))) -> d3; d3

library(patchwork)
d1+d2+d3
