
library(tidyverse)
library(patchwork)
library(MetBrewer)
library(MoMAColors)
library(vroom)

# setwd("~/Documents/Github/straightlinewasalie/Data")

# Access virus data
# library(virionData)
# get_versioned_data(version = "17636397", dir_path = "./Temp")

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

# Read in metagenomic data
macaca <- read_csv("./Data/Metagenomics/macaca_estimates_by_viral_family.csv")
ptero <- read_csv("./Data/Metagenomics/pteropus_estimates_by_viral_family.csv")

# Set up a structure in which to slot things

macaca %>%
  select(Site, Observed) %>%
  mutate(Site = recode(Site, !!!c("MmAdV" = "adenoviridae",
                                  "MmBoV" = "parvoviridae",
                                  "MmCoV" = "coronaviridae",
                                  "MmPicorna" = "picornaviridae",
                                  "MmPosaV" = "posavirus",
                                  "MmAaV" = "parvoviridae",
                                  "MmHV" = 'orthoherpesviridae',
                                  "MmPyV" = "polyomaviridae",
                                  "MmRbV" = "rhabdoviridae",
                                  "MmAstV" = "astroviridae",
                                  "MmPbV" = "picobirnaviridae",
                                  "MmRota" = "sedoreoviridae"))) %>%
  group_by(Site) %>% summarize(Observed = sum(Observed)) %>%
  rename(VirusFamily = Site, EstM = Observed) -> macaca

ptero %>%
  select(Site, Observed) %>%
  mutate(Site = recode(Site, !!!c("AdV" = "adenoviridae",
                                  "BoV" = "parvoviridae",
                                  "CoV" = "coronaviridae",
                                  "HV" = 'orthoherpesviridae',
                                  "PyV" = "polyomaviridae",
                                  "AstV" = "astroviridae",
                                  "PmV" = "paramyxoviridae"))) %>%
  group_by(Site) %>% summarize(Observed = sum(Observed)) %>%
  rename(VirusFamily = Site, EstP = Observed) -> ptero

# Excluding posaviruses because their rank is unknown
#
# adeno-associated viruses (AaV) and bocaparvovirus (BoV) are both
# in the parvoviridae
#
# rotavirus is in the sedoreoviridae

virus_groups <-
  data.frame(VirusFamily = c("adenoviridae","coronaviridae","picornaviridae","orthoherpesviridae","polyomaviridae","rhabdoviridae","astroviridae","picobirnaviridae","sedoreoviridae","paramyxoviridae","parvoviridae")) %>%
  arrange(VirusFamily)

virus_groups %>%
  left_join(macaca) %>%
  left_join(ptero) %>%
  mutate(EstM = replace_na(EstM, 0)) %>%
  mutate(EstP = replace_na(EstP, 0)) %>%
  mutate(Est = (EstM + EstP)/2) %>%
  select(-EstM) %>% select(-EstP) %>%
  rename(VirusesPerHost = Est) -> virus_groups

virus_groups %>%
  mutate(HostRange = NA,
         HostRangeICTV = NA) -> virus_groups

# Get host range out of VIRION

for (i in 1:nrow(virus_groups)){
  virion %>%
    filter(VirusFamily == virus_groups$VirusFamily[i]) %>%
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>% pull(HostRange) %>% mean() -> HostRange
  virus_groups$HostRange[i] <- HostRange
  virion %>%
    filter(ICTVRatified == TRUE) %>%
    filter(VirusFamily == virus_groups$VirusFamily[i]) %>%
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>% pull(HostRange) %>% mean() -> HostRange
  virus_groups$HostRangeICTV[i] <- HostRange
}

virus_groups %>%
  mutate(Scaling = VirusesPerHost/HostRange,
         ScalingICTV = VirusesPerHost/HostRangeICTV) -> virus_groups

MammalSpecies = 6759
VirusFamilies = 25
est0 <- MammalSpecies*VirusFamilies*mean(virus_groups$VirusesPerHost)
est1 <- MammalSpecies*VirusFamilies*mean(virus_groups$Scaling)
est2 <- MammalSpecies*VirusFamilies*mean(virus_groups$ScalingICTV)

virus_groups %>%
  rename(Family = 'VirusFamily') %>%
  mutate(Family = str_to_sentence(Family)) -> virus_groups

virus_groups %>%
  mutate_at(c('VirusesPerHost', 'HostRange','HostRangeICTV','Scaling','ScalingICTV'), function(x){round(x,digits = 1)})


#### graphs


set.seed(721)

straightLine <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")

met <- c("#F4C40F", "#FE9B00", "#D8443C", "#DE597C", "#E87B89", "#9F5691", "#633372", "#1F6E9C", "#2B9B81", "#92C051")

cols <- c(met, 'gray90')

virus_groups_nopico <- virus_groups
virus_groups_nopico$ScalingICTV[virus_groups$Family=='Picobirnaviridae'] <- 0

virus_groups$estVirus <- virus_groups$VirusesPerHost
virus_groups_nopico$estVirus <- virus_groups_nopico$ScalingICTV

virus_groups_plot <- rbind(
  data.frame(virus_groups, picoOut="Carroll et al. (2018)\nRaw metagenomic samples\nuncorrected for host range"),
  data.frame(virus_groups_nopico, picoOut="This study (2025)\nGlobal estimates broken down\nfor vertebrate-infective groups")
  ) %>%
  mutate(`Virus family` = factor(Family, levels = c("Adenoviridae","Astroviridae","Coronaviridae","Orthoherpesviridae","Paramyxoviridae","Parvoviridae","Picornaviridae","Polyomaviridae","Rhabdoviridae","Sedoreoviridae","Picobirnaviridae")))


pies <- ggplot(virus_groups_plot, aes(x="", y=estVirus, fill=`Virus family`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=cols) +
  facet_wrap("picoOut", nrow=1, scale="free") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5),
      text = element_text(size = 10),
      plot.subtitle = element_text(hjust=0.5),
      strip.text = element_text(size=10, margin=margin(0.01, 0.01, 0.05, 0.01, unit="cm")),
      panel.spacing = unit(-1, "cm"),
      legend.key.size = unit(0.5, "cm"))

# pies

ggplot(aes(x="", y=VirusesPerHost, fill=`Virus family`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + theme_void() +
  ggtitle("Carroll et al. (2018)", subtitle="Raw metagenomic samples \nuncorrected for host range") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        plot.subtitle = element_text(hjust=0.5)) +
  scale_fill_manual(values=cols) -> g1
#virus_groups %>%
#  ggplot(aes(x="", y=ScalingICTV, fill=Family)) +
#  geom_bar(stat="identity", width=1) +
#  coord_polar("y", start=0) + theme_void()  +
#  ggtitle("Corrected for host range") +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  scale_fill_manual(values=cols) -> g2

virus_groups_nopico %>%
  mutate(`Virus family` = factor(Family, levels = c("Adenoviridae","Astroviridae","Coronaviridae","Orthoherpesviridae","Paramyxoviridae","Parvoviridae","Picornaviridae","Polyomaviridae","Rhabdoviridae","Sedoreoviridae","Picobirnaviridae"))) %>%
  ggplot(aes(x="", y=ScalingICTV, fill=`Virus family`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + theme_void()  +
  ggtitle("This study (2025)", subtitle="Global estimates broken down\nfor vertebrate-infective groups") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        plot.subtitle = element_text(hjust=0.5)) +
  scale_fill_manual(values=cols) -> g3
#g1 + g2 + g3 + plot_layout(guides = 'collect')

### no picobirnas

virus_groups_nopb <- virus_groups %>% filter(!(Family == 'Picobirnaviridae'))
est0 <- MammalSpecies*VirusFamilies*mean(virus_groups_nopb$VirusesPerHost)
est1 <- MammalSpecies*VirusFamilies*mean(virus_groups_nopb$Scaling)
est2 <- MammalSpecies*VirusFamilies*mean(virus_groups_nopb$ScalingICTV)

virus_groups_nopb %>%
  mutate(Baltimore = c("DNA","RNA","RNA","DNA","RNA","DNA","RNA","DNA","RNA","RNA")) -> split_groups

split_groups %>% filter(Baltimore == "DNA") -> dna_viruses
split_groups %>% filter(Baltimore == "RNA") -> rna_viruses

est1dna <- MammalSpecies*VirusFamilies*mean(dna_viruses$Scaling)
est2dna <- MammalSpecies*VirusFamilies*mean(dna_viruses$ScalingICTV)
est1rna <- MammalSpecies*VirusFamilies*mean(rna_viruses$Scaling)
est2rna <- MammalSpecies*VirusFamilies*mean(rna_viruses$ScalingICTV)
est1 <- est1dna+est1rna
est2 <- est2dna+est2rna

############## Comparisons!

library(tidyverse)

df <- data.frame(Study = c("Carlson et al.\n(2019)",
                           "This study \n(conservative)",
                           "This study \n(liberal)",
                           "Carroll et al.\n(2018)"),
                 Diversity = c(40878,
                               est2,
                               est1,
                               1531745),
                 Lower = c(33357,
                           NA,
                           NA,
                           640211),
                 Upper = c(74354,
                           NA,
                           NA,
                           2423278))

df %>%
  mutate(Study = factor(Study, levels = df$Study)) %>%
  ggplot(aes(x = Study, y = Diversity/1000, fill = Study)) +
  geom_bar(stat="identity",
           position=position_dodge(),
           width = 0.7) +
  #geom_errorbar(aes(ymin=Lower/1000, ymax=Upper/1000), width=.1,
  #              position=position_dodge(.9)) +
  theme_bw(base_size = 9) +
  xlab(NULL) +
  ggtitle("Estimated global diversity\nof mammalian viruses") +
  ylab("Log species richness (thousands)") +
  scale_fill_manual(values = met.brewer("Juarez",6)[c(3,2,6,1)]) +
  theme(legend.position = 'n') -> g4

#####

#g4/plot_spacer()/(g1 + g3 + plot_layout(guides = 'collect')) + plot_layout(heights = c(1, 0.05, 1.5), widths = c(1,1.3,0.6))

library("cowplot")

{cairo_pdf("Figures/virus-diversity-wide.pdf", width=9.5, height=4)

ggdraw() + draw_plot(g4, 0, 0.05, 0.35, 0.9) + draw_plot(pies, 0.35, 0, 0.65, 1)

}
dev.off()

####

