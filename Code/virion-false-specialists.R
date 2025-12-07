library(tidyverse)
library(patchwork)
library(MetBrewer)
library(vroom)
library(cowplot)
library(MoMAColors)

# setwd("~/Github/straightlinewasalie/Data")

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
clo %>%
  select(Host, Virus) %>% distinct() %>%
  group_by(Virus) %>% summarize(HostRange = n()) %>%
  filter(HostRange==1) %>% pull(Virus) -> true_specialists

##### Charting it through real-time

params <- data.frame(
  year = c(1960:2024),
  hosts = NA, viruses = NA,
  prop_false = NA
)

for(year in 1960:2024) {

  clo %>% filter(Year <= year) -> clo_year

  clo_year %>% pull(Host) %>% unique() %>% length() -> hosts
  clo_year %>% pull(Virus) %>% unique() %>% length() -> viruses

  clo_year %>%
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>%
    filter(HostRange==1) %>% pull(Virus) -> obs_specialists

  prop_false <- 1-(sum(as.numeric(obs_specialists %in% true_specialists))/length(obs_specialists))

  params[params$year == year,] <- c(year, hosts, viruses, prop_false)
  print(year)
}

##### Charting it through sampling

params2 <- data.frame(
  hosts = NA, viruses = NA,
  prop_false = NA
)

hosts <- unique(clo$Host)

for(i in 1:1000) {

  n_hosts_i <- sample(c(5:total_hosts),1)
  hosts_i <- sample(hosts, n_hosts_i)

  clo %>% filter(Host %in% hosts_i) -> clo_i
  n_viruses_i <- length(unique(clo_i$Virus))
  clo_i %>%
    select(Host, Virus) %>% distinct() %>%
    group_by(Virus) %>% summarize(HostRange = n()) %>%
    filter(HostRange==1) %>% pull(Virus) -> obs_specialists
  prop_false <- 1-(sum(as.numeric(obs_specialists %in% true_specialists))/length(obs_specialists))

  params2[i,] <- c(n_hosts_i, n_viruses_i, prop_false)
  print(i)
}

### Stash the results for later
write.table(params, "Output/time_sampling_false_specialists.csv", col.names=TRUE, row.names=FALSE, sep=",")
# params <- read.csv("Output/time_sampling_false_specialists.csv")

write.table(params2, "Output/sim_sampling_false_specialists.csv", col.names=TRUE, row.names=FALSE, sep=",")
# params2 <- read.csv("Output/sim_sampling_false_specialists.csv")


##### Make plots
beths <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")


# AND ...

# P(\hat{k}=1|H,h,k) = {h\binom{H-h}{k-1}}/{\binom{H}{k}}
# P(\hat{k}=0|H,h,k) = {\binom{H-h}{k}}/{\binom{H}{k}}
# check these out against in-text calcs
H <- 6e3
h <- 600
k <- 10
h*choose(H-h, k-1)/choose(H, k)
choose(H-h, k)/choose(H, k)

h*choose(H-h, k-1)/choose(H, k) / (1 - choose(H-h, k)/choose(H, k))

falseSpec <- expand.grid(h=seq(1,5001,by=100),k=2:26)

H <- 1e4

falseSpec$PrFalseSpec <- mapply(
  function(h,k) h*choose(H-h, k-1)/choose(H, k),
  falseSpec$h,
  falseSpec$k
)
falseSpec$PrNotSeen <- mapply(
  function(h,k) choose(H-h, k)/choose(H, k),
  falseSpec$h,
  falseSpec$k
)
falseSpec$PrSeen <- 1-falseSpec$PrNotSeen
falseSpec$PrSeenButFalseSpec <- falseSpec$PrFalseSpec / falseSpec$PrSeen

glimpse(falseSpec)
hist(falseSpec$PrFalseSpec)
hist(falseSpec$PrSeenButFalseSpec)
hist(falseSpec$PrNotSeen)

falseSpecSurf <- ggplot(falseSpec) +
  geom_tile(aes(x=h/H, y=k, fill=PrSeenButFalseSpec)) +
  geom_contour(aes(x=h/H, y=k, z=PrSeenButFalseSpec), color="white", breaks=c(0.1,0.5,0.9)) +
  annotate("label", label="0.1", x=0.47, y=8, border.color=NA, text.color="white", fill=NA, size=3) +
  annotate("label", label="0.5", x=0.47, y=4, border.color=NA, text.color="white", fill=NA, size=3) +
  annotate("label", label="0.9", x=0.22, y=2.5, border.color=NA, text.color="white", fill=NA, size=3) +
  scale_fill_gradient(low=beths[6], high=beths[10], name="Pr(seen, as false specialist)") +
  scale_y_continuous(breaks=c(2,10,20)) +
  labs(x = "Fraction of hosts sampled", y="True symbiont host breadth") +
  theme_bw() +
  theme(legend.position = "none", legend.key.height = unit(0.25,"cm"), legend.key.width=unit(1,"cm"), legend.title.position = "top", plot.margin = margin(0.25, 0.2, 0.4, 0.4, "cm"))

# and then from data, above:
falseSpecSurf # test viewing

falseSpecData <- ggplot() +
  geom_point(data=params2, aes(x = hosts, y = prop_false), color=beths[5], alpha=0.5) +
  geom_point(data=params, aes(x = hosts, y = prop_false, color = year)) +
  scale_color_gradient(low=beths[3], high=beths[7], name="Year of detection") +
  labs(x = "Hosts sampled", y = "Proportion of specialists later\nrevealed to be generalists") +
  annotate("text", label="Simulated\nsubsampling", x=250, y=0.07, color=beths[5], lineheight=1) +
  annotate("text", label="Empirical host-\nvirus networks", x=1150, y=0.35, color=beths[7], lineheight=1) +
  theme_bw() +
  theme(legend.position="inside",
        legend.position.inside = c(0.725,0.85),
        legend.direction="horizontal",
        legend.key.height = unit(0.25, "cm"),
        legend.title.position = "top",
        plot.margin = margin(0.25, 0.2, 0.4, 0.4, "cm"))

falseSpecData # preview


{cairo_pdf("Figures/false-specialists_analytic_data.pdf", width=6.5, height=3)

ggdraw() + draw_plot(falseSpecSurf, 0, 0, 0.45, 1) + draw_plot(falseSpecData , 0.45, 0, 0.55, 1) + draw_plot_label(label=c("A", "B"), x=c(0,0.45), y=1)

}
dev.off()


