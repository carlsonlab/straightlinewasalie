
library(tidyverse)
library(patchwork)
library(MetBrewer)
library(cowplot)


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
} else {
  df1$type = 'Idealized (complete)'
  df2$type = 'Realistic (incomplete)'
  df3 <- bind_rows(df1, df2)

}

}

# Save for later --------------------------------

write.table(df, "Output/false_specialists_mammalian.csv", sep=",", col.names=TRUE, row.names=FALSE)
# df <- read.csv("Output/false_specialists_mammalian.csv")

write.table(df3, "Output/false_specialists_vertebrate.csv", sep=",", col.names=TRUE, row.names=FALSE)
# df3 <- read.csv("Output/false_specialists_vertebrate.csv")

# Figures ---------------------------------------

straightLine <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")


g1 <- ggplot(df, aes(x = hosts, y = false_specialists, color = type, group = type)) +
  theme_bw() +
  geom_point(alpha = 0.05) +
  scale_color_manual(values = straightLine[c(4,10)]) +
  labs(
    title='Mammalian viruses',
    x="Hosts",
    y="Proportion of false specialists") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = c(0.25, 0.1), legend.title = element_blank())
# g1

g2 <- ggplot(df3, aes(x = hosts, y = false_specialists, color = type, group = type)) +
  theme_bw() +
  geom_point(alpha = 0.05) +
  scale_color_manual(values = straightLine[c(4,10)]) +
  ggtitle('Vertebrate parasites') + xlab("Hosts") + ylab("Proportion of false specialists") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = c(0.25, 0.1), legend.title = element_blank())

# OR...

df4 <- rbind(
  data.frame(df, set="Mammalian viruses"),
  data.frame(df3, set="Vertebrate helminth parasites")
  )


falseSpecData <- ggplot(df4, aes(x = hosts, y = false_specialists, color = type, group = type)) +
  geom_point(alpha = 0.05) +
  facet_wrap("set", scale="free_x") +
  scale_color_manual(values = straightLine[c(4,9)], labels=c("Idealized (complete) sampling", "Realistic (incomplete)")) +
  labs(x="Hosts sampled", y="Proportion of false specialists") +
  guides(colour = guide_legend(override.aes = c(alpha = 1))) +
  theme_bw(base_size = 12) +
  theme(legend.position="inside", legend.position.inside = c(0.8, 0.9), legend.title = element_blank(), legend.key.size = unit(0.4, "cm"), plot.margin = margin(0.2, 0.2, 0.4, 0.4, "cm"))

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

falseSpecSurf <- ggplot(falseSpec, aes(x=h/H, y=k, fill=PrSeenButFalseSpec)) +
  geom_tile() +
  geom_contour(aes(z=PrSeenButFalseSpec), color="white", bins=4) +
  scale_fill_gradient(low=straightLine[6], high=straightLine[10], name="Pr(seen, as false specialist)") +
  labs(x = "Fraction of hosts sampled", y="True symbiont host breadth") +
  theme_bw() +
  theme(legend.position = "top", legend.key.height = unit(0.25,"cm"), legend.key.width=unit(1,"cm"), legend.title.position = "top", plot.margin = margin(0.2, 0.2, 0.4, 0.2, "cm"))



{cairo_pdf("Figures/false-specialists.pdf", width=6.5, height=3.5)

falseSpecData

}
dev.off()


{cairo_pdf("Figures/false-specialists-extra.pdf", width=10, height=3.5)

ggdraw() +
  draw_plot(falseSpecSurf, 0, 0, 0.3, 1) +
  draw_plot(falseSpecData, 0.3, 0, 0.7, 1) +
  draw_plot_label(label=c("A", "B"), x=c(0,0.3), y=1)

}
dev.off()

