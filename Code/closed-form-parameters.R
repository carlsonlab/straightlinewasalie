
library(tidyverse)
library(patchwork)
library(MetBrewer)
library(MoMAColors)
library(vroom)
library(cowplot)

# devtools::install_github("teunbrand/ggarrow")
library(ggarrow)

# setwd("~/Documents/Github/straightlinewasalie/Data")

# Access virus data
# library(virionData)
# get_versioned_data(version = "17636397", dir_path = "./Data/Temp")

# Set up virus data
virion <- vroom("./Data/Temp/17636397/virion.csv.gz")
virion %>%
  mutate(Host = str_to_lower(Host), Virus = str_to_lower(Virus)) %>%
  filter(DetectionMethod %in% c("PCR/Sequencing","Isolation/Observation")) %>%
  filter(HostClass=='mammalia') %>%
  filter(HostOrder=="chiroptera", VirusFamily=="coronaviridae") %>%
  filter(!(Host=='homo sapiens')) %>%
  filter(!is.na(Host)) %>% filter(!is.na(Virus)) %>%
  select(Host, Virus, VirusGenus, VirusFamily, ICTVRatified) %>%
  distinct() -> h

hosts <- unique(h$Host)
n.hosts <- length(hosts)
pars <- unique(h$Virus)
n.pars <- length(pars)
specs <- table(h$Virus) %>% data.frame() %>% filter(Freq == 1) %>% pull(Var1)
n.specs <- length(specs)
gens <- table(h$Virus) %>% data.frame() %>% filter(Freq > 1) %>% pull(Var1)
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
  h.spec.i <- h.i %>% filter(Virus %in% specs)
  h.gen.i <- h.i %>% filter(Virus %in% gens)

  exp.j$h[i] <- length(unique(h.i$Host))
  exp.j$p[i] <- length(unique(h.i$Virus))
  exp.j$ps[i] <- length(unique(h.spec.i$Virus))
  exp.j$pg[i] <- length(unique(h.gen.i$Virus))
  exp.j$phs[i] <- exp.j$ps[i]/exp.j$h[i]
  exp.j$phg[i] <- mean(table(h.gen.i$Host))
  exp.j$hpg[i] <- mean(table(h.gen.i$Virus))

  print(i)

}

if(j == 1) {
  exp <- exp.j
} else {
  exp <- bind_rows(exp, exp.j)
}
}

### MECHANISTIC MODEL

H <- n.hosts
Hg <- h %>% filter(Virus %in% gens) %>% pull(Host) %>% unique() %>% length()
alphag <-  h %>% filter(Virus %in% gens) %>% count(Host) %>% pull(n) %>% mean()
etag <- h %>% filter(Virus %in% gens) %>% count(Virus) %>% pull(n) %>% mean()
As <- n.specs

A <- alphag*Hg/H
B <- (etag - 1)/(H)
C <- As/H

est.s <- function(x) {C*x}
est.g <- function(x) {A*x/(1 + B*x)}
est.hpg <- function(x) {1+B*x}
est <- function(x) {A*x/(1 + B*x) + C*x}
est <- data.frame(h = 1:H,
                  p = sapply(c(1:H), est),
                  s = sapply(c(1:H), est.s),
                  g = sapply(c(1:H), est.g),
                  phs = C,
                  phg = alphag,
                  hpg = sapply(c(1:H), est.hpg))

### GRAPHS
straightLine <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")

cols <- met.brewer(name="Benedictus", n=11, type="discrete")
pink <- cols[[5]]
blue <- cols[[8]]

exp %>%
  ggplot(aes(x = h, y = p)) +
  geom_point(size = 3, alpha = 0.1, pch = 16, color = 'grey70') +
  geom_line(aes(x = h, y = p), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Total viruses") -> g1; g1
exp %>%
  ggplot(aes(x = h, y = ps)) +
  geom_point(size = 3, alpha = 0.1, pch = 16, color = pink) +
  geom_line(aes(x = h, y = s), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Specialist viruses")  -> g2; g2
exp %>%
  ggplot(aes(x = h, y = phs)) +
  geom_point(size = 3, alpha = 0.1, pch = 16, color = pink) +
  geom_line(aes(x = h, y = phs), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Specialist viruses per host")  -> g3; g3
exp %>%
  ggplot(aes(x = h, y = pg)) +
  geom_point(size = 1, alpha = 0.1, pch = 16, color = blue) +
  geom_line(aes(x = h, y = g), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Generalist viruses")  -> g4; g4
exp %>%
  ggplot(aes(x = h, y = hpg)) +
  geom_point(size = 3, alpha = 0.1, pch = 16, color = blue) +
  geom_line(aes(x = h, y = hpg), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Hosts per generalist virus")  -> g5; g5
exp %>%
  ggplot(aes(x = h, y = phg)) +
  geom_point(size = 3, alpha = 0.1, pch = 16, color = blue) +
  geom_line(aes(x = h, y = phg), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw()+ xlab("Hosts") + ylab("Generalist viruses per host") -> g6; g6


{cairo_pdf("Figures/closed-form.pdf", width=9.5, height=6)

(g1+g2+g3)/(g4+g5+g6) + plot_annotation(tag_levels = "A")

}
dev.off()

exp %>%
  ggplot(aes(x = h, y = p)) +
  geom_point(size = 3, alpha = 0.1, pch = 16, color = 'grey70') +
  geom_line(aes(x = h, y = p), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Total viruses") -> g1; g1
exp %>%
  ggplot(aes(x = h, y = ps)) +
  geom_point(size = 3, alpha = 0.1, pch = 16, color = pink) +
  geom_line(aes(x = h, y = s), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Specialist viruses")  -> g2; g2

exp %>%
  ggplot(aes(x = h, y = pg)) +
  geom_point(size = 1, alpha = 0.1, pch = 16, color = blue) +
  geom_line(aes(x = h, y = g), col = 'black', data = est, size = 0.8, lty = 5) +
  theme_bw() + xlab("Hosts") + ylab("Generalist viruses")  -> g4; g4

exp_plot <- exp %>% pivot_longer(c(p, ps, pg), values_to="symbionts", names_to="subset")

est_plot <- est %>% pivot_longer(c(p, s, g), values_to="symbionts", names_to="subset") %>% mutate(subset=factor(subset, c("p", "g", "s")))

# one panel conceptual figure
straightLine <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")

{cairo_pdf("Figures/symbiont-subsets.pdf", width=5, height=4.5)

ggplot(
  est_plot, aes(x=h, y=symbionts, group=subset, color=subset)
  ) +
  geom_line(size=1.2) +
  geom_arrow_curve(
    x=100,
    y=20,
    xend=76,
    yend=37,
    color="black",
    curvature=-0.55,
    arrow_head = arrow_head_halfwing(),
    linewidth=1
    ) +
  # ss = sg
  annotate("text",
           x=103,
           y=18,
           label=expression(S[s]==S[g]),
           hjust=0) +
  scale_color_manual(
    values=straightLine[c(4, 6, 7)],
    name="Symbiont set",
    labels=c("Total", "Generalists", "Specialists")) +
  labs(x="Hosts", y="Symbionts") +
  theme_bw(base_size=14) +
  theme(
    legend.position="inside",
    legend.position.inside=c(0.25,0.8)
  )

}
dev.off()


