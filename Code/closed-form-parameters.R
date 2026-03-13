
library(tidyverse)
library(patchwork)
library(MetBrewer)
#library(MoMAColors)
library(vroom)
library(cowplot)

# devtools::install_github("teunbrand/ggarrow")
library(ggarrow)

# setwd("~/Documents/Active_projects/straightlinewasalie")


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
           y=19.5,
           label=expression(s[s](h)==s[g](h)),
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


