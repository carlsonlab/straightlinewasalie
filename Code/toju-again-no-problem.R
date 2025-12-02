
library(reshape) # do NOT yell at me
library(codependent)
library(tidyverse)
library(patchwork)

setwd("~/Documents/Github/straightlinewasalie")

# 1. Robertson 1929: plant-pollinator interactions

rob1929_raw <- read.csv("./Data/Edgelists/robertson1929.csv")
rob1929 <- rob1929_raw[,c('plant','poll')]
colnames(rob1929) <- c("Plant","Pollinator")

# 2. Schleuning 2010: seed-disperser interactions

sch2010_raw <- read.csv("./Data/Edgelists/schleuning2010.csv")
sch2010_raw %>%
  select(-X) %>%
  dplyr::rename('Animal' = "Plant.species") %>%
  pivot_longer(-1, names_to = "Plant", values_to = "Interactions") %>%
  mutate(Presence = as.numeric(Interactions>0)) %>%
  select(-Interactions) %>%
  filter(Presence==1) %>%
  select(Plant, Animal) %>%
  data.frame() -> sch2010

# 3. Toju 2018: plant-arbuscular mycorrhizae

tojuraw.1 <- read.csv('./Data/Edgelists/toju descriptors.csv')
tojuraw.2 <- read.csv('./Data/Edgelists/toju.csv')
otu.list <- unique(tojuraw.1[tojuraw.1$Category.In.This.Study=='Arbuscular_Mycorrhizal',]$OTU.code)

tojuraw.2$Plant <- gsub('_TES', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_TOMA', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YKS', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_SGD', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YSD', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YAKU', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_YONA', '', tojuraw.2$Plant)
tojuraw.2$Plant <- gsub('_IRI', '', tojuraw.2$Plant)

toju <- aggregate(. ~ Plant, data=tojuraw.2, FUN=sum)
toju <- melt(toju, id=c("Plant"))
toju <- toju[toju$value>1,]
colnames(toju)[1:2] <- c("Plant","Microbe")
toju$Microbe <- gsub('X', '', toju$Microbe)
toju$Microbe <- gsub('\\.', ':', toju$Microbe)
toju <- toju[toju$Microbe %in% as.character(otu.list),]
toju <- toju[,c(1:2)]

# 4. Host-helminth relationships

helminths.raw <- read.csv('./Data/Edgelists/helminths.csv')
helminths <- helminths.raw[helminths.raw$group=='Nematoda',]
helminths <- helminths[helminths$hostgroup=='Mammalia',]
helminths <- na.omit(helminths)
helminths <- helminths[,c(2,3)]
helminths <- unique(helminths)

# We're gonna try to make these functions because that's the easiest thing to do, I guess

powerlaw <- function(df) { 
  model <- binera(df, iter = 100, plots = FALSE, subSample = NULL, koh = FALSE)
  b <- coefficients(model)['b']
  z <- coefficients(model)['z']
  df2 <- data.frame(n.host = c(1:length(unique(df[,1]))), n.par = NA)
  df2$n.par <- b*(df2$n.host^z)
  return(df2)
}

kohcolwell <- function(df) { 
  H <- length(unique(df[,1]))
  P <- length(unique(df[,2]))
  df2 <- data.frame(n.host = c(1:H, n.par = NA))
  
  degdist <- c(table(df[, 2]))
  sj <- function(num) { sum(degdist == num)}
  alpha <- function(j) {
    if (j + h <= H) { choose((H - j), h)/choose(H, h) } else { 0 }
  }
  tau.f <- function(j) { alpha(j) * sj(j) }
  
  for (h in 1:nrow(df2)) {
    tau <- sum(sapply(c(1:H), tau.f))
    df2$n.par[h] <- P - tau
  }
  
  return(df2)
}

closedform <- function(df) {
  colnames(df) <- c("Host", "Virus") # shut up. don't yell at me
  
  specs <- df %>% count(Virus) %>% filter(n==1) %>% pull(Virus) %>% unique()
  gens <- df %>% count(Virus) %>% filter(n>1) %>% pull(Virus) %>% unique()
  
  H <- df %>% pull(Host) %>% unique() %>% length()
  Hg <- df %>% filter(Virus %in% gens) %>% pull(Host) %>% unique() %>% length()
  alphag <-  df %>% filter(Virus %in% gens) %>% count(Host) %>% pull(n) %>% mean()
  etag <- df %>% filter(Virus %in% gens) %>% count(Virus) %>% pull(n) %>% mean()
  As <- length(specs)
  
  a <- alphag*Hg/H
  b <- (etag - 1)/(H)
  c <- As/H
  
  est <- function(x) {a*x/(1 + b*x) + c*x}
  df2 <- data.frame(n.host = 1:H, n.par = sapply(c(1:H), est))
  return(df2)
}

# Make some graphs

alpha <- 0.1
beths <- c("#C4C3C1", "#7D523A","#F3CF26","#2D4E7B","#7079B2","#54B2D8","#2F7858","#6D8842","#8BAA81","#E23639","#DF6F85")

samples <- curve.df(sch2010, iter = 100)
power <- powerlaw(sch2010)
koh <- kohcolwell(sch2010)
closed <- closedform(sch2010)
samples %>%
  ggplot(aes(x = n.host, y = n.par)) + 
  theme_bw() + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, color = '#C4C3C1') + 
  geom_line(data = power, aes(x = n.host, y = n.par), lwd = 1, col = '#2F7858') + 
  geom_line(data = koh, aes(x = n.host, y = n.par), lwd = 1, col = '#7079B2') + 
  geom_line(data = closed, aes(x = n.host, y = n.par), lwd = 1, col = '#54B2D8') +
  xlab("Plant hosts") + ylab("Seed dispersers") -> left1 

samples %>% left_join(power %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#2F7858') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals")+
  geom_smooth(col = 'white') -> top2 
samples %>% left_join(koh %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#7079B2') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals")+
  geom_smooth(col = 'white') -> top3
samples %>% left_join(closed %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#54B2D8') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals") +
  geom_smooth(col = 'white') -> top4

samples <- curve.df(toju, iter = 100)
power <- powerlaw(toju)
koh <- kohcolwell(toju)
closed <- closedform(toju)
samples %>%
  ggplot(aes(x = n.host, y = n.par)) + 
  theme_bw() + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, color = '#C4C3C1') + 
  geom_line(data = power, aes(x = n.host, y = n.par), lwd = 1, col = '#2F7858') + 
  geom_line(data = koh, aes(x = n.host, y = n.par), lwd = 1, col = '#7079B2') + 
  geom_line(data = closed, aes(x = n.host, y = n.par), lwd = 1, col = '#54B2D8') +
  xlab("Plant hosts") + ylab("Mycorrhizal OTUs") -> left2 

samples %>% left_join(power %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#2F7858') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals")+
  geom_smooth(col = 'white') -> middle2 
samples %>% left_join(koh %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#7079B2') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals")+
  geom_smooth(col = 'white') -> middle3
samples %>% left_join(closed %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#54B2D8') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals") +
  geom_smooth(col = 'white') -> middle4

samples <- curve.df(rob1929, iter = 100)
power <- powerlaw(rob1929)
koh <- kohcolwell(rob1929)
closed <- closedform(rob1929)
samples %>%
  ggplot(aes(x = n.host, y = n.par)) + 
  theme_bw() + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, color = '#C4C3C1') + 
  geom_line(data = power, aes(x = n.host, y = n.par), lwd = 1, col = '#2F7858') + 
  geom_line(data = koh, aes(x = n.host, y = n.par), lwd = 1, col = '#7079B2') + 
  geom_line(data = closed, aes(x = n.host, y = n.par), lwd = 1, col = '#54B2D8') +
  xlab("Plant hosts") + ylab("Floral visitors") -> left3

samples %>% left_join(power %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#2F7858') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals")+
  geom_smooth(col = 'white') -> middleagain2 
samples %>% left_join(koh %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#7079B2') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals") +
  geom_smooth(col = 'white') -> middleagain3
samples %>% left_join(closed %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#54B2D8') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Plant hosts") + ylab("Residuals") +
  geom_smooth(col = 'white') -> middleagain4

samples <- curve.df(helminths, iter = 100)
power <- powerlaw(helminths)
koh <- kohcolwell(helminths)
closed <- closedform(helminths)
samples %>%
  ggplot(aes(x = n.host, y = n.par)) + 
  theme_bw() + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, color = '#C4C3C1') + 
  geom_line(data = power, aes(x = n.host, y = n.par), lwd = 1, col = '#2F7858') + 
  geom_line(data = koh, aes(x = n.host, y = n.par), lwd = 1, col = '#7079B2') + 
  geom_line(data = closed, aes(x = n.host, y = n.par), lwd = 1, col = '#54B2D8') +
  xlab("Mammal hosts") + ylab("Nematode parasites") -> left4

samples %>% left_join(power %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#2F7858') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Mammal hosts") + ylab("Residuals") +
  geom_smooth(col = 'white') -> bottom2 
samples %>% left_join(koh %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#7079B2') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Mammal hosts") + ylab("Residuals") +
  geom_smooth(col = 'white') -> bottom3
samples %>% left_join(closed %>% dplyr::rename(n.par.est = 'n.par')) %>%
  mutate(residuals = n.par - n.par.est) %>%
  ggplot(aes(x = n.host, y = residuals)) + 
  geom_point(alpha = alpha, pch = 16, size = 1.5, col = '#54B2D8') + 
  geom_hline(yintercept = 0, linetype = 'dashed', lwd = 1, col = 'black') + 
  theme_bw() + xlab("Mammal hosts") + ylab("Residuals") +
  geom_smooth(col = 'white') -> bottom4

left1 + top2 + top4 + top3 + 
  left2 + middle2 + middle4 + middle3 + 
  left3 + middleagain2 + middleagain4 + middleagain3 + 
  left4 + bottom2 + bottom4 + bottom3 + 
  plot_layout(ncol = 4, nrow = 4)

ggsave("Figures/curve-fitting.pdf", width = 13, height = 12)
ggsave("Figures/curve-fitting.jpg", width = 13, height = 12, dpi = 600)
