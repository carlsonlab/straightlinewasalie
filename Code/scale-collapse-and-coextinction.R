
library(tidyverse)
library(MetBrewer)
library(cowplot)
library(patchwork)

h <- read_csv("./straightlinewasalie/Data/helminths.csv")

# Total values

hosts <- unique(h$Host)
n.hosts <- length(hosts)
pars <- unique(h$Parasite)
n.pars <- length(pars)
specs <- table(h$Parasite) %>% data.frame() %>% filter(Freq == 1) %>% pull(Var1)
n.specs <- length(specs)
gens <- table(h$Parasite) %>% data.frame() %>% filter(Freq > 1) %>% pull(Var1)
n.gens <- length(gens)

# Fit full degree distribution
h %>% 
  group_by(Parasite) %>%
  summarize(Degree = n()) %>%
  group_by(Degree) %>%
  summarize(Frequency = n()) -> dd

#Prove to yourself a power law is real
#dd %>% ggplot(aes(x = log(Degree), y = log(Frequency))) + geom_point() + geom_smooth(method = "lm")

exp.true <- lm(log(dd$Frequency) ~ log(dd$Degree))$coefficients

### EXPERIMENTAL DATA ON........  DEGREE DISTRIBUTIONS!

# Set up repetitions
reps <- 100
minhosts <- 10

# Set up data frame
degree_df <- data.frame(Proportion = rep(seq(0.01:0.99,by=0.01),reps), Intercept = NA, Exponent = NA)
degree_df <- degree_df %>%
  mutate(hosts_i = round(Proportion*n.hosts))

for(i in 1:nrow(degree_df)) {
  
  # Generate host-based subsample 
  hosts.i <- sample(hosts, degree_df$hosts_i[i])
  n.hosts.i <- length(hosts.i)
  h.i <- h %>% filter(Host %in% hosts.i)
  
  # Estimate degree distribution
  h.i %>% 
    group_by(Parasite) %>%
    summarize(Degree = n()) %>%
    group_by(Degree) %>%
    summarize(Frequency = n()) -> dd.i
  exp.i <- lm(log(dd.i$Frequency) ~ log(dd.i$Degree))$coefficients
  degree_df$Intercept[i] <- exp.i[[1]]
  degree_df$Exponent[i] <- exp.i[[2]]
  
  print(i)
}

degree_df %>%
  group_by(Proportion) %>%
  summarize(Intercept = mean(Intercept), Exponent = mean(Exponent)) -> summary_df

cols <- met.brewer("Hokusai1", n=100)

summary_df %>% 
  ggplot() + 
  xlab("\n log(Degree) \n") + 
  ylab("\n log(Frequency) \n") + 
  theme_bw() +
  xlim(0,4.5) + 
  ylim(0,9) + 
  geom_abline(intercept = exp.true[[1]], slope = exp.true[[2]], size = 1.5) -> g

for(i in 1:nrow(summary_df)) {
  g <- g + geom_abline(intercept = summary_df$Intercept[[i]], slope = summary_df$Exponent[[i]], color = cols[[summary_df$Proportion[[i]]*100]], size = 1, alpha = 0.5)
}
g -> panel1

# 
# # Silly little thing that makes a legend 
# summary_df %>% ggplot(aes(x = Proportion, y = Proportion, color = Proportion)) + 
#   geom_point() +
#   scale_color_gradientn(colors = cols) -> g2
# leg <- cowplot::get_legend(g2)
# 
# g + inset_element(leg, 0.8, 0.6, 0.9, 0.9) -> panel1


### EXPERIMENTAL DATA ON........  COEXTINCTIONS!!

# Set up repetitions
reps <- 100
minhosts <- 10
reps2 <- 5

# Set up data frame
proportions <- seq(0.01:0.99,by=0.01)
samples = 1:reps
extinction_highlevel <- expand.grid(Proportion = proportions, Iteration = samples)

for(i in 1:nrow(extinction_highlevel)) {
  prop_i <- extinction_highlevel$Proportion[i]
  
  # Generate host-based subsample 
  hosts.i <- sample(hosts, round(prop_i*n.hosts))
  n.hosts.i <- length(hosts.i)
  h.i <- h %>% filter(Host %in% hosts.i)
  n.symbionts.i <- length(unique(h.i$Parasite))
  
  # Simulate coextinction
  df_extinction_i <- data.frame(E_h = rep(seq(0.05:1, by = 0.05),reps2), E_s = NA)
  for(j in 1:nrow(df_extinction_i)) {
    # Generate another host-based subsample 
    hosts.j <- sample(hosts.i, round((1-df_extinction_i$E_h[j])*n.hosts.i))
    h.j <- h.i %>% filter(Host %in% hosts.j)
    symbionts.j <- length(unique(h.j$Parasite))
    df_extinction_i$E_s[j] <- 1-(symbionts.j/n.symbionts.i)
  }
  df_summary_i <- df_extinction_i %>%
    group_by(E_h) %>% summarize(E_s = mean(E_s)) %>%
    mutate(Proportion = extinction_highlevel$Proportion[i],
           Iteration = extinction_highlevel$Iteration[i])
  if (i == 1) { 
  extinction_df <- df_summary_i
    } else {
  extinction_df <- bind_rows(extinction_df, df_summary_i)
    }
  print(i)
}

extinction_df %>%
  group_by(Proportion, E_h) %>%
  summarize(E_s = mean(E_s)) %>%
  ggplot(aes(x = E_h, y = E_s, color = Proportion, group = Proportion)) + 
  geom_line(size = 1, alpha = 0.5) + 
  theme_bw() + 
  xlab("\n Host extinction (prop.) \n") + ylab("\n Symbiont extinction (prop.) \n")+
  scale_color_gradientn(colors = cols) -> panel2

panel1 + panel2

