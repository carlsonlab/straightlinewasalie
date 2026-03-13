
library(tidyverse)
library(MetBrewer)
library(cowplot)
library(patchwork)

helminths.raw <- read.csv('./Data/Edgelists/helminth associations cleaned v2.csv')
helminths <- helminths.raw[helminths.raw$hostgroup %in% c("Amphibia","Aves","Chondrichthyes","Osteichthyes","Reptilia","Mammalia"),]
helminths <- helminths[,c(2,3)]
helminths <- na.omit(helminths)
h <- unique(helminths)

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
degree_df <- data.frame(Proportion = rep(c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95),reps), Intercept = NA, Exponent = NA)
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
  xlab("Degree (log)") + 
  ylab("Frequency (log)") + 
  theme_bw() +
  xlim(0,4.2) + 
  ylim(0,8.5) + 
  geom_abline(intercept = exp.true[[1]], slope = exp.true[[2]], size = 1.5) -> g

for(i in 1:nrow(summary_df)) {
  g <- g + geom_abline(intercept = summary_df$Intercept[[i]], slope = summary_df$Exponent[[i]], color = cols[[summary_df$Proportion[[i]]*100]], size = 1, alpha = 0.9)
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
reps <- 500
minhosts <- 10
reps2 <- 20

# Set up data frame
proportions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
samples = 1:reps
extinction_highlevel <- expand.grid(Proportion = proportions, Iteration = samples)

for(i in 1:nrow(extinction_highlevel)) {
  prop_i <- extinction_highlevel$Proportion[i]
  
  # Generate host-based subsample 
  hosts.i <- sample(hosts, round(prop_i*n.hosts))
  n.hosts.i <- length(hosts.i)
  h.i <- h %>% filter(Host %in% hosts.i)
  n.symbionts.i <- length(unique(h.i$Parasite))
  pars.local <- h.i %>% pull(Parasite) %>% unique()
  specs.local <- table(h.i$Parasite) %>% data.frame() %>% filter(Freq == 1) %>% pull(Var1)
  specs.false <- setdiff(specs.local, specs)
  
  # Simulate coextinction
  df_extinction_i <- data.frame(E_h = rep(seq(0.05:1, by = 0.05),reps2), E_s = NA, ext_spec_local = NA, ext_spec_global = NA)
  for(j in 1:nrow(df_extinction_i)) {
    # Generate another host-based subsample 
    hosts.j <- sample(hosts.i, round((1-df_extinction_i$E_h[j])*n.hosts.i))
    h.j <- h.i %>% filter(Host %in% hosts.j)
    symbionts.j <- length(unique(h.j$Parasite))
    df_extinction_i$E_s[j] <- 1-(symbionts.j/n.symbionts.i)
    
    pars.j <- unique(h.j$Parasite)
    pars.extinct <- setdiff(pars.local, pars.j)
    df_extinction_i$ext_spec_local[j] <- sum(pars.extinct %in% specs.local)/length(pars.extinct)
    df_extinction_i$ext_spec_global[j] <- sum(pars.extinct %in% specs)/length(pars.extinct)
    df_extinction_i$ext_spec_false[j] <- sum(pars.extinct %in% specs.false)/length(pars.extinct)
  }
  df_summary_i <- df_extinction_i %>%
    group_by(E_h) %>% summarize(E_s = mean(E_s),
                                ext_spec_local = mean(ext_spec_local),
                                ext_spec_global = mean(ext_spec_global),
                                ext_spec_false = mean(ext_spec_false)) %>%
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
  geom_line(size = 1, alpha = 0.9) + 
  theme_bw() + 
  xlab("Host extinction rate") + ylab("Symbiont extinction rate")+
  scale_color_gradientn(colors = cols, name = "Proportion\nof hosts \nsampled\n") -> panel2

panel1 + panel2

extinction_df %>%
  group_by(Proportion, E_h) %>%
  summarize(ext_spec_false = mean(ext_spec_false, na.rm = TRUE)) %>%
  ggplot(aes(x = E_h, y = ext_spec_false, color = Proportion, group = Proportion)) + 
  geom_line(size = 1, alpha = 0.9) + 
  theme_bw() + 
  ylim(0,1) + 
  xlab("Host extinction rate") + ylab("Proportion of coextinctions by false specialists") +
  scale_color_gradientn(colors = cols, name = "Proportion\nof hosts \nsampled\n") -> panel3

extinction_df %>%
  group_by(Proportion, E_h) %>%
  summarize(ext_spec_global = mean(ext_spec_global, na.rm = TRUE)) %>%
  ggplot(aes(x = E_h, y = ext_spec_global, color = Proportion, group = Proportion)) + 
  geom_line(size = 1, alpha = 0.9) + 
  theme_bw() + 
  ylim(0,1) + 
  xlab("Host extinction rate") + ylab("Specialist proportion of coextinctions") +
  scale_color_gradientn(colors = cols, name = "Proportion\nof hosts \nsampled\n") -> panel4

# (panel1 + panel2) / (panel3 + panel4) + plot_layout(guides = 'collect')

panel1 + panel2 + panel4 + plot_layout(guides = 'collect') + plot_annotation(tag_levels = c("A", "B", "C")) & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave("./Figures/coextinction.pdf", width = 11.5, height = 3.5)
