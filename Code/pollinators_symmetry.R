# Floral symmetry and specialization
# from code via CJC
# assumes local environment
# last used/modified jby, 2026.01.07

# setwd('~/Documents/Active_projects/straightlinewasalie')

library(bipartite)
library(codependent)
library(igraph)
library(reshape)
library(stringr)

library(tidyverse)

source("../shared/Rscripts/base.R")

#--------------------------------------------------------------------
# load and parse data

# floral symmetry annotations
symm <- read.csv("Data/PlantPollinator/spp_symmetry.csv", h=TRUE)

table(symm$symmetry)

spp.z <- unique(subset(symm, symmetry=="zygomorphic")$user_supplied_name)
spp.a <- unique(subset(symm, symmetry=="actinomorphic")$user_supplied_name)

intersect(spp.z, spp.a) # good to check

# list of source files from Web of Life and IWDB
files.wol <- list.files(path='Data/PlantPollinator/web-of-life_2019-02-06_013943', pattern="M_PL", full=TRUE)
files.iwdb <- list.files(path='Data/PlantPollinator/iwdb_20190422/matrices', pattern=".csv", full=TRUE)

# load the matrices as a list
webs.wol <- lapply(files.wol, function(x) read.csv(x, row.names=1, strip.white=TRUE))
names(webs.wol) <- gsub("Data/PlantPollinator/web-of-life_.+/(.+)\\.csv", "\\1", files.wol)

webs.iwdb <- lapply(files.iwdb, function(x) read.csv(x, row.names=1, strip.white=TRUE))
names(webs.iwdb) <- gsub("[ ,]", "-", gsub("Data/PlantPollinator/iwdb_.+/matrices/(.+)\\.csv", "\\1", files.iwdb))

webs <- c(webs.wol, webs.iwdb)

webs <- lapply(webs, function(w){ row.names(w) <- gsub("(.+) $", "\\1", row.names(w)); return(w)})

webs[["Robertson-1928"]] <- read.csv("Data/PlantPollinator/Robertson-1928_Carlinville.csv", h=TRUE, sep=",", row.names=1)

length(webs) # 159!


# split matrices by symmetry
webs.z <- lapply(webs, function(x) x[intersect(spp.z, rownames(x)),])
webs.a <- lapply(webs, function(x) x[intersect(spp.a, rownames(x)),])
webs.n <- lapply(webs, function(x) x[setdiff(rownames(x),c(spp.a,spp.z)),])
# issue here with trailing spaces in species names (oy) needs resolution at the source but it's turning out to be tricky ... might gain me some zyg species, tho

# what's not getting binned?
lapply(webs.n[lapply(webs.n, nrow)>0], rownames)

# deal with too-small sub-matrices
length(which(lapply(webs.z, nrow) >= 5)) # woah. Not a lot left! 42 LOL
webs.z <- webs.z[lapply(webs.z, nrow) >= 5]
webs.a <- webs.a[names(webs.z)]

#--------------------------------------------------------------------
# coextinction curve slope function

coex <- function(net, part="higher") {
ex <- second.extinct(net, participant=part, method="random", nrep=50, details=FALSE)
s <- slope.bipartite(ex, plot.it=FALSE)
return(s[[1]])
}

#--------------------------------------------------------------------


#--------------------------------------------------------------------
# estimate summary stats for sub-matrices

# one big happy table
data.df <- data.frame(matrix(0,0,18))
names(data.df) <- c('web', 'symmetry', 'gen_only', 'web.asymmetry', 'links.per.species', 'nestedness', 'n.aff', 'n.hosts', 'host.sharing', 'aff.sharing', "modularity", "z.est", 'coex.plant', 'coex.poll')

# LOOP over webs --------------------------------
for(w in names(webs.z[40:42])){

# w <- names(webs.z)[1]

# grab the relevant datasets
web.z <- as.matrix(webs.z[[w]])
web.a <- as.matrix(webs.a[[w]])

# make sure we're in 0,1 format
web.z[web.z > 1] <- 1
web.a[web.a > 1] <- 1

# make versions with only non-specialist associates
web.z.g <- web.z[,apply(web.z, 2, sum)>1]
web.a.g <- web.a[,apply(web.a, 2, sum)>1]

# bipartite matrix stats
bip.z <- bipartite::networklevel(web.z,index=c("connectance", "web asymmetry", "links per species", "nestedness", "number of species", "mean number of shared partners"))
bip.z.g <- bipartite::networklevel(web.z.g,index=c("connectance", "web asymmetry", "links per species", "nestedness", "number of species", "mean number of shared partners"))

bip.a <- bipartite::networklevel(web.a,index=c("connectance", "web asymmetry", "links per species", "nestedness", "number of species", "mean number of shared partners"))
bip.a.g <- bipartite::networklevel(web.a.g,index=c("connectance", "web asymmetry", "links per species", "nestedness", "number of species", "mean number of shared partners"))

stats <- data.frame(rbind(bip.z, bip.z.g, bip.a, bip.a.g))

colnames(stats) <- c('connectance', 'web.asymmetry', 'links.per.species', 'nestedness', 'n.aff', 'n.hosts', 'host.sharing', 'aff.sharing')

# modularity
mod.z <- modularity(cluster_walktrap(graph_from_biadjacency_matrix(web.z)))
mod.z.g <- modularity(cluster_walktrap(graph_from_biadjacency_matrix(web.z.g)))

mod.a <- modularity(cluster_walktrap(graph_from_biadjacency_matrix(web.a)))
mod.a.g <- modularity(cluster_walktrap(graph_from_biadjacency_matrix(web.a.g)))

stats$modularity <- c(mod.z, mod.z.g, mod.a, mod.a.g)

# codependence!
# first the point estimate
web.z.m <- melt(web.z) %>% filter(value>0) %>% select(1:2)
web.z.g.m <- melt(web.z.g) %>% filter(value>0) %>% select(1:2)
web.a.m <- melt(web.a) %>% filter(value>0) %>% select(1:2)
web.a.g.m <- melt(web.a.g) %>% filter(value>0) %>% select(1:2)
  
z.z <- tryCatch(copredict(web.z.m,100,100)[[2]][[2]], error = function(e) {NA})
z.z.g <- tryCatch(copredict(web.z.g.m,100,100)[[2]][[2]], error = function(e) {NA})

z.a <- tryCatch(copredict(web.a.m,100,100)[[2]][[2]], error = function(e) {NA})
z.a.g <- tryCatch(copredict(web.a.g.m,100,100)[[2]][[2]], error = function(e) {NA})

stats$z.est <- c(z.z, z.z.g, z.a, z.a.g)

# then the resampling ci (not enough large matrices for this --- add columns if resume)
# zstar.z <- tryCatch(zstar(web.z.m,10,100, plot=FALSE), error = function(e) {NA})
# zstar.a <- tryCatch(zstar(web.a.m,10,100, plot=FALSE), error = function(e) {NA})

# stats$zstar <- c(zstar.z[[1]], zstar.a[[1]])
# stats$zstarLO <- c(zstar.z[[2]][1], zstar.a[[2]][1])
# stats$zstarHI <- c(zstar.z[[2]][2], zstar.a[[2]][2])

# coextinction --- not tested in full loop; may glitch on certain matrices
coex.plant.z <- coex(web.z)
coex.plant.z.g <- coex(web.z.g)

coex.plant.a <- coex(web.a)
coex.plant.a.g <- coex(web.a.g)

coex.poll.z <- coex(web.z, part="lower")
coex.poll.z.g <- coex(web.z.g, part="lower")

coex.poll.a <- coex(web.a, part="lower")
coex.poll.a.g <- coex(web.a.g, part="lower")

stats$coex.plant <- c(coex.plant.z, coex.plant.z.g, coex.plant.a, coex.plant.a.g)
stats$coex.poll <- c(coex.poll.z, coex.poll.z.g, coex.poll.a, coex.plant.a.g)

# add stats to the data frame
data.df <- rbind(data.df, data.frame(web=w, symmetry=rep(c("Zygomorphic", "Actinomorphic"), each=2), gen_only=c(FALSE, TRUE, FALSE, TRUE), stats))

# write it all out
write.table(data.df, "Output/PlantPollinator_network_stats.csv", col.names=TRUE, row.names=FALSE, sep=",", quote=FALSE)

cat("DONE with web ", w, "\n") 

}
# END loop --------------------------------------

# read back in if it's already run and saved
data.df <- read.csv("Output/PlantPollinator_network_stats.csv")

# useful reformat
data.df$gen_only <- as.numeric(data.df$gen_only)
data.df$gen_only[data.df$gen_only==1] <- "Generalists only"
data.df$gen_only[data.df$gen_only==0] <- "All visitors"

glimpse(data.df)
table(data.df$web)

beths <- c("#C4C3C1", "#7D523A", "#F3CF26", "#2D4E7B", "#7079B2", "#54B2D8", "#2F7858", "#6D8842", "#8BAA81", "#E23639", "#DF6F85")

{cairo_pdf("Figures/PlantPollinator_plantcoex.pdf", width=3.5, height=4)
ggplot() +
	geom_line(data=data.df, aes(y=coex.plant, x=symmetry, group=web), alpha=0.75, linewidth=0.75, color="gray30") + 
	geom_point(data=data.df, aes(y=coex.plant, x=symmetry, color=symmetry), alpha=0.75, size=2) + 
	facet_wrap("gen_only", nrow=1) +
	scale_color_manual(values=beths[c(6,7)], guide="none") +
	scale_x_discrete(labels=c("Act.", "Zyg.")) +
	labs(x="Floral symmetry", y="Plant coext. robustness") +
	theme_bw(base_size=14)
}	
dev.off()

wilcox.test(coex.plant~symmetry, data=filter(data.df, gen_only=="All visitors"), alt="g") # p = 4.2e-15
wilcox.test(coex.plant~symmetry, data=filter(data.df, gen_only!="All visitors"), alt="g") # p = 0.003


{cairo_pdf("Figures/PlantPollinator_pollcoex.pdf", width=3.5, height=4)
ggplot() +
	geom_line(data=data.df, aes(y=coex.poll, x=symmetry, group=web), alpha=0.75, linewidth=0.75, color="gray30") + 
	geom_point(data=data.df, aes(y=coex.poll, x=symmetry, color=symmetry), alpha=0.75, size=2) + 
	facet_wrap("gen_only", nrow=1) +
	scale_color_manual(values=beths[c(6,7)], guide="none") +
	scale_x_discrete(labels=c("Act.", "Zyg.")) +
	labs(x="Floral symmetry", y="Visitor coext. robustness") +
	theme_bw(base_size=14)
}	
dev.off()

wilcox.test(coex.poll~symmetry, data=filter(data.df, gen_only=="All visitors"), alt="g") # p = 4.9e-15
wilcox.test(coex.poll~symmetry, data=filter(data.df, gen_only!="All visitors"), alt="g") # p = 0.0002

	


