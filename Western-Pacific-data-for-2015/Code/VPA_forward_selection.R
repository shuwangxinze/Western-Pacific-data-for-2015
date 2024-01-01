otu <- read.delim(file.choose(), row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
otu_hel <- decostand(otu, method = 'hellinger')

env <- read.delim(file.choose(), row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
library(adespatial)
library(ape)
library(spdep)
library(vegan)
library(ade4)
install.packages("packfor", repos="http://R-Forge.R-project.org")
library(packfor)

install.packages("AEM", repos="http://R-Forge.R-project.org")
library(AEM)
install.packages("PCNM", repos="http://R-Forge.R-project.org")
library(PCNM)

# Calculate the Bray-curtis distance
dis_bray <- vegdist(otu_hel, method = 'bray')

# Forward selection of environment variables
env <- env[c("Salinity","Temperature","SSTA","DO","PO4","SiO3","NO2","NO3","NH4")]
mite.env.rda <-rda(otu_hel, env)

mite.env.R2a <-RsquareAdj(mite.env.rda)$adj.r.squared

mite.env.fwd <-forward.sel(otu_hel, env, adjR2thresh=mite.env.R2a,nperm=99)

env.sign <-sort(mite.env.fwd$order)

env.red <-env[,c(env.sign)]

colnames(env.red)


# variance decomposition based on MRM # assessing how a given geography, SSTA, and environmental factors affect a community 
# geographic distance between sampling points based on latitude and longitude using the geosphere package function distm () 
# default distance in meters, divide by 1000 to convert to kilometers and back again by LN (x + 1)
library(geosphere)

geo_dis <- distm(env[c('Lon', 'Lat')]) / 1000
geo_dis <- log(geo_dis+1, exp(1))
rownames(geo_dis) <- rownames(env)
colnames(geo_dis) <- rownames(env)
geo_dis <- as.dist(geo_dis)

# Significant environmental factors derived from forward selection
env1 <- env[c("SiO3","Temp","PO4","NO3")]

# For environment variables, calculate Euclidean distance
env_dis <- vegdist(env1, method = 'euclidean')
ENSO_dis <- vegdist(env['SSTA'], method = 'euclidean')


The variance decomposition of MRM was used to analyze the effects of geographic, climatic and environmental factors on community composition
# MRM variance decomposition code varpart4 () from: https://github.com/csdambros/BioGeoAmazonia
source('RFunctions.R')

r2part <- varpart4(
  spe_dis,
  list(geo_dis),
  list(env_dis),
  list(ENSO_dis)
)

r2part
summary(r2part)
plot(r2part, bg = c("#3C5488FF","#de572d","#FFBE7A"), Xnames = c('Geo','Env','ENSO'), cutoff = -1, digits = 2)
