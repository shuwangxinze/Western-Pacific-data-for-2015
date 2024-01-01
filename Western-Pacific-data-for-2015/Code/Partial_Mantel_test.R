# Environmental Data
env <- read.delim(file.choose(),sep = "\t", row.names = 1, check.names = FALSE)
env <- env[c("Salinity","Temperature","DO","PO4","SiO3","NO2","NO3","NH4")]# or select the environment variable from the forward selection
ENSO <- env[c("SSTA")]


spe <- read.delim(file.choose(), sep = "\t", check.names = FALSE)
spe <- data.frame(t(spe))
spe <- decostand(spe,method = 'hellinger')
# Bray-curtis diversity of species composition between communities
library(vegan)

comm_dis <- vegdist(spe, method = "bray")

library(geosphere)
site_dis <- distm(site_env[c("Lon", "Lat")]) / 1000  
site1_dis <- log(site_dis+1, exp(1))
rownames(site_dis) <- rownames(site_env)
colnames(site_dis) <- rownames(site_env)
site_dis <- as.dist(site_dis)

#For each environment variable, the Euclidean distance measure is calculated
ENV_dis <- vegdist(env, method = "euclidean")  
ENSO_dis <- vegdist(site_env["SSTA"], method = "euclidean") 

#MuMIn package stdize () can achieve the function of the distance measure of geography or environment standardization, reference the: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4961435/
library(MuMIn)
ENV_dis <- stdize(ENV_dis)

site_dis <- stdize(site_dis)

ENSO_dis <- stdize(ENSO_dis)

library(ecodist)
#Correlation with environmental factors, control for ENSO and geographical location
fit_MRM5 <- MRM(comm_dis~ENV_dis + ENSO_dis+ site_dis,
                nperm = 999, method = "linear")

#Correlation with geographical distance, control for ENSO and environmental variables
fit_MRM5 <- MRM(comm_dis~site_dis + ENSO_dis+ ENV_dis,
                nperm = 999, method = "linear")


#Correlation with ENSO, controlling for geographical distance and environmental variables
fit_MRM5 <- MRM(comm_dis~ENSO_dis + site_dis+ ENV_dis,
                nperm = 999, method = "linear")