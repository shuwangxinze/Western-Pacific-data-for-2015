library(geosphere)
#the geographic location data of the sampling sites
site <- read.delim("D:/1225/WP_data/Western-Pacific-data-for-2015/data/Environment_variable.txt", sep = '\t', row.names = 1, check.names = FALSE)
#Calculate the geographic distance between sampling points
#Calculate geographical distance based on latitude and longitude (default distance unit, meters)
#distm() requires two columns of data, the first for longitude and the second for latitude
site_dis <- geosphere::distm(site[c('Lon', 'Lat')])
rownames(site_dis) <- rownames(site)
colnames(site_dis) <- rownames(site)
#The geographic distance matrix of sampling points is converted into a data frame structure with pairwise corresponding values
site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)
library(vegan)
#vegdist() calculated the Bray-curtis heterogeneity matrix of species composition between communities
#And the Bray-curtis similarity is obtained by 1-Bray-curtis heterogeneity
library(vegan)
library(permute)
library(lattice)
otu <- read.delim("D:/1225/WP_data/Western-Pacific-data-for-2015/data/16S.txt",row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE) 
spe <- data.frame(t(otu))
comm_dis <- vegdist(spe, method = "bray")
comm_sim <- log(1 - comm_dis, exp(1))
#The matrix is transformed into a data frame structure corresponding to the values of pairwise communities
diag(comm_sim) <-  0  #Remove the diagonal values in the community similarity matrix, which are the self-similarity of the sample
comm_sim[upper.tri(comm_sim)] <- 0  #The community similarity matrix is symmetric, so only the values of the semi-triangular (below triangle) region can be selected
comm_sim <- reshape2::melt(comm_sim)
comm_sim <- subset(comm_sim, value != 0)

##Sampling site distance and community similarity data were combined
comm_dis <- merge(comm_sim, site_dis, by = c('Var1', 'Var2'))
names(comm_dis) <- c('site1', 'site2', 'comm_sim', 'site_dis')
comm_dis$site_dis_km <- comm_dis$site_dis/1000
DDR <- lm(comm_sim~site_dis_km, data = comm_dis)
summary(DDR)
#Draw a DDR curve
library(ggplot2)
p <- ggplot(data=comm_dis, aes(x=site_dis_km, y=comm_sim)) +
  geom_point() + 
  scale_color_manual(values = c('#e57373'),) +
  geom_smooth(data = comm_dis, aes(x = site_dis_km, y = comm_sim), method = 'lm', se = FALSE) +
  theme_bw()+
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
  # legend.position = 'none', plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limit = c(0, 1)) +
  labs(x = 'Distance (km)', y = 'Bray-curtis similarity', title = 'EM')

