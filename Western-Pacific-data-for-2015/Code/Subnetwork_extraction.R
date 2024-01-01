library(preprocessCore)
library(impute)
library(WGCNA)
otu <- read.delim(file.choose(),row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)#每个样本中的代表节点
otu1 <- t(otu)
#Network edge file triplet table
otu_net<- read.delim(file.choose(), row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
detach("package:Hmisc", unload = TRUE)
otu_net$from <-as.factor(otu_net$from)
otu_net$to <-as.factor(otu_net$to)
c(as.character(otu_net$from), as.character(otu_net$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> otu_nodes
colnames(otu_nodes) <- c("name", "degree")
library(Hmisc)
library(igraph)
otu_graph <- graph_from_data_frame(otu_net, vertices = otu_nodes, directed = FALSE )
#Subnetworks are extracted according to the representative otu of each sample
otu_subNet <- extract_subNet(graph = otu_graph, data = otu1, node = otu_nodes)

otu_subNet$sub_networkParameters#Topology parameters for each subnetwork
otu_subNet$node#
Each subnetwork contains nodes and degrees
otu_subNet$edge#
The edge that each subnetwork
#Output topology properties of each subnetwork, Point and edge files
write.csv(otu_subNet1$sub_networkParameters,"Topological properties of subnetworks.csv")
write.csv(otu_subNet$node$S0701,"Subnetwork_noed_01.csv")
#And so on until the last subnetwork
...
write.csv(otu_subNet$node$S0732,"Subnetwork_noed_32.csv")


write.csv(otu_subNet$edge$S0701,"Subnetwork_edge_01.csv")
#And so on until the last subnetwork
...
write.csv(otu_subNet$edge$S0732,"Subnetwork_edge_32.csv")



#The subnetwork matrix is extracted from the global matrix and the global efficiency of the subnetwork is calculated
#global matrix 
neetwork_adj <- read.delim(file.choose(), row.names = 1, sep = '\t', check.names = FALSE)
#The row and column names that need to be extracted. Based on the node of each subnetwork
sub <- read.delim(file.choose(), sep = '\t', check.names = FALSE)
A <- as.character(sub$otu_ID)
#And so on until the last subnetwork
...
i <- as.character(sub$otu_ID)
neetwork_adj01 <- neetwork_adj[A,A]
#And so on until the last subnetwork
...
neetwork_adj32 <- neetwork_adj[i,i]
#The adjacency matrix, the adjacency list, gets the weighted undirected network
g_01 <- graph_from_adjacency_matrix(as.matrix(neetwork_adj32), mode = 'undirected', weighted = TRUE, diag = FALSE)
#And so on until the last subnetwork
...
g_32 <- graph_from_adjacency_matrix(as.matrix(neetwork_adj32), mode = 'undirected', weighted = TRUE, diag = FALSE)
#Calculate global efficiency and record
global_efficiency(g)

#The correlation of subnetwork topological properties with environment variables
topo_cor <- read.delim(file.choose(),row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#Calculate the correlation matrix, using Pearson's correlation coefficient here
cor_pearson <- cor(otu, method = 'pearson')
#Visualization
library(corrplot)
col <- colorRampPalette(c('#3C5488FF','white','#E64B3599'))
#Use cor.mtest for significance testing
res1 <- cor.mtest(cor_pearson, conf.level = .95)
#Extract the p-value matrix
p.mat = res1$p
corrplot(cor_pearson, method = 'circle', col = col(100),p.mat = res1$p, sig.level = c(.001, .01, .05),insig = "label_sig",
         pch.cex = 0.5,type = "upper",diag = F)