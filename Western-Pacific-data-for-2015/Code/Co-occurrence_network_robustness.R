#Robustness after removing keystone nodes
#Direct import of existing adjacency matrices (weighted)
cormatrix<-read.delim(file.choose(),row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#Run the following after importing the existing adjacency matrix
cormatrix[is.na(cormatrix)]<-0
diag(cormatrix)<-0    #no links for self-self    
sum(abs(cormatrix)>0)/2  #this should be the number of links.
sum(colSums(abs(cormatrix))>0)
network.raw<-cormatrix[colSums(abs(cormatrix))>0,colSums(abs(cormatrix))>0]
sp.ra2 <-sp.ra[colSums(abs(cormatrix))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched
#get the keystone species list
node.attri<-read.delim(file.choose(),row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#module.hub<-as.character(node.attri$Name[node.attri$Zi > 2.5 & node.attri$Pi <= 0.62])
module.hub<-as.character(node.attri$Name)
#consider cascade effects: removed species will further influence the remaining nodes



#Reference from: https://github.com/Mengting-Maggie-Yuan/warming-network-complexity-stability
rand.remov2.once<-function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=T){
  rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
  id.rm<-sample(keystonelist, rm.num2)
  net.Raw=netRaw #don't want change netRaw
  
  net.new=net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
  if (nrow(net.new)<2){
    0
  } else {
    sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]
    
    if (abundance.weighted){
      net.stength= net.new*sp.ra.new
    } else {
      net.stength= net.new
    }
    
    sp.meanInteration<-colMeans(net.stength)
    
    
    while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
      id.remain<- which(sp.meanInteration>0) 
      net.new=net.new[id.remain,id.remain]
      sp.ra.new=sp.ra.new[id.remain]
      
      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }
      
      if (length(net.stength)>1){
        sp.meanInteration<-colMeans(net.stength)
      } else{
        sp.meanInteration<-0
      }
      
    }
    
    remain.percent<-length(sp.ra.new)/length(sp.ra)
    
    remain.percent}
}

rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}
#Each group is calculated separately
Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=cormatrix, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

write.csv(Unweighted.simu,"Unweighted_Robustness.csv")

write.csv(Weighted.simu,"Weighted_Robustness.csv")

# Visualization
library(ggplot2)
library(ggsignif)
library(ggpubr)

# Take the value when 50% nodes are removed to make a bar chart
plot_data1 <- read.delim(file.choose(),row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
# Remove the values of all nodes to make a scatter plot
plot_data2 <- read.delim(file.choose(), sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
p<- ggplot()+ 
  geom_bar(data=plot_data1,mapping=aes(x=group,y=mean), 
           fill = "white",
           size = 1.5,
           color = c("#3d4b96","#76c180","#e7a721","#e4826f"),
           position="dodge", 
           stat="identity", 
           width = 0.6)+  
  geom_jitter(data=plot_data2, 
             mapping=aes(x=group,y=mean,fill = group,colour = group),
             size = 2,
              height = 0.05,
              
              width = 0.1)+ 
  scale_color_manual(values = c("#3d4b96","#76c180","#e7a721","#e4826f"))+ 
  
 
  geom_errorbar(data=plot_data1,mapping=aes(x = group,ymin = mean-sd, ymax = mean+sd), 
                width = 0.3,
                color = c("#3d4b96","#76c180","#e7a721","#e4826f"), 
                size=0.8)+ 
  scale_y_continuous(limits =c(0, 1) ,expand = c(0,0))+ 
  theme_classic(  
    base_line_size = 1
  )+
  labs(title="Group",x="",y="Robustness")+ 
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", 
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  
                                   color = "black", 
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0), 
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  ) +
  geom_signif(data=plot_data2,
              comparisons = list(c("A", "B"), 
                                 c("A", "C"),
                                 c("A", "D"),
                                 c("B", "C"),
                                 c("B", "D"),
                                 c("C", "D")),
              annotation=c("**"), 
              map_signif_level=T, 
              tip_length=c(0,0,0,0,0,0), 
              #y_position = c(46,54,49),
              size=1, 
              textsize = 7, 
              test = "t.test") 