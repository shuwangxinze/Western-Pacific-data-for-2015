#CAP
library(vegan)
library(tidyverse)
library(Hmisc)
#species data
otu <- read.delim(file.choose(), row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
otu.hell <- decostand(otu, "hellinger")

#Read environment data
env <- read.delim(file.choose(), row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
env <- env[c("Salinity","Temperature","SSTA","DO","PO4","SiO3","NO2","NO3","NH4")]
db_rda <- capscale(otu~., env, distance = 'bray', add = TRUE)
summary(db_rda)

anova(db_rda,by = "term")#Significance of each variable

anova(db_rda,by = "axis")#Significance of each axis

#Extraction canonical coefficient
rda_coef <- coef(db_rda)

#R2 correction
r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared #Original R2
db_rda_adj <- r2$adj.r.squared  #The corrected R2

#permutation test
#The permutation test for all constrained axes, the global test, is based on 999 permutations
cap_term <- anova.cca(cap, permutations = 999,by = "term")

#Each constraint axis is checked one by one, based on 999 permutations
db_rda_test_axis <- anova.cca(cap, by = 'axis', permutations = 999)

#P-value correction 
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'bonferroni')
#Full model P-value correction
cap_test$`Pr(>F)` <- p.adjust(cap_test$`Pr(>F)`, method = 'bonferroni')

#visualization
library(ggplot2)

scrs<-scores(cap)

names(scrs)

#Calculate the arrow multiple and calculate the length of the arrow

multiplier <- vegan:::ordiArrowMul(scrs$biplot)

multiplier

df_arrows<- scrs$biplot*multiplier#Calculate the length of the arrow

df_arrows

colnames(df_arrows)

rownames(df_arrows)

#Add quadrat coordinates

site<-scrs$sites %>%
  
  as_tibble(rownames = "sample")

site

group <- read.delim(file.choose(), sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

site2 <- site %>%
  
  left_join(group, by = "sample")

df_arrows=as.data.frame(df_arrows)


ggplot(data=site2,aes(CAP1,CAP2,color = group)) +
  
  scale_color_manual(values = c("#3C5488FF","#de572d","#FFBE7A","#6dc889"),)+#Dot color
  
  geom_hline(aes(yintercept=0),colour="#d8d6d6",linetype=5)+
  
  geom_vline(aes(xintercept=0),colour="#d8d6d6",linetype=5)+
  
  geom_point(size=3.5,shape=16) +
  geom_text(aes(x = CAP1-0.05,y = CAP2+0.05 ,label = rownames(spe)))+
  
  geom_segment(data=as.data.frame(df_arrows), aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               
               arrow = arrow(angle= 15,length = unit(0.25, "cm")),color="#000000",alpha=0.5)+#Draw arrows
  
  geom_text(data=as.data.frame(df_arrows*1.1),aes(CAP1, CAP2, label = rownames(df_arrows)),
            
            color="#000000",alpha=0.7)+#Arrow information
  theme_bw() +
  
  theme(panel.grid=element_blank())+
  
  labs(x = "CAP1 (XX%)", y = "CAP2 (XX%)")