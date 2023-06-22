##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

##R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##Copyright (C) 2021 The R Foundation for Statistical Computing
##Platform: x86_64-apple-darwin17.0 (64-bit)


#Heatmap.2 visualization of Genes -- hierarchically cluster genes in the 3 program overlap list (396), visualize genes expression clusters with a heatmap,
#visualize each cluster using a multi-line graph, and run gene ontology analysis on the 3 program overlap list and each respective cluster from the hclust)    

## TPMS = Transcripts per million 

#####################
library(ggplot2)

## EDGE R 
ED_DEG <- data.frame( 
  Time_Interval = rep(c('U-00','00-01','01-03','03-06','06-12','12-24','24-48'), each = 2),
  DEG_count = c(27,0,38,-4,63,-175,127,-120,4,-1,0,0,2,-1),
  Direction = rep(c("up","down"),7))
  
ggplot(ED_DEG, aes(x = Time_Interval, y = DEG_count))+
  geom_col(aes(fill = Direction))+
  theme_classic()+
   scale_x_discrete(limits=ED_DEG$Time_Interval)+
scale_y_continuous(breaks = seq(-300,300, by=20))


#NOISEQBIO  
NOIS_DEG <- data.frame( 
  Time_Interval = rep(c('U-00','00-01','01-03','03-06','06-12','12-24','24-48'), each = 2),
  NOIS_DEG_count = c(36,-34,34,-26,1561,-1257,3699,-2223,99,-69,45,-67,127,-119),
  Direction = rep(c("up","down"),7))

ggplot(NOIS_DEG, aes(x = Time_Interval, y = NOIS_DEG_count))+
  geom_col(aes(fill = Direction) )+
  theme_classic()+
  scale_x_discrete(limits=NOIS_DEG$Time_Interval)+
  scale_y_continuous(limits = c(-5000,5000), breaks = seq(-5000,5000, by=500 ))+
  scale_fill_discrete(breaks=c("up", "down"))
  

 =
