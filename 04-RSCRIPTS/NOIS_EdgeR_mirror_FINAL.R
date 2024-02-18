##Dorothy Mitchell ##
### Mapping with Bowtie2 Quantification with RSEM -- on Trimmed Reads ## 
## CydRegenSeq ## 

##R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##Copyright (C) 2021 The R Foundation for Statistical Computing
##Platform: x86_64-apple-darwin17.0 (64-bit)


#####################
library(ggplot2)

## EDGE R 
ED_DEG <- data.frame( 
  Time_Interval = rep(c('U-00','00-01','01-03','03-06','06-12','12-24','24-48'), each = 2),
  DEG_count = c(30,0,25,-3,25,-39,129,-115,5,-3,0,0,7,-2),
  Direction = rep(c("up","down"),7))
  
ggplot(ED_DEG, aes(x = Time_Interval, y = DEG_count))+
  geom_col(aes(fill = Direction))+
  theme_classic()+
   scale_x_discrete(limits=ED_DEG$Time_Interval)+
scale_y_continuous(breaks = seq(-300,300, by=20))


#NOISEQ  
NOIS_DEG <- data.frame( 
  Time_Interval = rep(c('U-00','00-01','01-03','03-06','06-12','12-24','24-48'), each = 2),
  NOIS_DEG_count = c(91,-57,50,-16,94,-340,3843,-2223,58,-79,45,-78,108,-62),
  Direction = rep(c("up","down"),7))

ggplot(NOIS_DEG, aes(x = Time_Interval, y = NOIS_DEG_count))+
  geom_col(aes(fill = Direction) )+
  theme_classic()+
  scale_x_discrete(limits=NOIS_DEG$Time_Interval)+
  scale_y_continuous(limits = c(-4000,4000), breaks = seq(-4000,4000, by=500 ))+
  scale_fill_discrete(breaks=c("up", "down"))
  


