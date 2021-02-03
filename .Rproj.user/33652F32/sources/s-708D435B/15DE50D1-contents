#############################################
### Assess effects of the blocking primer ###
###    Yue Shi,University of Washington   ###  
#############################################

rm(list=ls())
library(tidyverse)
library(RColorBrewer)
count.final.merge <- read.csv("./data/canid.diet.postFiltering_long.csv", header=T)
count.final.merge$Predator=factor(count.final.merge$Predator, levels=c("Prey","Predator"))
count.final.merge$Blocking=factor(count.final.merge$Blocking, levels=c("Without blocking primer","With blocking primer"))

# read freq comparison for Figure 2.
p_freq <- count.final.merge %>% 
  group_by(sampleName,Blocking) %>% 
  mutate(readCountSum=sum(readCount)) %>% 
  mutate(readfreq=readCount/readCountSum) %>%
  group_by(sampleName,Blocking,Predator) %>% 
  summarize(readfreq=sum(readfreq)) %>% 
  ggplot(aes(x=sampleName,y=readfreq,fill=Predator))+
  geom_bar(stat="identity", na.rm=TRUE)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1))+
  theme_classic(base_size=15)+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())+
  ylab("Proportion")+
  scale_fill_brewer(palette="Set1")+
  facet_grid(~Blocking)
#ggsave(plot=p_freq, "./figures/Fig2_blockingeffects.pdf")

#changes in read proportion of prey and predator
count.final.merge %>% 
  group_by(sampleName, Blocking, Predator) %>% 
  summarize(readCount=sum(readCount)) %>% 
  group_by(Blocking, Predator) %>% 
  summarize(sum=sum(readCount)) %>% 
  group_by(Blocking) %>% 
  mutate(prop=sum/sum(sum))

# Supplementary Table 1
# Changes in the number of positive PCR products
blk_rep <- count.final.merge %>%
  filter(readCount>0) %>% 
  group_by(Blocking, taxonomy) %>% 
  summarize(n=n()) %>% 
  spread(Blocking,n) 

# Changes in the average read count per positive PCR product
blk_avg <- count.final.merge %>% 
  filter(readCount>0) %>% 
  group_by(taxonomy, Blocking) %>% 
  summarise(avg=sum(readCount)/n()) %>% 
  spread(Blocking, avg)
  

### Turn read count data into presence/absence data
# a prey is considered as present if it occurred in at least 2 PCR replicates; 

presence.all <- count.final.merge %>%
  mutate(presence = ifelse(readCount>0, 1, 0)) %>% 
  group_by(sampleID,taxonomy) %>% 
  summarize(presence=sum(presence)) %>%
  mutate(presence=ifelse(presence>=2, 1, 0))

# distribution of the number of prey taxa across samples
presence.all%>% 
  filter(taxonomy!="Wolf or Coyote") %>% 
  group_by(sampleID) %>% 
  summarize(presence=sum(presence)) %>% 
  summary()

presence.all%>% 
  filter(taxonomy!="Wolf or Coyote") %>% 
  group_by(sampleID) %>% 
  summarize(presence=sum(presence)) %>% 
  group_by(presence) %>% 
  tally() 

#write.csv(presence.all, "./data/canid.diet.presence.long.csv", quote=F, row.names = F)
