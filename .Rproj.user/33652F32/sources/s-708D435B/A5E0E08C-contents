################################################################
### Dietary differences between sympatric wolves and coyotes ###
###            Yue Shi,University of Washington              ###  
################################################################
library(tidyverse)
library(RColorBrewer)
library(randomcoloR)
library(cowplot)

rm(list=ls())

###Load data
#meta
meta <- read.csv("./data/meta_N202.csv", header=T)
#presence data
presence.all <- read.csv("./data/canid.diet.presence.long.csv", header=T)

presence.all.meta.prey <- 
  left_join(presence.all,meta,by=c("sampleID"="SampleID")) %>% 
  filter(taxonomy!="Wolf or Coyote")

prey.order=c("Cow", "Moose", "Elk", "Pig", "Deer", 
             "Rabbit", "Snowshoe hare", "Muskrat", "Ground squirrel","Red squirrel", 
             "Flying squirrel", "Chipmunk", "Meadow vole", "Red-backed vole", "Deer mouse",
             "Wild turkey","Ruffed grouse", "Spruce grouse", "European starling")

### Diet comparison by species

#Turn presence/absence data into frequency of occurrence 
bysp <- presence.all.meta.prey %>% 
  group_by(Predator.species.confirmed,taxonomy) %>% 
  summarise(presence=sum(presence)) %>% 
  spread(Predator.species.confirmed, presence)

bysp.freq <- bysp
bysp.freq$Wolf=bysp.freq$Wolf/99
bysp.freq$Coyote=bysp.freq$Coyote/103

bysp.freq <- bysp.freq %>% 
  gather(key="predator",value="freq",-taxonomy) 


bysp.freq$taxonomy=as.factor(bysp.freq$taxonomy)
bysp.freq$taxonomy=factor(bysp.freq$taxonomy,levels=prey.order)
bysp.freq$predator=as.factor(bysp.freq$predator)
bysp.freq$predator=factor(bysp.freq$predator, levels=c("Wolf","Coyote"))


#Fig 3
## prepare 19 distinct colors
val=colorRampPalette(brewer.pal(11,"Spectral"))(29)
val.red=val[seq(1,9,2)]
val.yellow=val[seq(11,20,1)]
val.blue=val[seq(23,29,2)]
val.select=c(val.red,val.yellow,val.blue)


# plot it
p1 <-bysp.freq %>% 
  filter(predator=="Wolf") %>% 
  ggplot(aes(x=taxonomy,y=freq,fill=taxonomy))+
  geom_bar(stat="identity")+
  theme_classic(base_size=15)+
  scale_fill_manual(values=val.select)+
  scale_x_discrete(limits = rev(levels(bysp.freq$taxonomy)))+
  scale_y_reverse(expand=c(0,0),limits = c(0.7, 0), breaks = seq(0.7, 0, -0.1))+
  ylab("Frequency of Occurrence\nWolf N=99")+
  coord_flip()+
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none",
        legend.title = element_blank())


p2 <- bysp.freq %>% 
  filter(predator=="Coyote") %>% 
  ggplot(aes(x=taxonomy,y=freq,fill=taxonomy))+
  geom_bar(stat="identity")+
  theme_classic(base_size=15)+
  scale_fill_manual(values=val.select)+
  scale_x_discrete(limits = rev(levels(bysp.freq$taxonomy)))+
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,0.7,0.1), limits=c(0,0.7))+
  ylab("Frequency of Occurrence\nCoyote N=103")+
  coord_flip()+
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none",
        legend.title = element_blank())

p3 <- bysp.freq %>% 
  filter(predator=="Wolf") %>% 
  ggplot(aes(x=taxonomy))+
  geom_text(aes(y=0, label=taxonomy), size=5)+
  theme_classic(base_size=15)+
  scale_x_discrete(limits = rev(levels(bysp.freq$taxonomy)))+
  coord_flip()+
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        legend.position="none",
        legend.title = element_blank())

p=plot_grid(p1,p3,p2,lables=NULL, ncol = 3,nrow=1, rel_widths = c(4,1.5,4), align="hv")
p
#ggsave(plot=p,"./figures/Fig3_dietbyspecies_butterfly.pdf",width=12,height=8)


