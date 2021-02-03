########################################
###   PERMANOVA and SIMPER analyses  ###
### Yue Shi,University of Washington ###  
########################################

library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(vegan)
library(lubridate)

rm(list=ls())

###Load data
#meta
meta <- read.csv("./data/meta_N202.csv", header=T)
meta$Year <- as.numeric(format(as.Date(meta$Date, format="%m/%d/%Y"), "%Y")) + 2000
meta$Year_Season <- paste0(meta$Year, meta$Season)
table(meta$Year_Season, meta$Pack, meta$Predator.species.confirmed)

#presence data
presence.all <- read.csv("./data/canid.diet.presence.long.csv", header=T)

presence.all.meta.prey <- 
  left_join(presence.all,meta,by=c("sampleID"="SampleID")) %>% 
  filter(taxonomy!="Wolf or Coyote")

prey.order=c("Cow", "Moose", "Elk", "Pig", "Deer", 
             "Rabbit", "Snowshoe hare", "Muskrat", "Ground squirrel","Red squirrel", 
             "Flying squirrel", "Chipmunk", "Meadow vole", "Red-backed vole", "Deer mouse",
             "Wild turkey","Ruffed grouse", "Spruce grouse", "European starling")

### Prepare for the ecological analyses

#Check samples have no prey information, but only predator information. These samples need to be removed for the following analyses
presence.all.meta.prey.spread <- 
  presence.all.meta.prey %>% 
  group_by(sampleID) %>% 
  mutate(presenceSum=sum(presence)) %>% 
  filter(presenceSum>0) %>% 
  spread(taxonomy,presence) 

prey.names <- unique(presence.all.meta.prey$taxonomy)

prey.mat <- 
  presence.all.meta.prey.spread %>% 
  ungroup() %>% 
  select(prey.names) %>% 
  as.matrix

presence.all.meta.prey.spread.wolf <- 
  presence.all.meta.prey.spread %>% 
  ungroup() %>% 
  filter(Predator.species.confirmed=="Wolf") #94 samples

presence.all.meta.prey.spread.coyote <- 
  presence.all.meta.prey.spread %>% 
  ungroup() %>% 
  filter(Predator.species.confirmed=="Coyote")  #98 samples

### Diet ~ species + pack + season

prey.dist=vegdist(prey.mat,method="jaccard", binary=TRUE) #jaccardd is more suitable for presence/absence data
set.seed(29)
#test them together, by="margin" or by "terms"
adonis2(prey.dist ~ Predator.species.confirmed + Pack + Season, by="margin",
        data=presence.all.meta.prey.spread,
        permutations = 999, method="jaccard")

adonis2(prey.dist ~ Predator.species.confirmed + Pack + Season , by="terms",
        data=presence.all.meta.prey.spread,
        permutations = 999, method="jaccard")

#interspecific dietary difference is the primary source of variation, followed by pack and season. 

#SIMPER test to see which prey taxa contribute the most to the difference
rownames(prey.mat)=presence.all.meta.prey.spread$sampleID
prey.meta=as.data.frame(presence.all.meta.prey.spread[,1:9])

(sim.sp <- with(prey.meta, simper(prey.mat, Predator.species.confirmed, permutation=999)))
sim.sp$Coyote_Wolf$overall
summary(sim.sp, ordered=TRUE)


#### Spatialtemporal dietary difference in wolf

prey.wolf.mat <- 
  presence.all.meta.prey.spread.wolf %>% 
  select(prey.names) %>% 
  as.matrix
prey.wolf.dist <- vegdist(prey.wolf.mat,method="jaccard", binary=TRUE) #jaccardd is more suitable for presence/absence data

set.seed(29)
adonis2(prey.wolf.dist~Pack*Season, data=presence.all.meta.prey.spread.wolf, 
        permutations = 999, method="jaccard")

rownames(prey.wolf.mat)=presence.all.meta.prey.spread.wolf$sampleID
prey.wolf.meta=as.data.frame(presence.all.meta.prey.spread.wolf[,1:9])


(sim.wolf.pack <- with(prey.wolf.meta, simper(prey.wolf.mat, Pack, permutation=999)))
sim.wolf.pack$`Goodman Meadows_Smackout`$overall
sim.wolf.pack$`Goodman Meadows_Dirty Shirt`$overall
sim.wolf.pack$`Smackout_Dirty Shirt`$overall
sim.wolf.pack
summary(sim.wolf.pack, ordered=TRUE)


#### Spatialtemporal dietary difference in coyote
prey.coyote.mat <- 
  presence.all.meta.prey.spread.coyote %>% 
  select(prey.names) %>% 
  as.matrix
prey.coyote.dist=vegdist(prey.coyote.mat,method="jaccard", binary=TRUE) #jaccard is more suitable for presence/absence data

adonis2(prey.coyote.dist~Pack*Season, data=presence.all.meta.prey.spread.coyote, 
        permutations = 999, method="jaccard")

rownames(prey.coyote.mat)=presence.all.meta.prey.spread.coyote$sampleID
prey.coyote.meta=as.data.frame(presence.all.meta.prey.spread.coyote[,1:9])


#subset into packs;
prey.coyote.gm.meta <- prey.coyote.meta %>% 
  filter(Pack=="Goodman Meadows")
prey.coyote.ds.meta <- prey.coyote.meta %>% 
  filter(Pack=="Dirty Shirt")
prey.coyote.sm.meta <- prey.coyote.meta %>% 
  filter(Pack=="Smackout")

prey.coyote.gm.mat <- prey.coyote.mat[rownames(prey.coyote.mat) %in% prey.coyote.gm.meta$sampleID, ]
prey.coyote.ds.mat <- prey.coyote.mat[rownames(prey.coyote.mat) %in% prey.coyote.ds.meta$sampleID, ]
prey.coyote.sm.mat <- prey.coyote.mat[rownames(prey.coyote.mat) %in% prey.coyote.sm.meta$sampleID, ]

(sim.coyote.gm.season <- with(prey.coyote.gm.meta, simper(prey.coyote.gm.mat, Season, permutation=999)))
sim.coyote.gm.season$Spring_Fall$overall
summary(sim.coyote.gm.season, ordered=TRUE)

(sim.coyote.ds.season <- with(prey.coyote.ds.meta, simper(prey.coyote.ds.mat, Season, permutation=999)))
sim.coyote.ds.season$Spring_Fall$overall
summary(sim.coyote.ds.season, ordered=TRUE)

(sim.coyote.sm.season <- with(prey.coyote.sm.meta, simper(prey.coyote.sm.mat, Season, permutation=999)))
sim.coyote.sm.season$Spring_Fall$overall
summary(sim.coyote.sm.season, ordered=TRUE)

#subset into seasons;
prey.coyote.sp.meta <- prey.coyote.meta %>% 
  filter(Season=="Spring")
prey.coyote.fa.meta <- prey.coyote.meta %>% 
  filter(Season=="Fall")

prey.coyote.sp.mat <- prey.coyote.mat[rownames(prey.coyote.mat) %in% prey.coyote.sp.meta$sampleID, ]
prey.coyote.fa.mat <- prey.coyote.mat[rownames(prey.coyote.mat) %in% prey.coyote.fa.meta$sampleID, ]

(sim.coyote.sp.pack <- with(prey.coyote.sp.meta, simper(prey.coyote.sp.mat, Pack, permutation=999)))
sim.coyote.sp.pack$`Smackout_Goodman Meadows`$overall
sim.coyote.sp.pack$`Goodman Meadows_Dirty Shirt`$overall
sim.coyote.sp.pack$`Smackout_Dirty Shirt`$overall
summary(sim.coyote.sp.pack, ordered=TRUE)


(sim.coyote.fa.pack <- with(prey.coyote.fa.meta, simper(prey.coyote.fa.mat, Pack, permutation=999)))
sim.coyote.fa.pack$`Goodman Meadows_Smackout`$overall
sim.coyote.fa.pack$`Goodman Meadows_Dirty Shirt`$overall
sim.coyote.fa.pack$`Smackout_Dirty Shirt`$overall
summary(sim.coyote.fa.pack, ordered=TRUE)


### Figure 4
## prepare 19 distinct colors
val=colorRampPalette(brewer.pal(11,"Spectral"))(29)
val.red=val[seq(1,9,2)]
val.yellow=val[seq(11,20,1)]
val.blue=val[seq(23,29,2)]
val.select=c(val.red,val.yellow,val.blue)

# WOLF data, only do packs
wolf.pack <- 
  presence.all.meta.prey %>% 
  filter(Predator.species.confirmed=="Wolf") %>% 
  group_by(Pack, taxonomy) %>% 
  summarise(presence=sum(presence))

wolf.pack.count <- meta %>% 
  filter(Predator.species.confirmed=="Wolf") %>% 
  group_by(Pack) %>% 
  mutate(n=n()) %>% 
  select(Pack,n) %>% 
  unique()

wolf <- left_join(wolf.pack, wolf.pack.count, by="Pack") %>% 
  mutate(freq=presence/n)

wolf$taxonomy=as.factor(wolf$taxonomy)
wolf$taxonomy=factor(wolf$taxonomy,levels=prey.order)

p.wolf.pack <- wolf %>% 
  ggplot(aes(x=taxonomy, y=freq, fill=taxonomy))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand = c(0, 0), limits=c(0, 0.9))+
  scale_fill_manual(values=val.select)+
  facet_grid(cols=vars(Pack)) +
  guides(fill = guide_legend(ncol=1)) +
  theme_bw(base_size=15) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 15))+
  ylab("Frequency of Occurrence")+
  ggtitle("Wolf")

# COYOTE data
coyote <- 
  presence.all.meta.prey %>% 
  filter(Predator.species.confirmed=="Coyote") %>% 
  group_by(Pack, Season, taxonomy) %>% 
  summarise(presence=sum(presence))

count.coyote <- meta %>% 
  filter(Predator.species.confirmed=="Coyote") %>% 
  group_by(Pack, Season) %>% 
  mutate(n=n()) %>% 
  select(Pack, Season, n) %>% 
  unique()

coyote <- left_join(coyote, count.coyote, by=c("Pack", "Season")) %>% 
  mutate(freq=presence/n)

coyote$taxonomy=as.factor(coyote$taxonomy)
coyote$taxonomy=factor(coyote$taxonomy,levels=prey.order)
coyote$Season=as.factor(coyote$Season)
coyote$Season=factor(coyote$Season, levels=c("Spring","Fall"))

p.coyote.both <- coyote %>% 
  ggplot(aes(x=taxonomy, y=freq, fill=taxonomy))+
  geom_bar(stat="identity")+
  scale_y_continuous(expand = c(0, 0), limits=c(0, 0.9))+
  scale_fill_manual(values=val.select)+
  facet_grid(cols=vars(Pack), rows=vars(Season)) +
  guides(fill = guide_legend(ncol=1)) +
  theme_bw(base_size=15) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 15))+
  ylab("Frequency of Occurrence")+
  ggtitle("Coyote")

#pdf("./figures/Fig4_spatialTemp_wolf_coyote.pdf", width = 12, height=8)
p.wolf.pack / p.coyote.both + 
  plot_layout(heights = c(1,2))+
  plot_layout(guides="collect") & theme(legend.position = "right")
#dev.off()

#Summary
wolf %>% 
  select(-presence, -n) %>% 
  spread(Pack, freq) %>% 
  arrange(taxonomy) 

coyote %>% 
  select(-presence, -n) %>% 
  filter(Pack=="Dirty Shirt") %>% 
  spread(Season, freq) %>% 
  arrange(taxonomy)

coyote %>% 
  select(-presence, -n) %>% 
  filter(Season=="Spring") %>% 
  spread(Pack, freq) %>% 
  arrange(taxonomy)


