########################################
### Data filtering on the MOTU table ###
### Yue Shi,University of Washington ###  
########################################

rm(list=ls())
library(tidyverse)

#### Load the files
#read count table
count=read.table("./data/canid.diet.obitools.tab", header=T) #374 MOTUs in total.
count <- count %>% 
  select(-sequence)

#taxonomy assignment
taxa=read.csv("./data/taxonomy_MOTUs_N348.csv", header=T)
dim(taxa) #only 348 MOTUs were assigned taxonomy, 26 MOTUs were discarded. 
taxa <- taxa %>% 
  select(query, taxonomy)

#sample meta information
meta=read.csv("./data/meta_N202.csv",header=TRUE) #202 samples. 


# organize the data
ng.names <- count %>%
  select(starts_with("sample.pcr"), starts_with("sample.Extr")) %>% 
  names() 

sp.names <- count %>% 
  select(starts_with("sample.")) %>% 
  select(-starts_with("sample.pcr")) %>% 
  select(-starts_with("sample.Extr")) %>% 
  names() 

#combine count table and taxonomy information
count.taxa <- left_join(taxa, count, by=c("query"="id")) 

#####################################
### Step 1: Remove contamination ###
#####################################

# First remove any MOTU whose maximum read count occur in any negative control. 

count.taxa %>% 
  mutate(ngmax = pmax(!!!rlang::syms(ng.names))) %>% 
  mutate(spmax = pmax(!!!rlang::syms(sp.names))) %>% 
  select(query, taxonomy, count, ngmax, spmax) %>% 
  filter(ngmax>=spmax) %>% 
  summary()

count.taxa %>% 
  mutate(ngmax = pmax(!!!rlang::syms(ng.names))) %>% 
  mutate(spmax = pmax(!!!rlang::syms(sp.names))) %>% 
  select(query, taxonomy, count, ngmax, spmax) %>% 
  filter(ngmax>=spmax) %>% 
  group_by(taxonomy) %>% 
  tally()

# 13 MOTUs were removed. 

count.taxa.rm <- count.taxa %>% 
  mutate(ngmax = pmax(!!!rlang::syms(ng.names))) %>% 
  mutate(spmax = pmax(!!!rlang::syms(sp.names))) %>% 
  filter(ngmax<spmax)

# Second, among the remaining 335 MOTUs, subtract the highest read count among all negative controls from the abundance of each sample PCR replicate.
count.taxa.rm.noctm <- count.taxa.rm %>% 
  select(ngmax, starts_with("sample")) %>% 
  select(-starts_with("sample.pcr"),-starts_with("sample.Extr")) %>% 
  as.matrix #sweep only works with matrix

count.taxa.rm.noctm <- as.data.frame(sweep(count.taxa.rm.noctm[,2:ncol(count.taxa.rm.noctm)],1,count.taxa.rm.noctm[,1],"-"))
count.taxa.rm.noctm[count.taxa.rm.noctm<0]=0 #replace negative values with 0;

count.taxa.rm.noctm$id=count.taxa.rm$query
count.taxa.rm.noctm<- 
  count.taxa.rm.noctm %>% 
  select(id,everything()) #add id back

# Third, suppress read count to zero if its relative abundance (abundance in PCR/total abundance across PCRs) is too low, alpha=0.03%.

id=count.taxa.rm.noctm[,1]
count.taxa.rm.noctm=count.taxa.rm.noctm[,-1]
rownames(count.taxa.rm.noctm)=id
cutoff <- round(rowSums(count.taxa.rm.noctm)*0.0003)
count.taxa.rm.noctm.cutoff <- count.taxa.rm.noctm
count.taxa.rm.noctm.cutoff[count.taxa.rm.noctm.cutoff<=cutoff]=0
count.taxa.rm.noctm.cutoff=cbind(id,count.taxa.rm.noctm.cutoff)

# add taxonomy back 
count.taxa.rm.noctm.cutoff.taxa <- 
  left_join(count.taxa.rm.noctm.cutoff, taxa, by=c("id"="query")) %>% 
  select(id, taxonomy, everything())

# remove any MOTUs assigned to unlikely prey, and any pcr replicates with zero count after filtering
table(count.taxa.rm.noctm.cutoff.taxa$taxonomy)

count.taxa.rm.noctm.cutoff.taxa %>% 
  filter(taxonomy=="Human" | taxonomy=="Cougar or Bobcat") %>% 
  mutate(count=rowSums(select(., starts_with("sample.")))) %>% 
  select(id, taxonomy, count) 

 count.taxa.rm.noctm.cutoff.taxa %>% 
  filter(taxonomy!="Human") %>% 
  filter(taxonomy!="Cougar or Bobcat") %>% 
  select(-taxonomy) %>% 
  column_to_rownames("id") %>% 
  select_if(~ !is.numeric(.) || sum(.)==0) %>% 
  dim() #9 sample PCR replicates with zero read count after filtering were removed;

 count.final <- 
   count.taxa.rm.noctm.cutoff.taxa %>% 
   filter(taxonomy!="Human") %>% 
   filter(taxonomy!="Cougar or Bobcat") %>% 
   select(-taxonomy) %>% 
   column_to_rownames("id") %>% 
   select_if(~ !is.numeric(.) || sum(.)!=0) %>% 
   rownames_to_column("id")

#write.table(count.final, "./data/canid.diet.postFiltering.tab", quote=F, row.names = F)

# read count distribution among samples after filtering
sp.count <- count.final %>% 
  column_to_rownames("id") %>% 
  colSums() %>% 
  as.data.frame()
names(sp.count)="count"

sum(sp.count$count)
summary(sp.count)


#########################################
### Step 2: summarize by prey species ###
#########################################
#make a long table
count.final.long <- 
  count.taxa.rm.noctm.cutoff.taxa %>% 
  filter(taxonomy!="Human") %>% 
  filter(taxonomy!="Cougar or Bobcat") %>% 
  gather(key=sampleName, value=readCount,-c("id","taxonomy"))

count.final.long$sampleName = gsub("sample.", "",count.final.long$sampleName)

count.final.merge <- 
  count.final.long%>% 
  group_by(sampleName, taxonomy) %>% #sum read count from multiple MOTUs belonging to the same prey taxa
  summarise(readCount=sum(readCount))

# add additional information columns
count.final.merge$sampleName=gsub("redo","",count.final.merge$sampleName)
count.final.merge$sampleID=gsub("[^0-9]","",count.final.merge$sampleName)
count.final.merge$replicateID=gsub("[[:digit:]]","",count.final.merge$sampleName)
count.final.merge$Predator=count.final.merge$taxonomy=="Wolf or Coyote"
count.final.merge$Blocking=count.final.merge$replicateID %in% c("D","E","F")

count.final.merge$Predator=as.factor(count.final.merge$Predator)
levels(count.final.merge$Predator)=c("Prey","Predator")
count.final.merge$Blocking=as.factor(count.final.merge$Blocking)
levels(count.final.merge$Blocking)=c("Without blocking primer","With blocking primer")

# remove sample replicate with 0 read count
keep <- count.final.merge %>% 
  select(-taxonomy, -Predator) %>% 
  group_by(sampleName) %>% 
  mutate(readCountSum=sum(readCount)) %>% 
  filter(readCountSum>0) %>% 
  pull(sampleName) %>% 
  unique()

count.final.merge <- count.final.merge %>% 
  filter(sampleName %in% keep) #1037 sample PCR replicates remained;

#write.csv(count.final.merge, "./data/canid.diet.postFiltering_long.csv", quote=F, row.names = F)
