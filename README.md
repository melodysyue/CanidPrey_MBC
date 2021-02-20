# CanidPrey_MBC
 
# How to cite:
Shi, Y., Hoareau, Y., Reese, E.M. et al. Prey partitioning between sympatric wild carnivores revealed by DNA metabarcoding: a case study on wolf (Canis lupus) and coyote (Canis latrans) in northeastern Washington. Conserv Genet (2021). https://doi.org/10.1007/s10592-021-01337-2

# Data

**./nc_raw/**: this directory contains raw sequence data from PCR products of negative controls (N=100), including extraction and PCR negative controls; 

**canid.diet.obitools.tab**: read count table of 374 MOTUs after obitools filtering;

**MOTUS_N374.fasta**: fasta file of 374 MOTUs; 

**taxonomy_MOTUs_N348.csv**: taxonomy classification of 348 MOTUs. Note: 26 out of 374 MOTUs were removed because they could not be assigned to any taxon with at least 98% identity.

**meta_N202.csv**: meta information of 202 samples; 

**canid.diet.postFiltering.tab** : final MOTU read count table, consisting of 332 MOTUs and 1037 PCR replicates across 202 samples. 

**canid.diet.postFiltering_long.csv**: final MOTU read count table in the long format along with meta information; 

**canid.diet.presence.long.csv**: final presence/absence data;

**byspecies.freq.csv**: frequency of occurrence data of 19 prey items in wolf and coyote;

# Scripts

**1_obitools.sh**: sequence data processing using the OBITools package;

**2_filtering.r**: data filtering to control for potential contamination and generate the final MOTU read count table;

**3_blocking.r**: assess the effects of predator-specific blocking primer and convert read count data to presence/absence data. Output: Fig.2 and Supplementary Table 1;

**4_dietByspecies.r**: assess the interspecific dietary differences between wolves and coyotes. Output: Fig.3;

**5_diet_sig.r**: PERMANOVA and SIMPER tests. Output: Fig.4, Supplementary Table 2 and Supplementary Table 3;

**6_rra.r**: Comparison of dietary profiles using FOO vs. RRA. Output: Supplementary Fig. 1
