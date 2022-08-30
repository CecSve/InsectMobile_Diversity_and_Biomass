#### load libraries ####

# importing and data wrangling
library(tidyverse)
library(readr)
library(data.table)

# stats
library(iNEXT)
library(vegan)
library(phyloseq)
library(fossil)
library(lme4)
library(lmerTest)
library(MuMIn)

# visualisation
library(sjPlot)
library(ggplot2)
library(cowplot)
library(ggpubr)

#### load data ####

asvs <- readRDS("data/cleaned_data/asvs_filtered.rds")
taxonomy <- readRDS("data/cleaned_data/taxonomy_filtered.rds")
data <- readRDS("data/cleaned_data/combData_noSizeFrac.rds")

#### Explore the data ####

#how many routes
length(unique(data$RouteID_JB))

# how many samples
length(unique(data$SampleID)) # total/complete samples

# how many pilots?
length(unique(data$PID))

# How many unique in each taxonomic level? 

taxonomy %>% summarise_all(n_distinct) 

# order
taxonomy %>% drop_na(order) %>% summarise_all(n_distinct) 

taxonomy %>%
  group_by(order) %>%
  drop_na(order) %>% 
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>% 
  arrange(desc(n))

# family
taxonomy %>% drop_na(family) %>% summarise_all(n_distinct) 

taxonomy %>%
  group_by(family) %>%
  drop_na(family) %>% 
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>% 
  arrange(desc(n))

# genus
taxonomy %>% drop_na(genus) %>% summarise_all(n_distinct) 

taxonomy %>%
  group_by(genus) %>%
  drop_na(genus) %>% 
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>% 
  arrange(desc(n))

# species
taxonomy %>% drop_na(species) %>% summarise_all(n_distinct) 

taxonomy %>%
  group_by(species) %>%
  drop_na(species) %>% 
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>% 
  arrange(desc(n))

table(taxonomy$order)

# how frequent are the different match types
taxonomy %>%
  group_by(matchtype) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n))

#### Adding richness and diversity variables ####

pa_asvs <- decostand(asvs, method = "pa") # transform asv table into presence absence 

#Data has species as rows and sites as columns, but we need the opposite. Thus, transpose the original data
tpa_asvs <- t(pa_asvs) 

#### iNEXT (Hill numbers) ####

freq_asvs <- iNEXT::as.incfreq(pa_asvs)# transform incidence raw data (a species by sites presence-absence matrix) to incidence frequencies data (iNEXT input format, a row-sum frequencies vector contains total number of sampling units)

# set a series of sample sizes (m) for R/E computation
t <- seq(1, 4000, by=50)

#apply `iNEXT` main function
asvspa.inext <- iNEXT(freq_asvs, q = 0, datatype = "incidence_freq", size = t) 

asvspa.inext$DataInfo # summarizing data information, returns basic data information including the reference sample size (n), observed species richness (S.obs), a sample coverage estimate (SC), and the first ten frequency counts (f1‐f10)

asvspa.inext$iNextEst # showing diversity estimates along with related statistics for a series of rarefied and extrapolated samples
asvspa.inext$AsyEst # showing asymptotic diversity estimates along with related statistics

ChaoRichness(freq_asvs, datatype = "incidence_freq", conf = 0.95)

#look at the data
asvspa.inext

# Sample-size-based R/E curves without figure legend
ggiNEXT(asvspa.inext, type=1) +
  theme_bw(base_size = 18) + theme(legend.position="none")

# Sample completeness curves without figure legend
ggiNEXT(asvspa.inext, type=2) +ylim(c(0.9,1)) +
  theme_bw(base_size = 18) + theme(legend.position="none")

# Coverage-based R/E curves with legend placed at the bottom, where “Guides” and “Method” are left out
ggiNEXT(asvspa.inext, type=3) + xlim(c(0.9,1)) +
  theme_bw(base_size = 18) +
  theme(legend.position="bottom", legend.title=element_blank())
