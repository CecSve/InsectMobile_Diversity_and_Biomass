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
library(Rarefy)
library(breakaway) # for calculating estimated richness instead of rarefied richness

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

# how many bins are species names associated with (ordered descending)
taxonomy %>%
  drop_na(species) %>% 
  group_by(species) %>%
  summarise(bins = n_distinct(scientificName)) %>% 
  arrange(desc(bins))

table(taxonomy$order) # several wingless orders in there - potential prey/gut content?

# how frequent are the different match types across sequences 
taxonomy %>%
  group_by(matchType) %>%
  summarise(asvs = n()) %>%
  mutate(freq = asvs/sum(asvs))

# how frequent are the different match types across BINs
taxonomy %>%
  group_by(matchType) %>%
  summarise(bins = n_distinct(scientificName)) %>% 
  mutate(freq = bins/sum(bins))

#### Adding richness and diversity variables ####

pa_asvs <- decostand(asvs, method = "pa") # transform asv table into presence absence 

#Data has species as rows and sites as columns, but we need the opposite. Thus, transpose the original data
tpa_asvs <- t(pa_asvs) 
tpa_asvs[1:5,1:5] # samples as rows, asvs as columns
tasvs <- t(asvs) # read abundance transposed table
pa_asvs[1:5,1:5] # asvs are rows, samples as columns 

##### vegan ####

# get minimum number of reads for rarefying
quantile(rowSums(dtasvs), probs = seq(0, 1, 0.2)) # first see how many samples are associated with how many reads, setting it to min is quite low in this case, but at least we do not remove any more samples
min_n_seqs <- min(rowSums(dtasvs))
dtasvs <- as.data.frame(tasvs)
str(dtasvs)

# rarefy each sample to the minimum sequence depth  - sequences in the samples are now distributed across their associated ASVs summarizing the sample abundances to the minimum reads (min_n_seq)

###### rarefied richness per sample (number of taxa per sample) ####
richness <- dtasvs %>% 
  rarefy(min_n_seqs) %>% 
  as_tibble(rownames = "SampleID") %>% 
  select(SampleID, richness = value)

###### rarefied Shannon diversity ####
# make a funtion for iterating the Sannon diversity per sample calculation
shannon_iteration <- function(){
  dtasvs %>% 
    rrarefy(sample = min_n_seqs) %>% 
    diversity()
  }

# iterate a 1000 times and extract mean Shannon diversity
shannon <- replicate(100, shannon_iteration()) %>% 
  as_tibble(rownames = "SampleID", .name_repair = "unique") %>% 
  pivot_longer(-SampleID) %>% 
  group_by(SampleID) %>% 
  summarise(shannon = mean(value))

# rarecurve - this step takes a long time
rarecurve_data <- rarecurve(dtasvs, step = 100)

map_dfr(rarecurve_data, bind_rows) %>% 
  bind_cols(SampleID = rownames(dtasvs), .) %>% 
  pivot_longer(-SampleID) %>% 
  drop_na() %>% 
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name) %>% 
  ggplot(aes(x=n_seqs, y = value, group = SampleID)) + 
  geom_vline(xintercept = min_n_seqs, color = "gray") + # add min_n_seq indication
  geom_line() + 
  theme_classic()

##### add richness and diversity to data ####
data_richness <- inner_join(data, richness, by = "SampleID") %>% 
  inner_join(., shannon, by = "SampleID")

saveRDS(data_richness, file = "data/cleaned_data/data_richness.RDS")

#### iNEXT (Hill numbers) ####

freq_asvs <- iNEXT::as.incfreq(pa_asvs)# transform incidence raw data (a species by sites presence-absence matrix) to incidence frequencies data (iNEXT input format, a row-sum frequencies vector contains total number of sampling units)

# compute species diversity (Hill numbers with q = 0, 1 and 2) with a particular user-specified level of sample size or sample coverage.
out <- estimateD(pa_asvs, q = c(0,1,2), datatype = "incidence_raw", base="coverage", level=0.985, conf=NULL, nboot = 5)
out

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
