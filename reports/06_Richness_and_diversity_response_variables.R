#### load libraries ####

# importing and data wrangling
library(tidyverse)
library(readr)
library(data.table)

# stats
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

#### explore the data ####

#how many routes
length(unique(data$RouteID_JB))

# how many samples
length(unique(data$SampleID)) # total/complete samples

# how many pilots?
length(unique(data$PID))

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

# for the analysis, we will only use sequences with an exact match to the reference database
taxonomy_allmatches <- taxonomy
taxonomy <- taxonomy_allmatches %>% filter(matchType == "BLAST_EXACT_MATCH") 
keep <- taxonomy$occurrenceId

asvs_allmatches <- asvs
asvs <- asvs_allmatches %>% filter(rownames(asvs_allmatches) %in% keep)
min(rowSums(asvs)) # some very low read sequences in there...
min(colSums(asvs))

# keep unique BINs and family-level combinations

nrow(unique(taxonomy[,c('scientificName','genus')]))
n_distinct(taxonomy$scientificName) # no difference between the two, so no BINs are assigned to the same family

asvs_allBINs <- asvs # keep the dataframe
forsubset <- taxonomy_allBINs %>% select(occurrenceId, scientificName) # to keep the samples, we need to include BINs temporarily in the AVS table
asvs <- asvs_allBINs %>% rownames_to_column("occurrenceId") 
asvs <- left_join(asvs, forsubset, by = "occurrenceId") # include BINs in the asv table

# keep only unique BINs while retaining reads in samples
asvs <- asvs %>%
  pivot_longer(!c(occurrenceId, scientificName),
               names_to = "SampleID",
               values_to = "read_count") %>%
  filter(read_count > 0) %>%
  distinct(scientificName, SampleID, .keep_all = T) %>%
  pivot_wider(names_from = SampleID,
              values_from = read_count,
              values_fill = 0) %>% 
  select(-scientificName) %>% 
  column_to_rownames(var = "occurrenceId")

keep <- rownames(asvs)
taxonomy_allBINs <- taxonomy
taxonomy <- taxonomy_allmatches %>% filter(occurrenceId %in% keep)

min(rowSums(asvs)) # some very low read sequences in there...
min(colSums(test)) 

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

taxonomy %>%
  group_by(order) %>%
  drop_na(order) %>% 
  summarise(bins = n_distinct(scientificName)) %>% 
  mutate(freq = bins/sum(bins)) %>% 
  arrange(desc(bins))

# family
taxonomy %>% drop_na(family) %>% summarise_all(n_distinct) 

taxonomy %>%
  group_by(family) %>%
  drop_na(family) %>% 
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>% 
  arrange(desc(n))

taxonomy %>%
  group_by(family) %>%
  drop_na(family) %>% 
  summarise(bins = n_distinct(scientificName)) %>% 
  mutate(freq = bins/sum(bins)) %>% 
  arrange(desc(bins))

# genus
taxonomy %>% drop_na(genus) %>% summarise_all(n_distinct) 

taxonomy %>%
  group_by(genus) %>%
  drop_na(genus) %>% 
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>% 
  arrange(desc(n))

taxonomy %>%
  group_by(genus) %>%
  drop_na(genus) %>% 
  summarise(bins = n_distinct(scientificName)) %>% 
  mutate(freq = bins/sum(bins)) %>% 
  arrange(desc(bins))

# species
taxonomy %>% drop_na(species) %>% summarise_all(n_distinct) 

taxonomy %>%
  group_by(species) %>%
  drop_na(species) %>% 
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)) %>% 
  arrange(desc(n))

taxonomy %>%
  group_by(species) %>%
  drop_na(species) %>% 
  summarise(bins = n_distinct(scientificName)) %>% 
  mutate(freq = bins/sum(bins)) %>% 
  arrange(desc(bins))

# how many bins are species names associated with (ordered descending)
taxonomy %>%
  drop_na(species) %>% 
  group_by(species, order) %>%
  summarise(bins = n_distinct(scientificName)) %>% 
  arrange(desc(bins))

table(taxonomy$order) # insects, hitchhikers and by-catch of non-flying orders (or prey?)

#### adding richness and diversity variables ####

pa_asvs <- decostand(asvs, method = "pa") # transform asv table into presence absence 

#Data has species as rows and sites as columns, but we need the opposite. Thus, transpose the original data
tpa_asvs <- t(pa_asvs) 
tpa_asvs[1:5,1:5] # samples as rows, asvs as columns
tasvs <- t(asvs) # read abundance transposed table
pa_asvs[1:5,1:5] # asvs are rows, samples as columns 
dtasvs <- as.data.frame(tasvs)

##### vegan ####

# get minimum number of reads for rarefying
quantile(rowSums(dtasvs), probs = seq(0, 1, 0.2)) # first see how many samples are associated with how many reads, setting it to min is quite low in this case, but at least we do not remove any more samples
min_n_seqs <- min(rowSums(dtasvs))
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
#rarecurve_data <- rarecurve(dtasvs, step = 100)

#map_dfr(rarecurve_data, bind_rows) %>% 
#  bind_cols(SampleID = rownames(dtasvs), .) %>% 
#  pivot_longer(-SampleID) %>% 
#  drop_na() %>% 
#  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
#  select(-name) %>% 
#  ggplot(aes(x=n_seqs, y = value, group = SampleID)) + 
#  geom_vline(xintercept = min_n_seqs, color = "gray") + # add min_n_seq indication
#  geom_line() + 
#  theme_classic()

##### estimated richness #####
asvs_long <- asvs %>% 
  rownames_to_column(var = "seqID") %>% 
  pivot_longer(cols = starts_with("P"), names_to = "SampleID", values_to = "value")

# testing if zeros needs to be removed
asvs_long_nozero <- asvs_long %>% filter(value > 0) 

#detach("package:seqinr", unload=TRUE)
asvs_long_nozero %>% filter(SampleID == "P100.1B") %>% 
  count(SampleID, value) %>% arrange(desc(value)) # this means, for example for the first row, that 1 asv had 20,687 reads in the sample

get_breakaway <- function(x){
  ba <- breakaway(x)
  tibble(richness_est = ba$estimate, est_richness_lci = ba$interval[1], est_richness_uci = ba$interval[2], est_richness_model = ba$model)
}

est_richness <- asvs_long_nozero %>% 
  count(SampleID, value) %>% 
  nest(data = -SampleID) %>% # creates a tibble for each sample
  mutate(n_reads = map_dbl(data, ~sum(.x$value * .x$n)), # number of reads in each sample
         obs_richness = map_dbl(data, ~sum(.x$n)), # (asv) richness per sample
         ba = purrr::map(data, ~get_breakaway(.x))) %>% # add breakaway calculations
  select(-data) %>% 
  unnest(ba)

# visualise the output
est_richness %>% 
  ggplot(aes(x = n_reads, y = richness_est, color=est_richness_model)) + geom_point() + geom_smooth() # slight increase in estimated diversity with increasing amount of reads, but some variation is clear from the different models

##### add richness and diversity to data ####
data_richness <- inner_join(data, richness, by = "SampleID") %>% 
  inner_join(., shannon, by = "SampleID") %>% inner_join(., est_richness, by = "SampleID")

# make the different richness estimate names more clear
data_richness <- data_richness %>%
  rename(richness_rarefied = richness,
         richness_rarefied_shannon = shannon)

# how different is the rarefied richness and the estimated richness (which tend to have quite narrow CIs)
data_richness %>% 
  mutate(diff = richness_rarefied - richness_est) %>% 
  summarize(mean(diff), sd(diff)) # so mean richness seems a bit lower for rarefied richness

# visualise it
data_richness %>% 
  select(SampleID, n_reads, obs_richness, richness_rarefied, richness_est) %>% 
  pivot_longer(-c(SampleID, n_reads)) %>% 
  ggplot(aes(x = n_reads, y = value, color=name)) + geom_point() + geom_smooth()

saveRDS(data_richness, file = "data/cleaned_data/data_richness_BINs.RDS") # the saved output without the BIN extension is for unique sequences 
