# standardising the sequencing data

### load libraries ######################
library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(seqinr)
library(zoo)
library(sp)
library(lubridate)

### fwh primer ##########################

# NB!!! we still need to make sure that all libraries are included - library 44 fails no matter how the analysis is carried out, but we should be able to get data from 40, 42 and 46(?) is we combine the different DADA2 runs. Also be aware that the data is from Novaseq6000 and HiSeq4000(?), and quality thresholds are much fewer for Novaseq data (4 instead of 20+) which is not handled well by the default parametres in DADA2, which we used (https://github.com/benjjneb/dada2/issues/791). The consequence appears to be that one fails to detect the rare sequences/taxa. 

#sequence data for the first sequenced libraries 1-9 (which are fwh?), 10, 13, 16, 19, 22, 25, 28, 34, 36, 38, 40, 42, 44, 46 (library 4 not included, mostly 2018 samples) 
lulufied_1 <- readRDS("data/sequencing_data/firstrun/lulified_nochim_firstrun.RDS")
otutable_1 <- lulufied_1[["curated_table"]] # extract the otutable
names(otutable_1) # contains blanks, negatives, samples from 2018 and 2019 in Denmark, and samples from 2018 in Germany

# the fasta file has been subset to only contain sequences with a length over 200 bp, so the asv table should only have asvs with the correct read length
fastas_1 <- read.fasta("data/sequencing_data/firstrun/lulufied_otus_lengthcorrected_firstrun.fasta") # notice that there are fewer fastas than asvs - this is because we subset the fastas to be >200 bp (target length is 205 bp for the fwh primer) and the short sequences are discarded
keep <- names(fastas_1)
asvs_1 <- otutable_1[(rownames(otutable_1) %in% keep), ] # this removes 2632 asvs (not the correct length)
names(asvs_1)

# rename the German sample names that were numbers and had an X added to the beginning. If we put 'S_' instead of X, then they will match the second sequencing run - and we will add the year of the sampling, since the numbering started over in 2019 so it would cause duplicate names.

names(asvs_1) <- asvs_1 %>% names() %>% str_replace("[X]", "S18_")

##### taxonomy data #####

# get taxonomy assigned with the GBIF sequence ID tool March 17th (taxonomy_1) and March 18th (the other two)
taxonomy_1 <- read.delim("data/sequencing_data/firstrun/blastresult(15).csv", sep = ",")
taxonomy_2 <- read.delim("data/sequencing_data/firstrun/blastresult(16).csv", sep = ",")
taxonomy_3 <- read.delim("data/sequencing_data/firstrun/blastresult(17).csv", sep = ",")

# merge taxonomy data
taxonomy <- rbind(taxonomy_1, taxonomy_2)
taxonomy_1 <- rbind(taxonomy, taxonomy_3)

taxa_1 <- taxonomy_1 # keep an output where taxstrings are not split up
taxonomy_1 <- taxonomy_1 %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package - there are a lot of warnings, but this is because NAs are put in where there taxonomy is not complete

#### second sequence run ####

lulufied_2 <- readRDS("data/sequencing_data/secondrun/lulified_nochim_secondrun.RDS")
otutable_2 <- lulufied_2[["curated_table"]] # extract the otutable
names(otutable_2)

# the fasta file has been subset to only contain sequences with a length over 200 bp, so the asv table should only have asvs with the correct read length
fastas_2 <- read.fasta("data/sequencing_data/secondrun/lulufied_otus_lengthcorrected_secondrun.fasta") # this file was created as a part of the diversity study
keep <- names(fastas_2)
asvs_2 <- otutable_2[(rownames(otutable_2) %in% keep), ] # 3916 asvs are removed because they were too short
names(asvs_2)

# rename the German samples so the year is included
names(asvs_2) <- asvs_2 %>% names() %>% str_replace("[S]", "S19_")

# NB!!! the main issue now is that both 2019 and 2019 German samples were sequenced in the second run - so the sample names probably need to be renamed in another way (based on sampling names or another logical way) instead of the solution above - otherwise it will be impossible to link the sequences to the correct samples 

##### taxonomy #####
taxonomy_1_1 <- read.delim("data/sequencing_data/secondrun/blastresult(7).csv", sep = ",")
taxonomy_2 <- read.delim("data/sequencing_data/secondrun/blastresult(8).csv", sep = ",")
taxonomy_3 <- read.delim("data/sequencing_data/secondrun/blastresult(9).csv", sep = ",")
taxonomy_4 <- read.delim("data/sequencing_data/secondrun/blastresult(10).csv", sep = ",")

# merge taxonomy data
taxonomy <- rbind(taxonomy_1_1, taxonomy_2)
taxonomy_2 <- rbind(taxonomy, taxonomy_3)

taxa_2 <- taxonomy_2
taxonomy_2 <- taxonomy_2 %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package

#### clean up the environment ####
rm(lulufied_1)
rm(lulufied_2)
rm(asvtable_1)
rm(asvtable_2)
rm(keep)

### minor edits to taxonomy ####

# the best assigned taxonomy needs to be specified in intrespecificEpiphet and taxonRank
# first assign the highest taxon rank
names(taxonomy)
taxonomy$taxonRank <- names(taxonomy[, c(8:14)])[max.col(!is.na(taxonomy[, c(8:14)]), "last")]

# add intraspecificEpithet  - not the most optimal solution perhaps, but it works
test <- taxonomy[, c(8:14)]
test$intraspecificEpithet <- NA

test2 <- test %>% 
  mutate(intraspecificEpithet = as.character(intraspecificEpithet)) %>%
  pmap_dfr(., ~ na.locf(c(...)) %>%
             as.list %>%
             as_tibble) %>% select(intraspecificEpithet)

taxonomy <- cbind(taxonomy, test2)