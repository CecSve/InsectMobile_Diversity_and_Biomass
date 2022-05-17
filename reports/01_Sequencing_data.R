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
library(data.table)

### fwh primer ##########################

# NB!!! Library 44 fails no matter how the analysis is carried out, but data from 40, 42 and 46 is included because we combine the different DADA2 runs (the three libraries also fail in some combinations of library processing with DADA2 and it is not due to tag overlap). Also be aware that the data is from Novaseq6000 and HiSeq4000, and quality thresholds are very different for Novaseq data (4 instead of 20+ fopr HiSeq) which is not handled well by the default parameters in DADA2, which we used (https://github.com/benjjneb/dada2/issues/791). The consequence appears to be that we may fail to detect the rare sequences/taxa with the NovaSeq data. 

# sequence data for the first sequenced libraries 1, 4, 7 (German samples - MISSING!), 10, 13, 16, 19, 22, 25, 28 (contains German samples), 31 (German samples), 34, 36, 38, 40, 42 (contains German samples), (44 - fails in the DADA2 pipeline so not included), 46  
lulufied_1 <- readRDS("data/sequencing_data/firstrun/lulified_nochim_rerun2022.RDS") # updated bioinformatics
otutable_1 <- lulufied_1[["curated_table"]] # extract the otutable
names(otutable_1) # contains blanks, negatives, samples from 2017 (IM17_*), 2018 (IM18*_) and 2019 (IM19_*) in Denmark, and samples from 2018 (X*) in Germany

# the fasta file has been subset to only contain sequences with a length over 200 bp, so the asv table should only have asvs with the correct read length
fastas_1 <- read.fasta("data/sequencing_data/firstrun/lulufied_otus_lengthcorrected_rerun2022.fasta") # notice that there are fewer fastas than asvs - this is because we subset the fastas to be >200 bp (target length is 205 bp for the fwh primer) and the short sequences are discarded
keep <- names(fastas_1)
asvs_1 <- otutable_1[(rownames(otutable_1) %in% keep), ] # this removes 3177 asvs (not the correct length - note to future self: could any of these be true Hymenoptera reads following Zizkas comment?)
names(asvs_1)

# check-up on sequences in samples
min(colSums(asvs_1)) # 0
mean(colSums(asvs_1)) # 91013.17
median(colSums(asvs_1)) # 75243.5
max(colSums(asvs_1)) # 444566

# rename the German sample names that were numbers and had an X added to the beginning. If we put 'S_' instead of X, then they will match the second sequencing run - and we will add the year of the sampling, since the numbering started over in 2019 so it would cause duplicate names.

names(asvs_1) <- asvs_1 %>% names() %>% str_replace("[X]", "S18_")
check_de <- asvs_1 %>% dplyr:: select(starts_with("S"))
names(check_de) # 237 samples out of 252 present (some negatives in there but most are discarded)

##### taxonomy data #####
# Data have been matched against a 99% clustered version of the BOLD Public Database v2022-02-22 public data (COI-5P sequences) All returned matches have then been matched against the GBIF backbone taxonomy by their identifier (e.g. BOLD:ADJ8357). These OTU identifiers can be used for publishing sequence based data to GBIF. The result can be downloaded as a csv with identifiers included.

### Match types ###
#Blast exact match: identity >= 99% and queryCoverage >= 80%. This is within the threshold of the OTU.
#Blast ambiguous match:	identity >= 99% and queryCoverage >= 80%, but there is at least one more match with similar identity
#Blast close match: identity < 99% but > 90% and queryCoverage >= 80%. It is something close to the OTU, maybe the same Genus.
#Blast weak match: there is a match, but with identity < 90% or/and queryCoverage < 80%. Depending on the quality of the sequence, bit score, identity and expect value, a higher taxon could be inferred from this.
#Blast no match: No match. 

# get taxonomy assigned with the GBIF sequence ID tool May 16th 2022 
taxonomy_1 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool/blastresult.csv", sep = ",") # 5,000 sequences, 99.2%with blast match, 85%, with identity > 99%, 96% with GBIF backbone match
taxonomy_2 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool/blastresult(1).csv", sep = ",") # 5,000 sequences, 97% with blast match, 73% with identity > 99%, 92% with GBIF backbone match
taxonomy_3 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool/blastresult(2).csv", sep = ",") #  5,000 sequences, 91% with blast match, 54% with identity > 99%, 84% with GBIF backbone match
taxonomy_4 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool/blastresult(3).csv", sep = ",") # #  2,485 sequences, 45% with blast match, 19% with identity > 99%, 42% with GBIF backbone match

# merge taxonomy data
taxonomy <- rbind(taxonomy_1, taxonomy_2)
taxonomy_1 <- rbind(taxonomy, taxonomy_3)
taxonomy_1 <- rbind(taxonomy_1, taxonomy_4)

taxa_1 <- taxonomy_1 # keep an output where taxstrings are not split up
taxonomy_1 <- taxonomy_1 %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package - there are a lot of warnings, but this is because NAs are put in where there taxonomy is not complete

# split fasta into 50 sequences
#subtax <- taxonomy_1[taxonomy_1$identity == '100', ] # only choose the sequences with a clear match
#new_DF <- subtax[is.na(subtax$species),] # only choose sequences not identified to species level
new_newDF <- subtax[is.na(subtax$kingdom),]

#keep <- new_DF$occurrenceId
#subfas <- fastas_1[(names(fastas_1) %in% keep)] # only keep those fastas that where defined above

#chunk <- function(x, n) (mapply(function(a, b) (x[a:b]), seq.int(from=1, to=length(x), by=n), pmin(seq.int(from=1, to=length(x), by=n)+(n-1), length(x)), SIMPLIFY=FALSE))

#fassplit <- chunk(subfas, 50) # output not saved in individual files yet - could be used for manual blast against BOLD (only 50 sequences per turn), to get the most updated taxonomy. But this would require 46 individual blast sessions just for the first sequencing run!

#### second sequence run ####
lulufied_2 <- readRDS("data/sequencing_data/secondrun/lulified_nochim_secondrun.RDS")
otutable_2 <- lulufied_2[["curated_table"]] # extract the otutable
names(otutable_2)

# the fasta file has been subset to only contain sequences with a length over 200 bp, so the asv table should only have asvs with the correct read length
fastas_2 <- read.fasta("data/sequencing_data/secondrun/lulufied_otus_lengthcorrected_secondrun.fasta") # this file was created as a part of the diversity study
keep <- names(fastas_2)
asvs_2 <- otutable_2[(rownames(otutable_2) %in% keep), ] # 3916 asvs are removed because they were too short
names(asvs_2)

# check-up on sequences in samples - generally much higher numbers than the first sequencing run - this is probably due to that most libraries were run on NovaSeq
min(colSums(asvs_2)) # 0
mean(colSums(asvs_2)) # 142182.1
median(colSums(asvs_2)) # 123150
max(colSums(asvs_2)) # 770966

# rename the German samples so the year is included
names(asvs_2) <- asvs_2 %>% names() %>% str_replace("[S]", "S19")
check_de <- asvs_2 %>% dplyr:: select(starts_with("S"))
names(check_de)

# NB!!! the main issue now is that both 2019 and 2019 German samples MAYBE were sequenced in the second run - so the sample names probably need to be renamed in another way (based on sampling names or another logical way) instead of the solution above - otherwise it will be impossible to link the sequences to the correct samples 

##### taxonomy #####

taxonomy_1_1 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool/blastresult(23).csv", sep = ",") #   5,000 sequences, 99.9% with blast match, 85% with identity > 99%, 97% with GBIF backbone match

taxonomy_2 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool/blastresult(24).csv", sep = ",") #   5,000 sequences, 99.3% with blast match, 73% with identity > 99%, 94% with GBIF backbone match

taxonomy_3 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool/blastresult(25).csv", sep = ",") #  5,000 sequences 94% with blast match, 53% with identity > 99%, 87% with GBIF backbone match

taxonomy_4 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool/blastresult(26).csv", sep = ",") #  1,121 sequences, 81% with blast match, 32% with identity > 99%, 74% with GBIF backbone match

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