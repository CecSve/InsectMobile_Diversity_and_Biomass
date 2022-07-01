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
library(rgbif)
library(taxize)

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

# remove the German samples from the asv table since we will not include them in the analysis
#asvs_1 <- asvs_1 %>% dplyr:: select(!starts_with("S"))

##### taxonomy data #####
# Data have been matched against a 99% clustered version of the BOLD Public Database v2022-02-22 public data (COI-5P sequences) All returned matches have then been matched against the GBIF backbone taxonomy by their identifier (e.g. BOLD:ADJ8357). These OTU identifiers can be used for publishing sequence based data to GBIF. The result can be downloaded as a csv with identifiers included.

### Match types ###
#Blast exact match: identity >= 99% and queryCoverage >= 80%. This is within the threshold of the OTU.
#Blast ambiguous match:	identity >= 99% and queryCoverage >= 80%, but there is at least one more match with similar identity
#Blast close match: identity < 99% but > 90% and queryCoverage >= 80%. It is something close to the OTU, maybe the same Genus.
#Blast weak match: there is a match, but with identity < 90% or/and queryCoverage < 80%. Depending on the quality of the sequence, bit score, identity and expect value, a higher taxon could be inferred from this.
#Blast no match: No match. 

# !TAX SHOULD BE UPDATED WITH THE STAGING VERSION USING CHECKLIST ANNOTATION INSTEAD OF THE BACKBONE

# get taxonomy assigned with the GBIF sequence ID tool May 16th 2022 
taxonomy_1 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool_nonBackbone/blastresult(27).csv", sep = ",") # 5,000 sequences, 99.2%with blast match, 85%, with identity > 99%, 92% with identity >95%
taxonomy_2 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool_nonBackbone/blastresult(28).csv", sep = ",") # 5,000 sequences, 97% with blast match, 73% with identity > 99%, 84% with identity >95%
taxonomy_3 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool_nonBackbone/blastresult(29).csv", sep = ",") #  5,000 sequences, 91% with blast match, 54% with identity > 99%, 66% with identity >95%
taxonomy_4 <- read.delim("data/sequencing_data/firstrun/GBIF_seq_id_tool_nonBackbone/blastresult(30).csv", sep = ",") # #  2,485 sequences, 45% with blast match, 19% with identity > 99%, 24% with identity >95%

# merge taxonomy data
taxonomy <- rbind(taxonomy_1, taxonomy_2)
taxonomy_1 <- rbind(taxonomy, taxonomy_3)
taxonomy_1 <- rbind(taxonomy_1, taxonomy_4)

taxa_1 <- taxonomy_1 # keep an output where taxstrings are not split up
taxonomy_1 <- taxonomy_1 %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package - there are a lot of warnings, but this is because NAs are put in where there taxonomy is not complete

#### second sequence run ####
# fwh libraries: 58 (DE samples 1-94), 59 (DE samples 95-188), 60 (+ DE samples 189-215), 62, 64, 65, 66, 67, 68, 69, 7, 71, 72, 73
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
names(check_de) # 206 samples with some blanks and negatives

# option to remove German samples
#asvs_2 <- asvs_2 %>% dplyr:: select(!starts_with("S"))

##### taxonomy #####
taxonomy_1_1 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool_nonBackbone/blastresult(31).csv", sep = ",") #   5,000 sequences, 99.9% with blast match, 85% with identity > 99%, 93% with identity >95%

taxonomy_2 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool_nonBackbone/blastresult(32).csv", sep = ",") #   5,000 sequences, 99.3% with blast match, 73% with identity > 99%, 86% with identity >95%

taxonomy_3 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool_nonBackbone/blastresult(33).csv", sep = ",") #  5,000 sequences 94% with blast match, 53% with identity > 99%, 66% with identity >95%

taxonomy_4 <- read.delim("data/sequencing_data/secondrun/GBIF_seq_id_tool_nonBackbone/blastresult(34).csv", sep = ",") #  1,121 sequences, 81% with blast match, 32% with identity > 99%, 43% with identity >95%

# merge taxonomy data
taxonomy <- rbind(taxonomy_1_1, taxonomy_2)
taxonomy_2 <- rbind(taxonomy, taxonomy_3)
taxonomy_2 <- rbind(taxonomy_2, taxonomy_4)

taxa_2 <- taxonomy_2
taxonomy_2 <- taxonomy_2 %>% separate(classification, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "_") #split the taxonomy string into ranks using the dplyr and tidyr package

#### clean up the environment ####
rm(lulufied_1)
rm(lulufied_2)
rm(keep)

### merge the two sequence runs ####
asvs_1 <- asvs_1 %>% rownames_to_column(var = "otuid") 
asvs_2 <- asvs_2 %>% rownames_to_column(var = "otuid") 

asvs <- merge(asvs_1, asvs_2, by = "otuid", all = TRUE) # combine the two asvtables and keep all the samples while merging identical asv IDs

asvs <- asvs %>% column_to_rownames(var = "otuid")
asvs[is.na(asvs)] <- 0 # replace na's with zeroes

# check sample size
min(colSums(asvs))
mean(colSums(asvs))
median(colSums(asvs))
max(colSums(asvs))

#### merge the taxonomy ####
taxonomy <- merge(taxonomy_1, taxonomy_2, all = TRUE) # combine the two asvtables and keep all the samples while merging identical asv IDs 

# notice that taxonomy obs. and asvs obs. match in number of obs.

### name parsing ####
# since we use the BOLD checklist for assigning names, we get a lot of names that are not true Linnean names, e.g. OTU names and sp. placeholders. They need to be identified and the taxonomy should be reassigned to the most correct name

# remove non-Linnean names (INFORMAL (a scientific name with some informal addition like "cf." or indetermined like Abies spec.) and OTU names) from genus rank names
genus_parsed <- taxonomy %>%
  mutate(parsed = rgbif::parsenames(genus)) %>%
  tidyr::unnest(cols = "parsed") %>% mutate_at(vars(genus), ~ replace(., type == "INFORMAL", NA)) %>% mutate_at(vars(genus), ~ replace(., type == "OTU", NA)) %>% mutate_at(vars(genus), ~ replace(., type == "NO_NAME", NA))

# for some reason the Janzen and Malaise names make it through the parser, although I don't know why - they need to be removed before we continue

genus_parsed <- genus_parsed %>%
  mutate(across(c(species),
                ~ replace(.,   str_detect(., "Janzen|Malaise"), NA)))

# the Janzens in the scientificName will be removed in the next step

# remove the parse columns so we can parse species after genus is parsed
genus_noparse <- genus_parsed[, c(1:15)]

genus_species_parsed <-
  genus_noparse %>% mutate(parsed = rgbif::parsenames(species)) %>% tidyr::unnest(cols = "parsed") %>% mutate_at(vars(species), ~ replace(., type == "INFORMAL", NA)) %>% mutate_at(vars(species), ~ replace(., type == "OTU", NA)) %>% mutate_at(vars(genus), ~ replace(., type == "NO_NAME", NA))

# summarize the difference
table(is.na(taxonomy$species)) # this many species had no name before cleaning
table(is.na(genus_species_parsed$species)) # this many species do not have a name after cleaning

#### remove unwanted IDs ####
# NB! I noticed wild boar in the sequences - remove everything but spiders and insects

unique(genus_species_parsed$class)
genus_species_parsed <-
  genus_species_parsed[genus_species_parsed$class %in% c('Arachnida', 'Insecta'),] # 22990 records

### minor edits to taxonomy ####

# the best assigned taxonomy needs to be specified in infraspecificEpithet  and taxonRank

# first assign the highest taxon rank
names(genus_species_parsed)
genus_species_parsed$taxonRank <-
  names(genus_species_parsed[, c(8:14)])[max.col(!is.na(genus_species_parsed[, c(8:14)]), "last")]

# add intraspecificEpithet  - not the most optimal solution perhaps, but it works
test <- genus_species_parsed[, c(8:14)]
test$infraspecificEpithet <- NA

epithet <- test %>% 
  mutate(infraspecificEpithet = as.character(infraspecificEpithet)) %>%
  pmap_dfr(., ~ na.locf(c(...)) %>%
             as.list %>%
             as_tibble) %>% select(infraspecificEpithet)

genus_species_parsed <- genus_species_parsed %>% select(!infraspecificepithet)

taxonomy_cleaned <- cbind(genus_species_parsed, epithet) # how to deal with the 'no match' sequences? Should Biota be added as the domain (taxonRank = kingdom)? This is the hacky solution that is possible in GBIF right now

# current output
#write.table(taxonomy_cleaned, file= "data/sequencing_data/taxonomy_cleaned.txt", sep="\t", col.names = T, row.names = F)

# save output
#saveRDS(taxonomy_cleaned, file = "data/sequencing_data/taxonomy_cleaned.rds")
#saveRDS(asvs, file = "data/sequencing_data/asvs.rds")
#write.table(asvs, file = "data/sequencing_data/asvtable.txt", col.names = NA, sep = "\t")
#write.table(taxonomy, file = "data/sequencing_data/taxonomy.txt", col.names = T, row.names = F, sep = "\t")


# MATCH NAMES TO CURRENT TAXONOMIC STATUS TO MAKE SURE WE HAVE THE ACCPTED NAMES (AND KNOW THE SYNONYMS) - 
# THERE MAY BE A MANUEL CHECK INCLUDED FOR THE NAMES THAT DO NOT MATCH 
# we will use taxize for this, Diana will carry out the check

#get all species names
allSpecies <- taxonomy_cleaned %>%
  filter(!is.na(species)) %>%
  pull(species) %>%
  unique() 

#check with gbif parse
checkSpecies <- allSpecies %>% map_dfr(gbif_parse) %>%
                  add_column(originalName = allSpecies)

mean(checkSpecies$parsed==TRUE) # proportion that pass

#which ones don't parse
#missingSpecies <- checkSpecies %>% filter(parsed==FALSE)

#bold_search(missingSpecies$scientificname) # this should be fixed in the section above

#rgbif - check names against catalogue of life

#check species
checkSpecies <- lapply(allSpecies, function(x){
  name_backbone(name = x, rank="species")
}) %>%
  reduce(bind_rows)

#how many synonymsn do we have?
table(checkSpecies$status)

synonym <- checkSpecies %>%
            filter(status!="ACCEPTED")

# save the synonyms so we can check if needed
#write.table(synonym, file= "data/sequencing_data/synonyms.txt", sep="\t", col.names = T, row.names = F)

#the species column contains the accepted names for these species
