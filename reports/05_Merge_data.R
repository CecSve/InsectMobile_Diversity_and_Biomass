# Merge data tables prior to analysis - this script should generate the final table(s) for analysis and the output to be shared in a data sharing repository (journal specific)

library(tidyverse) # data wrangling
library(readr) # import data
library(lubridate) # working with dates

### read in cleaned data tables ############

# sampling metadata
sampling_data <- readRDS("data/sampling_data/sampling_data_cleaned.rds")

# samples processed in the lab
lab_data <- readRDS("data/lab_data/lab_data_cleaned.rds")

# sequencing data
asv_table <- readRDS("data/sequencing_data/asvs.rds")

# taxonomy data - NB! prior to taxize edits!
taxonomy_data <- readRDS("data/sequencing_data/taxonomy_cleaned.rds") # less records than the asv table but this is because the taxonomy has been filtered to only arachnids and insects!

# environmental data
environData <- read_delim("data/environmental_data/covariate-data/DK_environData.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

# test <- sampling_data %>% filter(!is.na(Date)) # only retain samples with an associated date

### merging total samples from size sorted samples ###########
# copy the sampleID_size to the NA SampleID cells
lab_data <- lab_data %>% 
  mutate(SampleID = coalesce(SampleID_size, SampleID))

# remove size fraction delimiters from the sampleID
lab_data$SampleID<-gsub("S","",as.character(lab_data$SampleID))
lab_data$SampleID<-gsub("L","",as.character(lab_data$SampleID)) # there are still some ethanol samples in there that eventually should be removed, called for example P1.1AE

#### merge sampling data and lab data ####

# merge and keep blanks and negatives in
sampling_lab_data <- full_join(sampling_data, lab_data, by = "SampleID")

# retain unique samples with size fractioning 
data_unique <- sampling_lab_data %>% distinct(SampleID_size, .keep_all = TRUE)
examine <- setdiff(sampling_lab_data, data_unique) # seems to be removing the unsampled routes mostly, but not only. Two sampled routes have duplicate samples in the sampling_lab_data sheet for some reason.The issue with P27 is easy to fix:

sampling_lab_data_filter <- sampling_lab_data %>%
  arrange(SampleID_size, -DryMass_mg) %>%
  filter(!duplicated(SampleID_size))

examine <- setdiff(sampling_lab_data, sampling_lab_data_filter) # now only duplicates of the P27.1B and P112.1BS are removed (and all the unsampled routes)

keep <- sampling_lab_data_filter %>% dplyr::select(PCRID, SampleID) # we need to use the filter to avoid duplicate PCRIDs

# create a column for total biomass and mean DNA concentration
# first we need numeric concentrations
sampling_lab_data_filter$concentration = as.numeric(as.character(sampling_lab_data_filter$concentration))

sampling_lab_data_filter <- sampling_lab_data_filter %>% group_by(SampleID) %>% mutate(totalBiomass_mg = sum(DryMass_mg), meanDNAconc = mean(concentration, na.rm = TRUE)) %>% ungroup()

# prepare the asv table by making a column with PCRIDs
t.asvs <- t(asv_table)
t.asvs <- as.data.frame(t.asvs) %>% rownames_to_column(var = "PCRID") 
test <- full_join(keep, t.asvs, by = "PCRID") # merge sampleID, PCRID with the sequenced samples - this step removes the German samples
str(test)
test[, 3:26344][is.na(test[, 3:26344])] <- 0 # replace introduced NAs with zero

# summarize the reads per total sample instead of by size fraction
test2 <- test %>% select(-PCRID) %>% group_by(SampleID) %>% summarise_all(list(sum))

# see if any sample IDs are missing - when we carry it out on test instead we see it is samples from 2017 and bird feces samples and they should be excluded from the analysis anyways
check <- test2 %>% filter(is.na(SampleID))

# remove the NA sample
asvs_combined <- test2 %>% filter(!is.na(SampleID)) 
str(asvs_combined)
totsample_asvs <- asvs_combined %>% column_to_rownames(var = "SampleID") 
str(totsample_asvs)
# revert back to asv table format with samples as columns and asvs as rows
totsample_asvs <- as.data.frame(t(totsample_asvs))

# save output 
# saveRDS(totsample_asvs, file = "data/sequencing_data/asvtable_total_samples.rds") 

#### check-ups ####
# get summaries of how many samples there is for each variable and their levels
test <- anti_join(sampling_data, lab_data, by = "SampleID") # in sampling data but not in lab data, not sampled entries and pilots that did not sample after the instructions 

test2 <- anti_join(lab_data, sampling_data, by = "SampleID") # in lab data but not in sampling data, mostly blanks, negative, ethanol test and other lab test sub samples 

# how many collected samples?
investigate <- sampling_data %>%
  filter(str_detect(SampleID, "P")) %>% filter(!is.na(Date))

# how many lab processed samples?
investigate <- lab_data %>%
  filter(str_detect(SampleID, "P")) %>% filter(!is.na(DryMass_mg) & !is.na(PCRID) & !DryMass_mg < 0 & !biomassUncertainty == "high") %>% select(-biomassUncertainty) # technically, the negative biomass samples were processed but the should be removed as they do not make any sense incl. the biomass measured with high uncertainty

#### match asvtable data with lab meta data (not all samples were sequenced) ####
#keep <- colnames(asv_table)
#dplyr::setdiff(lab_data$PCRID, colnames(asv_table))
#missingsamples <- anti_join(lab_data, t.asvs) # samples in the lab data that is not in the sequencing output

# how many of the samples that should have been in the sequencing data are true samples and how are they arranged (will help us manually sort out which libraries are affected)
#investigate <- missingsamples %>% filter(str_detect(SampleID, "P")) %>% arrange(desc(PCRID))

# library 36, 44 is missing from the 2019 data (we expected an issue with 44, not 36). the rest seems to be few samples here and there and un-sequenced samples

##### combine sampling data, lab data (one table already) with environmental data #######

# to merge, we need a Year column in the sampling lab data
sampling_lab_data_filter$Year <- year(ymd(sampling_lab_data_filter$Date))

# merge tables
combData <- merge(sampling_lab_data_filter, environData, by.x = c("RouteID_JB", "Year"), by.y = c("routeID", "Year"), all = T) # still contains blanks, sample size fractions

# remove NA samples (samples that have environmental variables calculated but no data in sampling and lab tables)
combData <- combData %>% filter(!is.na(SampleID)) %>% filter(!is.na(RouteID_JB)) 
combData$concentrationUnit <- combData$concentrationUnit %>% replace_na('ng/µl')

# remove size fraction and data associated with size fractions and keep only unique rows/sampleIDs
combData_noSizeFrac <- combData %>% select(-SampleID_size, -DryMass_mg, -PCRID, -concentration) %>% distinct(.keep_all = TRUE)

### check samples between sequenced samples and samples from the combined data ####

keep <- colnames(totsample_asvs) # the sample IDs of the samples we have sequenced (incl. blanks and negatives)

check <- setdiff(colnames(totsample_asvs), combData_noSizeFrac$SampleID) # in asv table but not in the other data - blanks and tests
setdiff(combData_noSizeFrac$SampleID, colnames(totsample_asvs)) # the other way around, no match

# check read count of blanks before removing them
test <- totsample_asvs[, which((names(totsample_asvs) %in% check)==TRUE)]
colSums(test) # blanks are empty, tests are not - it is ok to remove them!

keep <- combData_noSizeFrac$SampleID

# remove the blanks from the asv table
asvtable <- totsample_asvs[ ,colnames(totsample_asvs) %in% keep]
rowSums(asvtable)

# check whether any samples are missing
setdiff(combData_noSizeFrac$SampleID, colnames(asvtable))
setdiff(colnames(asvtable), combData_noSizeFrac$SampleID)

asvs <-
  asvtable[apply(asvtable[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (ASVs that are not present in the subsetted samples)

keep <- rownames(asvs)

# remove taxonomy data for sequences that are not present in the subsetted asv table
taxonomy <- taxonomy_data %>% filter(occurrenceId %in% keep) 

# keep only asvs that passed the taxonomic assignment and filtering
20874-18707 # this step removes 2167 sequences with no taxonomic assignment to insecta and arachnids

# some samples do not have any reads after filtering and cleaning
asvs_subset = asvs[, !(colSums(asvs) == 0)]
1324-1180 # 144 samples do not have sequences - these were processed in lab, but never send for sequencing
asvs_examine = asvs[, (colSums(asvs) == 0)]

keep <- taxonomy$occurrenceId
asv <- asvs_subset[rownames(asvs_subset) %in% keep, ] # now the asv table and the taxonomy table match

# match samples between data and asvs
keep <- colnames(asvs)
combData_noSizeFrac <- combData_noSizeFrac %>% filter(SampleID %in% keep) 

# save output for analysis - now all tables are aligned
saveRDS(taxonomy, file = "data/cleaned_data/taxonomy_filtered.rds")
saveRDS(asv, file = "data/cleaned_data/asvs_filtered.rds") 
saveRDS(combData_noSizeFrac, file = "data/cleaned_data/combData_noSizeFrac.RDS") 
