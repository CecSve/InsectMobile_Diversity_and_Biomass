# Merge data tables prior to analysis - this script should generate the final table(s) for analysis and the output to be shared in a data sharing repository (journal specific)

### read in cleaned data tables ############

# sampling metadata
sampling_data <- readRDS("data/sampling_data/sampling_data_cleaned.rds")

# samples processed in the lab
lab_data <- readRDS("data/lab_data/lab_data_cleaned.rds")

# sequencing data
asv_table <- readRDS("data/sequencing_data/asvs.rds")

# NB! prior to taxize edits!
taxonomy_data <- readRDS("data/sequencing_data/taxonomy_cleaned.rds") # less records than the asv table but this is because the taxonomy has been filtered to only arachnids and insects!

test <- sampling_data %>% filter(!is.na(Date)) # only retain samples with an associated date

# get summaries of how many samples there is for each variable and their levels
length(unique(sampling_data[["RouteID_JB"]])) # how many routes  
length(unique(sampling_data[["PID"]])) # how many pilots
length(unique(sampling_data[["SampleID"]])) # how many samples
data.frame(table(sampling_data$Wind)) # how often were the different wind categories registered  - the empty entries is routes where we either do not have metadata or that they weren't sampled
data.frame(table(sampling_data$Temperature)) # how many samples were collected at different temperature intervals
data.frame(table(sampling_data$eventDate)) # how many samples per day

# NB GOT TO HERE - NEED TO MERGE SAMPLES IN THE ASV TABLE - this means that the code below is just notes, and needs to be finalised

### merging total samples from size sorted samples ###########
# 
data_unique <- data %>% distinct(SampleID, .keep_all = TRUE)
keep <- data %>% dplyr::select(PCRID, SampleID)

t.asvs <- t(asvs)
t.asvs <- as.data.frame(t.asvs) %>% rownames_to_column(var = "PCRID") 
test <- full_join(keep, t.asvs, by = "PCRID")

test2 <- test %>% select(-PCRID) %>% group_by(SampleID) %>% summarise_all(list(sum))

totsample_asvs <- test2 %>% column_to_rownames(var = "SampleID")
totsample_asvs <- as.data.frame(t(totsample_asvs))

#### match asvtable data with lab meta data (not all samples were sequenced) ####
keep <- colnames(asvs)
dplyr::setdiff(data$PCRID, colnames(asvs))
tasvs <- as.data.frame(t(asvs))
tasvs <- tasvs %>% rownames_to_column(var = "PCRID")
missingsamples <- anti_join(data, tasvs) # Rows in the data that do not have a match in the asv table. 266 samples are not in the ASV table and at least some of them were sequenced but apparently they were removed during the bioinformatics. It seems library 38 and 44 is largely missing - check this.
labdata <- data %>% filter(PCRID %in% keep)
nomatchsamples <- anti_join(tasvs, labdata) # rows in the asv table that don't have a match in the labdata - many bird droppings samples and some samples from lib 28
keep <- labdata$PCRID

asvtable <- asvs[ ,colnames(asvs) %in% keep]
rowSums(asvtable)

otus <-
  asvtable[apply(asvtable[,-1], 1, function(x)
    ! all(x == 0)),] # remove rows that contain only zeros (OTUs that are not present in the subsetted samples)

keep <- rownames(otus)

taxonomy <- taxonomy %>% filter(occurrenceId %in% keep)
min(colSums(otus))