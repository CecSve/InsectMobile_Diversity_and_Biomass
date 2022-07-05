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

# test <- sampling_data %>% filter(!is.na(Date)) # only retain samples with an associated date

# NB GOT TO HERE - NEED TO MERGE SAMPLES IN THE ASV TABLE - this means that the code below is just notes, and needs to be finalised

### merging total samples from size sorted samples ###########
# copy the sampleID_size to the NA SampleID cells
lab_data <- lab_data %>% 
  mutate(SampleID = coalesce(SampleID_size,SampleID))

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

# prepare the asv table by making a column with PCRIDs
t.asvs <- t(asv_table)
t.asvs <- as.data.frame(t.asvs) %>% rownames_to_column(var = "PCRID") 
test <- full_join(keep, t.asvs, by = "PCRID") # this step removes the German samples
test[, 3:26344][is.na(test[, 3:26344])] <- 0 # replace introduced NAs with zero

# sumarize the reads per total sample instead of by size fraction
test2 <- test %>% select(-PCRID) %>% group_by(SampleID) %>% summarise_all(list(sum))
# see if any sampleIDs are missing - when we carry it out on test instead we see it is samples from 2017 and bird feces samples and they should be excluded from the analysis anyways
check <- test2 %>% filter(is.na(SampleID))

# remove the NA sample
asvs_combined <- test2 %>% filter(!is.na(SampleID)) 

totsample_asvs <- asvs_combined %>% column_to_rownames(var = "SampleID")
totsample_asvs <- as.data.frame(t(totsample_asvs))

#### check-ups ####
# get summaries of how many samples there is for each variable and their levels
test <- anti_join(sampling_data, lab_data, by = "SampleID") # sampling data but not lab data, not sampled entries and pilots that did not sample after the instructions 
test2 <- anti_join(lab_data, sampling_data, by = "SampleID") # in lab data but not sampling data, mostly blanks, negative, ethanol test and other lab test subsamples 
length(unique(sampling_lab_data_filter[["RouteID_JB"]])) # how many routes  
length(unique(sampling_lab_data_filter[["PID"]])) # how many pilots
data.frame(table(sampling_lab_data_filter$Wind)) # how often were the different wind categories registered  - the empty entries is routes where we either do not have metadata or that they weren't sampled
data.frame(table(sampling_lab_data_filter$Temperature)) # how many samples were collected at different temperature intervals

# how many collected samples?
investigate <- sampling_data %>%
  filter(str_detect(SampleID, "P")) %>% filter(!is.na(Date)) 

# how many lab processed samples?
investigate <- lab_data %>%
  filter(str_detect(SampleID, "P")) %>% filter(!is.na(DryMass_mg) & !is.na(PCRID)) #& !DryMass_mg < 0) # techically, the negative biomass samples were processed but the should be removed as they do not make any sense

#### plotting sampling effort by year ####
test <- data.frame(table(sampling_lab_data_filter$Date)) # how many samples per day

# noticed two samples from May 20th 2019! 
test <- test %>% filter(!Var1 == "2019-05-20")

#Format dates
df <- test
df$Date <- as.Date(test$Var1)
#Create year
df$Year <- format(df$Date,'%Y')
#Create day and month of year
df$Day <- format(df$Date,'%d')
df$Month <- format(df$Date,'%m')
#Assign a dummy date
df$DayMonth <- as.Date(paste0(2019,'-',df$Month,'-',df$Day))
#Now sketch for plot
ggplot(data = df, aes(x = DayMonth, y = Freq, color = Year, group = Year)) + 
  geom_line(size = 2) +
  scale_x_date(date_breaks = "3 days", date_labels = "%d-%m") + theme_bw() + scale_color_manual(values=c("#CC6666", "#9999CC")) +
  xlab("sampling date") + ylab("# of samples collected")

#### match asvtable data with lab meta data (not all samples were sequenced) ####
keep <- colnames(asv_table)
dplyr::setdiff(lab_data$PCRID, colnames(asv_table))
missingsamples <- anti_join(lab_data, t.asvs) # samples in the lab data that is not in the sequencing output

# how many of the samples that should have been in the sequencing data are true samples and how are they arranged (will help us manually sort out which libraries are affected)
investigate <- missingsamples %>%
  filter(str_detect(SampleID, "P")) %>% arrange(desc(PCRID))

# library 36, 44 is missing from the 2019 data (we expected an issue with 44, not 36). the rest seems to be few samples here and there and un-sequenced samples

### got to here! #### data should be aligned for analysis from here

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