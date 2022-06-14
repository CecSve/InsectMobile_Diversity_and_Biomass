# Standardising lab data

### lab data from 2018 ####
sampleid_size_drymass <- read.delim("data/lab_data/sampleid_size_drymass.txt")
sampleid_size_pcrid_concentration <- read.delim("data/lab_data/sampleid_size_pcrid_concentration.txt")
#sampleid_small_large_biomass <- read.delim("2018/sampleid_small_large_biomass.txt") # not used
imagerecognition <- read.delim("data/lab_data/imagerecognition.txt")

# merge imagerecognition and sampleid_size_drymass
labdata1 <- sampleid_size_drymass %>% filter(!DryMass_mg == 0)
labdata2 <- imagerecognition %>% filter(!DryMass_mg == 0)

mergedData <- rbind(labdata1, labdata2)
duplicated(mergedData$SampleID_size) # check for duplicates

mergedData <- merge(mergedData, sampleid_size_pcrid_concentration, by = "SampleID_size", all = T) # add qubit measurement - NB concentration is DNA concentration of the purified DNA extract

### lab data from 2019 ####
# for the biomass data, lab data was sorted by 1) size fractions, 2) whether biomass was measured, 3) if different weights or only ZM weights were used, then the data was flagged for potential erronous biomass estimates (biomassuncertainty = high)

sampleid_small_drymass <- read.delim("data/lab_data/sampleid_size_drymass_small_2019.txt")
sampleid_large_drymass <- read.delim("data/lab_data/sampleid_size_drymass_large_2019.txt")
pcrid_qubit_19 <- read.delim("data/lab_data/sampleid_size_pcrid_qubit_2019.txt")

biomass19 <- rbind(sampleid_small_drymass, sampleid_large_drymass)
biomass19 <- biomass19 %>% select(SampleID_size, SampleID, DryMass_mg, biomassUncertainty)
emptylines <- drop_na(biomass19) # empty entries appeared upon merge and they should be deleted

str(sampleid_size_drymass)
str(emptylines)

labdata2019 <- merge(emptylines, pcrid_qubit_19, by = "SampleID_size", all = T) # add qubit measurement - NB concentration is DNA concentration of the purified DNA extract

labdata2019 <- labdata2019[!(labdata2019$SampleID_size == ""), ] # remove empty rows
labdata2019$SampleID <- ifelse(is.na(labdata2019$SampleID), labdata2019$SampleID_size, labdata2019$SampleID) # copy extraction blank sample IDs to the SampleID column

#### merge & format data ###########
mergedData$biomassUncertainty <- ""

# check if columns match
str(mergedData)
str(labdata2019)
combdata <- rbind(mergedData, labdata2019) # merging lab processed data from 2018 and 2019

# standardise qubit measurement values by changing too low to zero and too high readings to 100
test <- combdata %>% 
  mutate(concentration = ifelse(as.character(concentration) == "too low", "0", as.character(concentration))) 

labdata <- test %>% 
  mutate(concentration = ifelse(as.character(concentration) == "Too low", "0", as.character(concentration)))

labdata <- labdata %>% 
  mutate(concentration = ifelse(as.character(concentration) == "Too Low", "0", as.character(concentration)))

labdata <- labdata %>% 
  mutate(concentration = ifelse(as.character(concentration) == "too high", "100", as.character(concentration)))

# save output
write.table(labdata, file = "data/lab_data/lab_data_cleaned.txt", sep = "\t", row.names = F)
