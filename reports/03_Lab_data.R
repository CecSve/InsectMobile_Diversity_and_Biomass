# Standardising lab data

### lab data from 2018 ####
sampleid_size_drymass <- read.delim("data/lab_data/sampleid_size_drymass.txt")
sampleid_size_pcrid_concentration <- read.delim("data/lab_data/sampleid_size_pcrid_concentration.txt")
#sampleid_small_large_biomass <- read.delim("2018/sampleid_small_large_biomass.txt") # not used
imagerecognition <- read.delim("data/lab_data/imagerecognition.txt")

# merge imagerecognition and sampleid_size_drymass
test <- sampleid_size_drymass %>% filter(!DryMass_mg == 0)
test2 <- imagerecognition %>% filter(!DryMass_mg == 0)

mergedData <- rbind(test, test2)
duplicated(mergedData$SampleID_size) 

### lab data from 2019 ####
# for the biomass data, lab data was sorted by 1) size fractions, 2) whether biomass was measured, 3) if different weights or only ZM weights were used, then the data was flagged for potential erronous biomass estimates (biomassuncertainty = high)

sampleid_small_drymass <- read.delim("data/lab_data/sampleid_size_drymass_small_2019.txt")
sampleid_large_drymass <- read.delim("data/lab_data/sampleid_size_drymass_large_2019.txt")
pcrid_qubit_19 <- read.delim("data/lab_data/sampleid_size_pcrid_qubit_2019.txt")

biomass19 <- rbind(sampleid_small_drymass, sampleid_large_drymass)
biomass19 <- biomass19 %>% select(SampleID_size, SampleID, DryMass_mg, biomassUncertainty)
test <- drop_na(biomass19)

str(sampleid_size_drymass)
str(test)

#### merge & format data ###########
mergedData$biomassUncertainty <- ""
combdata <- rbind(mergedData, test) # merging biomass data fra 2018 and the 2019 sequencing run

names(metadata)
data <- merge(metadata, combdata, by = "SampleID")
names(data)

# reformat to fit occurrence core standard for biomass - will be changed to MeasurementOrFact later in the script
data <- data %>% rename(organismQuantity = DryMass_mg)
names(data)
data$organismQuantityType <- "mgbiomass"

### merge with PCRID ############
test <- rbind(sampleid_size_pcrid_concentration, pcrid_qubit_19)
data <- merge(data, test, by = "SampleID_size")