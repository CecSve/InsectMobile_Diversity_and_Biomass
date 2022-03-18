# Standardising sampling data

### sampling data ####
metadata_1 <- read.delim("data/sampling_data/SamplingEvent.txt")
metadata_2 <- read.delim("data/sampling_data/SampleReg_20042020.txt")

# make sure the columns are identical between the two dataframes
names(metadata_1) # remove exmpty and unneccesary columns
metadata_1 <- metadata_1 %>% select(SampleID, PID, subLandUseType, Date, StartTime, EndTime, Wind, Temperature)

names(metadata_2)
metadata_2 <- metadata_2 %>% select(PilotTripID, PID, SubLandUseType, Date, StartTime, EndTime, Wind, Temperature)
metadata_2 <- metadata_2 %>% rename(SampleID = PilotTripID)
metadata_2 <- metadata_2 %>% rename(subLandUseType = SubLandUseType)

#### merge metadata ####
metadata <- merge(metadata_1, metadata_2, all = TRUE)

#### coordinates and other locality data ####

labeldata_1 <- read.delim("data/sampling_data/sampling_event_label_data.txt")
str(labeldata_1)
labeldata_2 <- read.delim("data/sampling_data/ethanol_labels_placenames_landuse_2019.txt")
str(labeldata_2)

##### align label data and merge #####
# first select the same columns as in label data 1
labeldata_2 <- labeldata_2 %>% select(Sample.ID, LandUseType, StartTime.y, EndTime.y, DOFAtlasQuadrantID, name, utm_x, utm_y)

labeldata_1 <- labeldata_1 %>% select(SampleID_size, LandUSeType, StartTime, EndTime, DOFAtlasQuadrantID, Location, utm_x, utm_y)

# rename columns to match
labeldata_2 <- labeldata_2 %>% rename(SampleID_size = Sample.ID)
labeldata_1 <- labeldata_1 %>% rename(LandUseType = LandUSeType)
labeldata_2 <- labeldata_2 %>% rename(Location = name)
labeldata_2 <- labeldata_2 %>% rename(StartTime = StartTime.y)
labeldata_2 <- labeldata_2 %>% rename(EndTime = EndTime.y)
str(labeldata_1)
str(labeldata_2) # now they match

#merge
labeldata <- merge(labeldata_1, labeldata_2, all = T)
#labeldata is the data used to generate the labels for the physical samples in the Natural History Museum of Denmarks ethanol collection (label within each sample tube). The samples are stored under ACC.NO. 2018-EN-001

###### transform utm into decimal degrees ######

# the few lines below used to work but due to updates on the rgdal package it crashes

#utmcoor<-SpatialPoints(cbind(labeldata$utm_x,labeldata$utm_y), proj4string=CRS("+init=epsg:25832"))
#longlatcoor<-sp::spTransform(utmcoor,CRS("+proj=longlat"))
#labeldata$decimalLongitude <- coordinates(longlatcoor)[,1]
#labeldata$decimalLatitude <- coordinates(longlatcoor)[,2]

# work-around solution
test <- labeldata
coordinates(test) <- c("utm_x", "utm_y")
proj4string(test) <- CRS("+init=epsg:25832")
summary(test)

test2 <- spTransform(test, CRS("+proj=longlat"))
summary(test2)

labeldata$decimalLongitude <- coordinates(test2)[,1]
labeldata$decimalLatitude <- coordinates(test2)[,2]

### standardise the sublanduse category ####
unique(data$subLandUseType) # 89 unique

test <- data %>%
  mutate(subLandUseType=replace(subLandUseType, subLandUseType=="strandeng, mose, eng, sø", "marsh, bog, meadow, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="hede", "heath")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="mose, eng, sø", "bog, meadow, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="eng, mose", "meadow, bog")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="hede, overdrev", "heath, pasture")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="intensive, markblok", "intensive, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="intensive, semi-intensive, markblok", "intensive, semi-intensive, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="overdrev", "pasture")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="mose, sø", "bog, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="overdrev, intensive, extensive", "pasture, intensive, extensive")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="tør, mark, våd", "dry, agriculture, wet")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="sø", "lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="intesive, extensive, markblok", "intensive, extensive, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="strandeng, eng", "marsh, meadow")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="eng, mose, sø", "bog, meadow, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="eng", "meadow")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="sø, mose, eng", "bog, meadow, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="sø, mose", "bog, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="strandeng, sø, mose", "marsh, lake, bog")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="mose, eng", "meadow, bog")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="eng, mose, strandeng, sø", "marsh, bog, meadow, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="pasture, agricultureblock", "pasture, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="intensive, markblock", "intensive, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="strandeng, sø", "marsh, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="intensive, agricultureblock", "intensive, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="intensive, extensive, agricultureblock", "intensive, extensive, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="agricultureblock, extensive, intensive, semi_intensive", "agricultural field, extensive, intensive, semi_intensive")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="0", "")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="hede/overdrev", "heath, pasture")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="strandeng, mose, sø", "marsh, bog, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="overdrev, strandeng, sø", "pasture, marsh, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="overdrev, markblok", "pasture, agricultural field")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="eng, mose, sø, strandeng", "marsh, bog, meadow, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="intesive", "intensive")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="sø, strandeng", "lake, marsh")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="strandeng", "marsh")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="eng, sø", "meadow, lake")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="eng, overdrev", "meadow, pasutre")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="mose", "bog")) %>% mutate(subLandUseType=replace(subLandUseType, subLandUseType=="markblok, extensive, intensive, semi-intensive", "agricultural field, extensive, intensive, semi_intensive"))

### time & date formatting ##########################
test <- data
test <- test %>% separate(Date, c("day", "month", "year"))
test <- test %>%
  mutate(year=replace(year, year=="2010", "2019")) # correct typo
test <- test %>% separate(StartTime, c("hour", "minute", "second"))
test <- test %>% separate(EndTime, c("hour2", "minute2", "second2"))

str(test)
names(test)
test[, c(4:12)] <- test[, c(4:12)] %>% mutate_if(is.character,as.integer) # change values into integers instead of characters
str(test)

Sys.timezone(location = TRUE)

eventStart <- test %>% 
  select(year, month, day, hour, minute) %>% 
  mutate(eventStart = make_datetime(year, month, day, hour, minute, tz = "Europe/Paris")) %>% select(eventStart)

eventEnd <- test %>% 
  select(year, month, day, hour2, minute2) %>% 
  mutate(eventEnd = make_datetime(year, month, day, hour2, minute2, tz = "Europe/Paris")) %>% select(eventEnd)

test <- cbind(test, eventStart, eventEnd)
str(test)

# create eventTime interval
data$eventTime <- paste(test$eventStart, test$eventEnd, sep="/")
data <- data %>% rename(eventDate = Date)