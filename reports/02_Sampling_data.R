# Standardising sampling data

### sampling data ####
metadata_1 <- read.delim("data/sampling_data/SamplingEvent.txt")
metadata_2 <- read.delim("data/sampling_data/SampleReg_20042020.txt") # update sheet with PilotNotes!

# make sure the columns are identical between the two dataframes
names(metadata_1) # remove exmpty and unneccesary columns
metadata_1 <- metadata_1 %>% select(SampleID, PID, DOFAtlasQuadrantID, subLandUseType, Date, StartTime, EndTime, Wind, Temperature, Notes, PilotNotes) # we will not include google map links since route information will be added later from GIS

names(metadata_2)
metadata_2 <- metadata_2 %>% select(PilotTripID, PID, DOFAtlasQuadrantID, SubLandUseType, Date, StartTime, EndTime, Wind, Temperature)
metadata_2 <- metadata_2 %>% rename(SampleID = PilotTripID)
metadata_2 <- metadata_2 %>% rename(subLandUseType = SubLandUseType)

#### merge metadata ####
metadata <- merge(metadata_1, metadata_2, all = TRUE)

# add routeid and sampleid link for both years
routesample2018 <- read_delim("data/sampling_data/pilotTripIdToRouteID2018.txt", 
                          delim = ";", escape_double = FALSE, 
                          trim_ws = TRUE)

routesample2019 <- read_delim("data/sampling_data/pilotTripIdToRouteID2019.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

names(routesample2019)

routesample2019 <- routesample2019 %>% rename(RouteID_JB = `Rute_adresse-ID`) %>% select(RouteID_JB, PIDRouteID)

# add a PIDRouteID to metadata
metadata <- merge(metadata, routesample2018, by.x = "SampleID", by.y = "PilotTripID", all.x = T)

metadata$PIDRouteID <- metadata$SampleID 
metadata <- metadata %>% mutate(PIDRouteID = gsub("A", "", PIDRouteID)) %>% mutate(PIDRouteID = gsub("B", "", PIDRouteID))

metadata <- merge(metadata, routesample2019, by = c("PIDRouteID"), all.x = T)

# make a combined RouteID_JB column
metadata <- metadata %>%
  mutate(RouteID_JB = coalesce(RouteID_JB.x,RouteID_JB.y))

# remove the redundant columns
metadata <- metadata %>% select(-RouteID_JB.x, -RouteID_JB.y)

###### add coordinates for route centroids ######
coord_1 <- read_delim("data/sampling_data/DK_ruter2018_pkt_koordinater.txt", 
                      delim = "\t", escape_double = FALSE, 
                      locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)
coord_2 <- read_delim("data/sampling_data/ruter2019_koordinater_DOFKvadr.txt", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE) # the route IDs a re currently wrong - the way they are added here, we have no idea

# check which columns match
names(coord_1)
names(coord_2)

# align columns
coord_1 <- coord_1 %>% select(routeID, utm_x, utm_y)
coord_2 <- coord_2 %>% select(routeID, utm_x, utm_y)

# combine coordinates
coords <- merge(coord_1, coord_2, all = T)
str(coords)

###### transform utm into decimal degrees ######

# the few lines below used to work but due to updates on the rgdal package it crashes

#utmcoor<-SpatialPoints(cbind(coords$utm_x,coords$utm_y), proj4string=CRS("+init=epsg:25832"))
#longlatcoor<-sp::spTransform(utmcoor,CRS("+proj=longlat"))
#labeldata$decimalLongitude <- coordinates(longlatcoor)[,1]
#labeldata$decimalLatitude <- coordinates(longlatcoor)[,2]

# work-around solution
test <- coords
coordinates(test) <- c("utm_x", "utm_y")
proj4string(test) <- CRS("+init=epsg:25832")
summary(test)

test2 <- spTransform(test, CRS("+proj=longlat"))
summary(test2)

coords$decimalLongitude <- coordinates(test2)[,1]
coords$decimalLatitude <- coordinates(test2)[,2]

### examine differences in samples ##########################
setdiff(coords$SampleID, metadata$routeID) # all metadata are in labeldata
setdiff(metadata$SampleID, coords$SampleID) # not all metadata are in labeldata

### time & date formatting ########################## NOT DONE! still also need to merge it all in the end and make sure no data is missing
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