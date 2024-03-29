# Standardising environmental data

# preparing DK landuse data prepared by JB in GIS. README file in H:\Documents\Insektmobilen\Data\Arealanvendelse_Århus\2018_bufferzones_data\Final_buffers_2018 (Cecilie Svenningsens work drive). README should be included in final git submission. 

# load libraries required for reformatting and merging data
library(tidyverse) # data wrangling
library(readr) # import data
library(stringr) # character string/text/natural language processing tools for pattern searching
library(data.table) # aggregation of large data
library(ggpubr) # vizualisation
library(tidyr) # data wrangling
library(vegan) # for calculating heterogeneity

### checking we have all data ####################

sampling_data_cleaned <- readRDS("data/sampling_data/sampling_data_cleaned.rds")
sampling_data_cleaned$Year <- lubridate::year(sampling_data_cleaned$Date)
table(sampling_data_cleaned$Year) # contains data for both 2018 and 2019

#get rid of NAs - empty rows
sampling_data_cleaned <- sampling_data_cleaned %>% filter(!is.na(Year))
1582-1358 # removes 224 samples without metadata information
(1358/1582)*100 # 85% have samling metadata for further analysis

#get routes sampled in each year
routes2018 <- unique(sampling_data_cleaned$PIDRouteID[sampling_data_cleaned$Year==2018])
routes2018_ID_JB <- unique(sampling_data_cleaned$RouteID_JB[sampling_data_cleaned$Year==2018])# Jespers ID
routes2019 <- unique(sampling_data_cleaned$PIDRouteID[sampling_data_cleaned$Year==2019])
routes2019_ID_JB <- unique(sampling_data_cleaned$RouteID_JB[sampling_data_cleaned$Year==2019])

#read in an example land use file for 2018
df <- read_delim("data/environmental_data/covariate-data/ruter2018buf1000_areas.txt",
                 "\t", escape_double = FALSE, trim_ws = TRUE)
#check overlap
routes2018[!routes2018_ID_JB %in% df$routeID] # no missing routes between data tables

setdiff(df$routeID, sampling_data_cleaned$RouteID_JB) # in environmental data but not in the sampling data

# setdiff(sampling_data_cleaned$RouteID_JB, df$routeID) # in sampling data but not in the environmental data - routes from 2019, not included yet

#read in an example land use file for 2019
df <- read_delim("data/environmental_data/covariate-data/ruter2019buf1000_areas.txt",
                 "\t", escape_double = FALSE, trim_ws = TRUE) # the proportion column is not being read in right
#check overlap
routes2019[!routes2019 %in% df$routeID]# all there. Good!

#check also against the traffic lights
df <- read_delim("data/environmental_data/covariate-data/DK_TrafficLightsCount.txt","\t", escape_double = FALSE, trim_ws = TRUE)
routes2018[!routes2018_ID_JB %in% df$routeID]
df <- read_delim("data/environmental_data/covariate-data/ruter2019_countStops.txt","\t", escape_double = FALSE, trim_ws = TRUE)
routes2019[!routes2019 %in% df$PIDRouteID]

# many are missing but I guess these are zero stops?? - YES

### make file to relate JB route ID to standard route ID

sampling_data_cleaned <- readRDS("data/sampling_data/sampling_data_cleaned.rds") # this includes the samples that are not sampled

idRelate <- sampling_data_cleaned %>%
              select(PIDRouteID, RouteID_JB) %>%
              unique()

### 2018 data #####

# code below from script 02_DK_environDaat_processing.R from the Biomass git

#### load buffer zone files #### 
#oeko <- read.delim("covariate-data/DK_ruter2018_OekoAreas.txt") # oeko is now part of the buffer zone data
hedge <- read.delim("data/environmental_data/covariate-data/DK_ruter2018_hegnAreas.txt")

urbangreen <- read.delim("data/environmental_data/covariate-data/DK_ruter2018_urbGreenAreas.txt")

nep <- read_delim("data/environmental_data/covariate-data/ruter2018buf500_NEP.txt","\t", escape_double = FALSE, trim_ws = TRUE)

buf_50m <- read_delim("data/environmental_data/covariate-data/ruter2018buf50_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

buf_250m <- read_delim("data/environmental_data/covariate-data/ruter2018buf250_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

buf_500m <- read_delim("data/environmental_data/covariate-data/ruter2018buf500_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

buf_1000m <- read_delim("data/environmental_data/covariate-data/ruter2018buf1000_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

##### WIDE format with tidyr: reformatting environmental data #######
# transform oekodata from long to wide format prior to merging 
#oekocast <- oeko %>% pivot_wider(names_from = bufferDist, values_from = propOeko, names_prefix = "propOeko_")

# rename the land use categories to _organic if it is registered as organic
setDT(buf_50m)
buf_50m[oeko == "1", type := paste0(type, "_", "organic")]
setDT(buf_250m)
buf_250m[oeko == "1", type := paste0(type, "_", "organic")]
setDT(buf_500m)
buf_500m[oeko == "1", type := paste0(type, "_", "organic")]
setDT(buf_1000m)
buf_1000m[oeko == "1", type := paste0(type, "_", "organic")]

#test <- dcast(melt(oeko, id.vars=c("routeID", "bufferDist", "propOeko")), routeID~bufferDist+propOeko)

# transform hedgedata and urbangreen from long to wide format prior to merging - note! multiple value columns
hedgecast <-
  hedge %>% pivot_wider(
    names_from = bufferDist,
    values_from = c(hegnLength, byHegnLength, hegnMeterPerHa, byHegnMeterPerHa)
  )

urbangreencast <-
  urbangreen %>% pivot_wider(
    names_from = bufferDist,
    values_from = c(urbGreenAreaHa, urbGreenPropArea)
  )

# buffer zone data - include column with buffer distance for each dataset
buf_50m$bufferDist <- 50
buf_250m$bufferDist <- 250
buf_500m$bufferDist <- 500
buf_1000m$bufferDist <- 1000

#add on land_use data (following 02_DE script to create DK_environData)
buf_50m$propLand_use <- NA
buf_250m$propLand_use <- NA
buf_500m$propLand_use <- NA
buf_1000m$propLand_use <- NA

# define wetland cover
buf_50m$propLand_use[buf_50m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"
buf_250m$propLand_use[buf_250m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"
buf_500m$propLand_use[buf_500m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"
buf_1000m$propLand_use[buf_1000m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"

# define forest cover
buf_50m$propLand_use[buf_50m$type %in% c("Skov","Skov_organic")] <- "Forest"
buf_250m$propLand_use[buf_250m$type %in% c("Skov","Skov_organic")] <- "Forest"
buf_500m$propLand_use[buf_500m$type %in% c("Skov","Skov_organic")] <- "Forest"
buf_1000m$propLand_use[buf_1000m$type %in% c("Skov","Skov_organic")] <- "Forest"

# define agriculture (farmland) cover - NB there is recreative purposes and some afforestation codes in there that maybe should be removed?
buf_50m$propLand_use[buf_50m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"
buf_250m$propLand_use[buf_250m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"
buf_500m$propLand_use[buf_500m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"
buf_1000m$propLand_use[buf_1000m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"

# define urban cover
buf_50m$propLand_use[buf_50m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"
buf_250m$propLand_use[buf_250m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"
buf_500m$propLand_use[buf_500m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"
buf_1000m$propLand_use[buf_1000m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"

# define open uncultivated land cover
#temp <- subset(buf_1000m,type=="Hede")
#nrow(temp)
#max(temp$areaProportion)

# based on the code above we will make a separate category for heathland
buf_50m$propLand_use[buf_50m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"
buf_250m$propLand_use[buf_250m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"
buf_500m$propLand_use[buf_500m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"
buf_1000m$propLand_use[buf_1000m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"

buf_50m$propLand_use[buf_50m$type %in% c("Hede", "Hede_organic")] <- "Heathland"
buf_250m$propLand_use[buf_250m$type %in% c("Hede", "Hede_organic")] <- "Heathland"
buf_500m$propLand_use[buf_500m$type %in% c("Hede", "Hede_organic")] <- "Heathland"
buf_1000m$propLand_use[buf_1000m$type %in% c("Hede", "Hede_organic")] <- "Heathland"

#unspecified category
buf_50m$propLand_use[buf_50m$type %in% "andet"] <- "Unspecified land cover"
buf_250m$propLand_use[buf_250m$type %in% "andet"] <- "Unspecified land cover"
buf_500m$propLand_use[buf_500m$type %in% "andet"] <- "Unspecified land cover"
buf_1000m$propLand_use[buf_1000m$type %in% "andet"] <- "Unspecified land cover"

#subset to the above land-uses for each buffer
# 50
table(buf_50m$propLand_use)
output50 <- subset(buf_50m,!is.na(propLand_use))

#250
table(buf_250m$propLand_use)
output250 <- subset(buf_250m,!is.na(propLand_use))

#500
table(buf_500m$propLand_use)
output500 <- subset(buf_500m,!is.na(propLand_use))

#1000
table(buf_1000m$propLand_use)
output1000 <- subset(buf_1000m,!is.na(propLand_use))

#cast the data
outputCast50 <- reshape2::dcast(output50,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)
outputCast250 <- reshape2::dcast(output250,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)
outputCast500 <- reshape2::dcast(output500,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)
outputCast1000 <- reshape2::dcast(output1000,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)

#merge buffer zone cast outputs
outputCast <- merge(outputCast50, outputCast250, by = "routeID")
outputCast <- merge(outputCast, outputCast500, by = "routeID")
outputCast <- merge(outputCast, outputCast1000, by = "routeID")

#merge with oeko, urbangreen, and hedge cast outputs
#outputCast <- merge(outputCast, oekocast, by = "routeID")
outputCast <- merge(outputCast, hedgecast, by = "routeID")
outputCast <- merge(outputCast, urbangreencast, by = "routeID")

#add on traffic light information
trafficDF <- read_delim("data/environmental_data/covariate-data/DK_TrafficLightsCount.txt","\t", escape_double = FALSE, trim_ws = TRUE)
outputCast$Num_trafficLights <- trafficDF$Num_trafficLights[match(outputCast$routeID,
                                                                  trafficDF$routeID)]
#missing values are zeros
outputCast$Num_trafficLights[is.na(outputCast$Num_trafficLights)] <- 0
#check associated with urban cover
qplot(Urban_250, Num_trafficLights, data=outputCast)

write.table(outputCast,file="data/environmental_data/covariate-data/environData_2018_DK.txt",sep="\t")

### 2019 data ####

# code below from script 02_DK_environData_processing.R from the Biomass git

#### load buffer zone files #### 

#in 2019, we have a separate oeko file
oeko <- read.delim("data/environmental_data/covariate-data/ruter2019_OekoAreas.txt")

hedge <- read.delim("data/environmental_data/covariate-data/ruter2019_hegnAreas.txt")

urbangreen <- read.delim("data/environmental_data/covariate-data/ruter2019_urbGreenAreas.txt")

buf_50m <- read_delim("data/environmental_data/covariate-data/ruter2019buf50_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

buf_250m <- read_delim("data/environmental_data/covariate-data/ruter2019buf250_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

buf_500m <- read_delim("data/environmental_data/covariate-data/ruter2019buf500_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

buf_1000m <- read_delim("data/environmental_data/covariate-data/ruter2019buf1000_areas.txt","\t", escape_double = FALSE, trim_ws = TRUE)

##### WIDE format with tidyr: reformatting environmental data #######
# transform oekodata from long to wide format prior to merging 

oekocast <- oeko %>% pivot_wider(names_from = bufferDist, values_from = oekoPropArea, names_prefix = "propOeko_")

# transform hedgedata and urbangreen from long to wide format prior to merging - note! multiple value columns
hedgecast <-
  hedge %>% pivot_wider(
    names_from = bufferDist,
    values_from = c(hegnLength, byHegnLength, hegnMeterPerHa, byHegnMeterPerHa)
  )

urbangreencast <-
  urbangreen %>% pivot_wider(
    names_from = bufferDist,
    values_from = c(urbGreenAreaHa, urbGreenPropArea)
  )

# buffer zone data - include column with buffer distance for each dataset
buf_50m$bufferDist <- 50
buf_250m$bufferDist <- 250
buf_500m$bufferDist <- 500
buf_1000m$bufferDist <- 1000

#add on land_use data
buf_50m$propLand_use <- NA
buf_250m$propLand_use <- NA
buf_500m$propLand_use <- NA
buf_1000m$propLand_use <- NA

# define wetland cover
buf_50m$propLand_use[buf_50m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"
buf_250m$propLand_use[buf_250m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"
buf_500m$propLand_use[buf_500m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"
buf_1000m$propLand_use[buf_1000m$type %in% c("Sø", "Mose", "Sø_organic", "Mose_organic")] <- "Wetland"

# define forest cover
buf_50m$propLand_use[buf_50m$type %in% c("Skov","Skov_organic")] <- "Forest"
buf_250m$propLand_use[buf_250m$type %in% c("Skov","Skov_organic")] <- "Forest"
buf_500m$propLand_use[buf_500m$type %in% c("Skov","Skov_organic")] <- "Forest"
buf_1000m$propLand_use[buf_1000m$type %in% c("Skov","Skov_organic")] <- "Forest"

# define agriculture (farmland) cover - NB there is recreative purposes and some afforestation codes in there that maybe should be removed?
buf_50m$propLand_use[buf_50m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"
buf_250m$propLand_use[buf_250m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"
buf_500m$propLand_use[buf_500m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"
buf_1000m$propLand_use[buf_1000m$type %in% c("Ekstensiv", "Ekstensiv_organic", "Intensiv", "Intensiv_organic", "Markblok", "Markblok_organic", "Semi-intensiv", "Semi-intensiv_organic")] <- "Agriculture"

# define urban cover
buf_50m$propLand_use[buf_50m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"
buf_250m$propLand_use[buf_250m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"
buf_500m$propLand_use[buf_500m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"
buf_1000m$propLand_use[buf_1000m$type %in% c("Lav bebyggelse", "Erhverv", "Høj bebyggelse", "Bykerne")] <- "Urban"

# define open uncultivated land cover
#temp <- subset(buf_1000m,type=="Hede")
#nrow(temp)
#max(temp$areaProportion)

# based on the code above we will make a separate category for heathland
buf_50m$propLand_use[buf_50m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"
buf_250m$propLand_use[buf_250m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"
buf_500m$propLand_use[buf_500m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"
buf_1000m$propLand_use[buf_1000m$type %in% c("Overdrev", "Overdrev_organic", "Eng", "Eng_organic", "Strandeng", "Strandeng_organic")] <- "Open uncultivated land"

buf_50m$propLand_use[buf_50m$type %in% c("Hede", "Hede_organic")] <- "Heathland"
buf_250m$propLand_use[buf_250m$type %in% c("Hede", "Hede_organic")] <- "Heathland"
buf_500m$propLand_use[buf_500m$type %in% c("Hede", "Hede_organic")] <- "Heathland"
buf_1000m$propLand_use[buf_1000m$type %in% c("Hede", "Hede_organic")] <- "Heathland"

#unspecified category
buf_50m$propLand_use[buf_50m$type %in% "Andet"] <- "Unspecified land cover"
buf_250m$propLand_use[buf_250m$type %in% "Andet"] <- "Unspecified land cover"
buf_500m$propLand_use[buf_500m$type %in% "Andet"] <- "Unspecified land cover"
buf_1000m$propLand_use[buf_1000m$type %in% "Andet"] <- "Unspecified land cover"

#subset to the above land-uses for each buffer
# 50
table(buf_50m$propLand_use)
output50 <- subset(buf_50m,!is.na(propLand_use))

#250
table(buf_250m$propLand_use)
output250 <- subset(buf_250m,!is.na(propLand_use))

#500
table(buf_500m$propLand_use)
output500 <- subset(buf_500m,!is.na(propLand_use))

#1000
table(buf_1000m$propLand_use)
output1000 <- subset(buf_1000m,!is.na(propLand_use))

#cast the data
outputCast50 <- reshape2::dcast(output50,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)
outputCast250 <- reshape2::dcast(output250,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)
outputCast500 <- reshape2::dcast(output500,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)
outputCast1000 <- reshape2::dcast(output1000,routeID~propLand_use+bufferDist,value.var="areaProportion",fun=sum,na.rm=T)

#merge buffer zone cast outputs
outputCast <- merge(outputCast50, outputCast250, by = "routeID")
outputCast <- merge(outputCast, outputCast500, by = "routeID")
outputCast <- merge(outputCast, outputCast1000, by = "routeID")

#merge with oeko, urbangreen, and hedge cast outputs
#outputCast <- merge(outputCast, oekocast, by = "routeID")?????
outputCast <- merge(outputCast, hedgecast, by = "routeID")
outputCast <- merge(outputCast, urbangreencast, by = "routeID")

#add on traffic light information
trafficDF <- read_delim("data/environmental_data/covariate-data/ruter2019_countStops.txt","\t", escape_double = FALSE, trim_ws = TRUE)
outputCast$Num_trafficLights <- trafficDF$COUNT_STOPS[match(outputCast$routeID,
                                                                  trafficDF$PIDRouteID)]
#missing values are zeros
outputCast$Num_trafficLights[is.na(outputCast$Num_trafficLights)] <- 0
#check associated with urban cover
qplot(Urban_250, Num_trafficLights, data=outputCast)

write.table(outputCast,file="data/environmental_data/covariate-data/environData_2019_DK.txt",sep="\t")

### combine all #########################################

#get data for each year
environData2018 <-read.delim("data/environmental_data/covariate-data/environData_2018_DK.txt")%>%
                      add_column(Year = 2018)

environData2019 <-read.delim("data/environmental_data/covariate-data/environData_2019_DK.txt")%>%
                      add_column(Year = 2019)

#check all present - the following should be empty characters
#CS should check
names(environData2018)[!names(environData2018) %in% names(environData2019)] # empty
names(environData2019)[!names(environData2019) %in% names(environData2018)] # empty

#route IDs are different!!!!

# add routeid and sampleid link for both years
routesample2018 <- read_delim("data/sampling_data/pilotTripIdToRouteID2018.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

routesample2019 <- read_delim("data/sampling_data/pilotTripIdToRouteID2019.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

names(routesample2019)

# merge data to have all ID columns in the same table
routedata2019 <- merge(environData2019, routesample2019, by.x = "routeID", by.y = "PIDRouteID")

# keep only the true route ID column
environData2019 <-
  routedata2019 %>% rename(RouteID_JB = `Rute_adresse-ID`) %>% select(-routeID,-RouteID,-kvnr,-`Adresse_adresse-ID`,-PID) %>% select(RouteID_JB, everything()) %>% rename(routeID = RouteID_JB)

#combine all together and get unique values per route 
environData <- bind_rows(environData2018, environData2019) %>%
                  unique()

# in some cases there are two route IDs because the land cover has changed between years. Be mindful of this when merging with sampling data
table(environData$routeID)

### heterogeneity calculation #########################################

#### shannon diversity #####

routeDiversity <- environData %>%
                    select(routeID,
                           starts_with("Forest"),
                           starts_with("Heathland"),
                           starts_with("Open.uncultivated"),
                           starts_with("Agriculture"),
                           starts_with("Urban"),
                           starts_with("Wetland")) %>%
                    pivot_longer(contains("_"),
                                 names_to = "land_use",
                                 values_to = "prop") %>%
                    separate(land_use, c("land_use", "buffer"), sep="_") %>%
                    group_by(routeID, buffer) %>%
                    summarise(diversity = diversity(prop)) %>%
                    ungroup() %>%
                    pivot_wider(everything(),
                                names_from = "buffer",
                                values_from = "diversity",
                                names_prefix = "Diversity_")
  
#add onto main file
environData <- left_join(environData, routeDiversity,
                         by = "routeID")
                    
write.table(environData,file="data/environmental_data/covariate-data/DK_environData.txt", sep="\t", col.names = T, row.names = F)            

                  