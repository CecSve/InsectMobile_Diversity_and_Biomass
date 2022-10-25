
#### Load required libraries ################################################################
library(cowplot) # for visuals
library(ggplot2) # for visuals
library(ggpubr) # for visuals
library(scales) # for visuals
library(ggpmisc) # ggplot extensions, for visuals
library(grid) # for visuals
library(gridExtra) # for visuals
library(car) # Companion to Applied Regression
library(lme4) # Fit linear and generalized linear mixed-effects models
library(lmerTest) # Provides p-values in type I, II or III anova and summary tables for lmer model fits (cf. lme4) 
library(nlme) # Fit and compare Gaussian linear and nonlinear mixed-effects models
library(effects) # Graphical and tabular effect displays, e.g., of interactions, for various statistical models with linear predictors
library(MuMIn)#AIC, R2
library(multcomp) # generalized linear hypothesis, to compare between effects of land cover
library(sjPlot) # to plot nice tables
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(tidyverse) # for data wrangling
library(lubridate) # for date and time wrangling

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified") # if including unspecified/other category

### Load data ###################################################
# NB! choose which country to run the analysis on and load the corresponding data

# load Danish data
allInsects <- read_rds("data/cleaned_data/data_richness_BINs.RDS")

# load German data
# allInsects <- read.delim("cleaned-data/DE_allInsects.txt")

# change the five land covers to be 0-100 instead of 0-1
allInsects_trans_landcovers <- allInsects
allInsects[,33:60] <- allInsects[,33:60]*100

#### General check-ups and adding variables for analysis ############################################

# remove samples with high uncertainty for the biomass estimation, NA biomass samples, samples without a date
# allInsects <- allInsects %>% filter(!biomassUncertainty == "high") # DK: 77 samples with high biomass uncertainty (different scales used in the lab)

##### Land cover check ################

#check land covers within a buffer of the same size
data1000m <- allInsects[,grepl("_1000",names(allInsects))]
summary(apply(data1000m[,1:5],1,sum))#does not exceed 100
data500m <- allInsects[,grepl("_500",names(allInsects))]
summary(apply(data500m[,1:5],1,sum))#does not exceed 100
data250m <- allInsects[,grepl("_250",names(allInsects))]
summary(apply(data250m[,1:5],1,sum))
data50m <- allInsects[,grepl("_50",names(allInsects))]
data50m <- data50m[,!grepl("_500",names(data50m))]
summary(apply(data50m[,1:5],1,sum))#does not exceed 100

# sumary of biomass and richness (estimated) by year
allInsects %>%
  group_by(Year) %>%
  summarize(min_biomass = min(totalBiomass_mg),
            max_biomass = max(totalBiomass_mg),
            mean_biomass = mean(totalBiomass_mg),
            median_biomass = median(totalBiomass_mg),
            min_richness = min(richness_est),
            max_richness = max(richness_est)) 

# summary statistics for time band - not used for anything, but to visualise the data
allInsects %>%
  group_by(Time_band, Temperature) %>% 
  summarize(min = min(totalBiomass_mg),
            q1 = quantile(totalBiomass_mg, 0.25),
            median = median(totalBiomass_mg),
            mean = mean(totalBiomass_mg),
            q3 = quantile(totalBiomass_mg, 0.75),
            max = max(totalBiomass_mg),
            sd = sd(totalBiomass_mg))

### Analysis ##############################################

##### Supplementary figure: the data ##########################

###### richness #####
qU <- ggplot(allInsects,aes(x=Urban_1000,y=(richness_est)))+
  geom_point(col=landuseCols[1])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70") + scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Urban cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qF <- ggplot(allInsects,aes(x=Agriculture_1000,y=(richness_est)))+
  geom_point(col=landuseCols[2])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Farmland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qD <- ggplot(allInsects,aes(x=Open.uncultivated.land_1000,y=(richness_est)))+
  geom_point(col=landuseCols[3])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    #limits = c(0, 0.16),
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("Estimated richness") + labs(subtitle = "Grassland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 20))

qW <- ggplot(allInsects,aes(x=Wetland_1000,y=(richness_est)))+
  geom_point(col=landuseCols[4])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Wetland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qFo <- ggplot(allInsects,aes(x=Forest_1000,y=(richness_est)))+
  geom_point(col=landuseCols[5])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Forest cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

supp_fig <- cowplot::plot_grid(qU,qF,qD,qW,qFo,ncol=1)

#cowplot::save_plot("plots/DK_Landcover_percent.tiff", supp_fig, base_width = 5, base_height = 14, dpi = 800)

cowplot::save_plot(
  "plots/DK_Landcover_richness_percent.png",
  supp_fig,
  base_width = 5,
  base_height = 14,
  dpi = 800
)

###### biomass #####
qU <- ggplot(allInsects,aes(x=Urban_1000,y=(totalBiomass_mg)))+
  geom_point(col=landuseCols[1])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70") + scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Urban cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qF <- ggplot(allInsects,aes(x=Agriculture_1000,y=(totalBiomass_mg)))+
  geom_point(col=landuseCols[2])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Farmland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qD <- ggplot(allInsects,aes(x=Open.uncultivated.land_1000,y=(totalBiomass_mg)))+
  geom_point(col=landuseCols[3])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    #limits = c(0, 0.16),
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("Estimated biomass") + labs(subtitle = "Grassland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 20))

qW <- ggplot(allInsects,aes(x=Wetland_1000,y=(totalBiomass_mg)))+
  geom_point(col=landuseCols[4])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Wetland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qFo <- ggplot(allInsects,aes(x=Forest_1000,y=(totalBiomass_mg)))+
  geom_point(col=landuseCols[5])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Forest cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

supp_fig <- cowplot::plot_grid(qU,qF,qD,qW,qFo,ncol=1)

#cowplot::save_plot("plots/DK_Landcover_percent.tiff", supp_fig, base_width = 5, base_height = 14, dpi = 800)

cowplot::save_plot(
  "plots/DK_Landcover_biomass_percent.png",
  supp_fig,
  base_width = 5,
  base_height = 14,
  dpi = 800
)

##### Figure: pie chart#####################################
###### richness ####
routeMeans <- allInsects %>% 
  group_by(RouteID_JB) %>%
  dplyr::summarise(meanEstRichness = mean(richness_est))

allInsects <- inner_join(allInsects,routeMeans,by="RouteID_JB") 

#remove duplicates
allInsects_pie <- allInsects %>%
  dplyr::select(RouteID_JB,meanEstRichness,Agriculture_1000,
                Forest_1000,Open.uncultivated.land_1000,
                Urban_1000,Wetland_1000) %>%
  distinct()

#fill in missing column
allInsects_pie$totalLand <- apply(allInsects_pie[,3:7],1,sum)
allInsects_pie$Other_1000 <- 1-allInsects_pie$totalLand

#divide up biomass into quantiles
allInsects_pie$BiomassCats <- cut_number(allInsects_pie$meanEstRichness,n=5)

#mean land cover per biomass cats
allInsects_cat <- allInsects_pie %>%
  group_by(BiomassCats) %>%
  dplyr::summarise(Agriculture_1000 = mean(Agriculture_1000),
                   Forest_1000 = mean(Forest_1000),
                   Open.uncultivated.land_1000 = mean(Open.uncultivated.land_1000),
                   Urban_1000 = mean(Urban_1000),
                   Wetland_1000 = mean(Wetland_1000),
                   Other_1000 = mean(Other_1000))

#plot data by biomass categories
allInsects_melt <- gather(allInsects_cat, Land_cover, value, -BiomassCats)

#plot
allInsects_melt <- allInsects_melt %>% mutate(
  Land_cover = fct_relevel(
    Land_cover,
    "Urban_1000",
    "Agriculture_1000",
    "Open.uncultivated.land_1000",
    "Wetland_1000",
    "Forest_1000",
    "Other_1000"))  

biomass.labs <- c("[8.5,75.1]"="9-75 estimated species", "(75.1,106]"="75-106 estimated species", "(106,130]"="106-130 estimated species", "(130,156]"="130-156 estimated species", "(156,255]"="156-255 estimated species")

levels(allInsects_melt$BiomassCats)
forplot <- allInsects_melt %>% filter(BiomassCats %in% c("[8.5,75.1]", "(156,255]"))

fig_pie <- ggplot(forplot,aes(x="",y=value,fill=Land_cover, order = Land_cover))+
  geom_bar(stat="identity")+
  facet_wrap(~BiomassCats, labeller = labeller(BiomassCats=biomass.labs[c(1,5)]), ncol = 1)+
  coord_polar("y")+
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position="left", strip.text = element_text(size = 12), legend.text = element_text(size = 12), legend.key.size = unit(1, 'cm')) + scale_fill_manual(values = landuseCols, name = "", labels = c("Urban", "Farmland", "Grassland", "Wetland", "Forest", "Other")) 

cowplot::save_plot("plots/DK_landcover_estRichness_proportions.png", fig_pie, base_width = 12, base_height = 6, dpi = 800) # cange to .tiff if included in publication

###### biomass ####
routeMeans <- allInsects %>% 
  group_by(RouteID_JB) %>%
  dplyr::summarise(meanEstRichness = mean(totalBiomass_mg))

allInsects <- inner_join(allInsects,routeMeans,by="RouteID_JB") 

#remove duplicates
allInsects_pie <- allInsects %>%
  dplyr::select(RouteID_JB,meanEstRichness,Agriculture_1000,
                Forest_1000,Open.uncultivated.land_1000,
                Urban_1000,Wetland_1000) %>%
  distinct()

#fill in missing column
allInsects_pie$totalLand <- apply(allInsects_pie[,3:7],1,sum)
allInsects_pie$Other_1000 <- 1-allInsects_pie$totalLand

#divide up biomass into quantiles
allInsects_pie$BiomassCats <- cut_number(allInsects_pie$meanEstRichness,n=5)

#mean land cover per biomass cats
allInsects_cat <- allInsects_pie %>%
  group_by(BiomassCats) %>%
  dplyr::summarise(Agriculture_1000 = mean(Agriculture_1000),
                   Forest_1000 = mean(Forest_1000),
                   Open.uncultivated.land_1000 = mean(Open.uncultivated.land_1000),
                   Urban_1000 = mean(Urban_1000),
                   Wetland_1000 = mean(Wetland_1000),
                   Other_1000 = mean(Other_1000))

#plot data by biomass categories
allInsects_melt <- gather(allInsects_cat, Land_cover, value, -BiomassCats)

#plot
allInsects_melt <- allInsects_melt %>% mutate(
  Land_cover = fct_relevel(
    Land_cover,
    "Urban_1000",
    "Agriculture_1000",
    "Open.uncultivated.land_1000",
    "Wetland_1000",
    "Forest_1000",
    "Other_1000"))  

levels(allInsects_melt$BiomassCats)
max(allInsects$totalBiomass_mg)
biomass.labs <- c("[2,59]"="2-59 mg", "(59,107]"="59-107 mg", "(107,154]"="107-154 mg", "(154,253]"="154-253 mg", "(253,6.83e+03]"="253-6831 mg")

forplot <- allInsects_melt %>% filter(BiomassCats %in% c("[2,59]", "(253,6.83e+03]"))

fig_pie <- ggplot(forplot,aes(x="",y=value,fill=Land_cover, order = Land_cover))+
  geom_bar(stat="identity")+
  facet_wrap(~BiomassCats, labeller = labeller(BiomassCats=biomass.labs[c(1,5)]), ncol = 1)+
  coord_polar("y")+
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position="left", strip.text = element_text(size = 12), legend.text = element_text(size = 12), legend.key.size = unit(1, 'cm')) + scale_fill_manual(values = landuseCols, name = "", labels = c("Urban", "Farmland", "Grassland", "Wetland", "Forest", "Other")) 

cowplot::save_plot("plots/DK_landcover_biomass_proportions.png", fig_pie, base_width = 12, base_height = 6, dpi = 800) # cange to .tiff if included in publication

##### Table: linear Mixed Effects Model - Land covers #################
# used cStops instead of cTL for DK data

# add average speed 
allInsects$avg_speed <- (allInsects$Distance_driven/1000)/(allInsects$Time_driven/60) # 1000 to get km and 60 to get hours instead of minutes
table(allInsects$avg_speed)
mean(allInsects$avg_speed, na.rm = T) # we will not use this variable, but it gives us an idea for how fast the routes were driven

###### richness ####
hist(allInsects$richness_est) # data is skewed
hist(sqrt(allInsects$richness_est)) # square root of richness gives normal distribution makes grassland significant, but the pattern is generally the same

# center the land covers to help model convergence if the covariates dont have large numbers - not the prettiest, but gets the job done 
allInsects <- allInsects %>% 
  mutate(cUrban_1000 = median(Urban_1000),
         cFarmland_1000 = median(Agriculture_1000),
         cGrassland_500 = median(Open.uncultivated.land_500),
         cWetland_50 = median(Wetland_50),
         cForest_1000 = median(Forest_1000),)

# testing if glm may be better since response is right skewed
test_model <- glmer(obs_richness ~ 
                      (Urban_1000-cUrban_1000) +
                      (Agriculture_1000-cFarmland_1000) + 
                      (Open.uncultivated.land_500-cGrassland_500) +
                      (Wetland_50-cWetland_50) +
                      (Forest_1000-cForest_1000) +
                      Time_band + 
                      Year +
                      Time_band:cnumberTime + cTL + cyDay + 
                      (1|RouteID_JB) + (1|PID), family = poisson, data = allInsects)
                  
summary(test_model)
anova(test_model)
shapiro.test(resid(test_model))
qqnorm(resid(test_model)) # plot of the normal theoretical quantiles against the exponential data quantiles
qqline(resid(test_model)) # the residuals should fall along this line if they are normally distributed
AICc(test_model)
#check variance inflation factor
vif(test_model)

tab_model(test_model, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Farmland (1000 m)","Grassland (500 m)", "Wetland (50 m)", "Forest (1000 m)", "Time band: evening vs midday", "Year", "Potential stops", "Day of year", "Time within midday", "Time within evening"), digits = 2)

#full and final model - the effect of each land cover on estimated flying insect richnes
full_model <- lmer(richness_est ~ 
                     Urban_1000 +
                     Agriculture_1000 + 
                     Open.uncultivated.land_500+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_1000 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) # data=subset(allInsects, c(Open.uncultivated.land_1000 < 20, Wetland_50 < 40))

summary(full_model)
anova(full_model)
shapiro.test(resid(full_model))
qqnorm(resid(full_model)) # plot of the normal theoretical quantiles against the exponential data quantiles
qqline(resid(full_model)) # the residuals should fall along this line if they are normally distributed
AICc(full_model)
#check variance inflation factor
vif(full_model)

tab_model(full_model, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Farmland (1000 m)","Grassland (500 m)", "Wetland (50 m)", "Forest (1000 m)", "Time band: evening vs midday", "Year", "Potential stops", "Day of year", "Time within midday", "Time within evening"), digits = 2) #, file = "plots/full_model_richness.html"

# a quick plot
sjPlot::plot_model(full_model)

# "Biomass (mg)", "mean DNA conc. (ng/ul)",

###### biomass ####

hist(allInsects$totalBiomass_mg)
hist(log(allInsects$totalBiomass_mg))

#full and final model - the effect of each land cover on estimated flying insect biomass. Notice the different buffer sizes from biomass!
full_model <- lmer(log(totalBiomass_mg+1) ~ 
                     Urban_1000 +
                     Agriculture_50 + 
                     Open.uncultivated.land_1000+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_500 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) # data=subset(allInsects, c(Open.uncultivated.land_1000 < 20, Wetland_50 < 40))

summary(full_model)
AICc(full_model)
#check variance inflation factor
vif(full_model)

tab_model(full_model, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Farmland (50 m)","Grassland (1000 m)", "Wetland (50 m)", "Forest (500 m)", "Time band: evening vs midday", "Year", "Potential stops", "Day of year", "Time within midday", "Time within evening"), digits = 2, file = "plots/full_model_biomass.html")

# a quick plot
sjPlot::plot_model(full_model)

# "Biomass (mg)", "mean DNA conc. (ng/ul)",

### predominant land cover: mean estimated richness ###
# Create new column that picks out samples that are at the extreme in representing a single land cover type (>60% or >80% of a single land cover type) at the buffer zone with the largest effect size

allInsects_totsample <- allInsects

allInsects_totsample$hab50 = 'Mix' # samples that do not have more than 50% of one specific land type
allInsects_totsample$hab50[allInsects_totsample$Agriculture_500>=0.50]<-'Agriculture50'
allInsects_totsample$hab50[allInsects_totsample$Forest_1000>=0.50]<-'Forest50'
allInsects_totsample$hab50[allInsects_totsample$Urban_1000>=0.50]<-'Urban50'
table(allInsects_totsample$hab50) # notice the variation in sample size

allInsects_totsample %>% 
  dplyr::select(hab50, richness_est) %>% 
  dplyr::filter(hab50 == "Agriculture50") %>% 
  dplyr::summarise(mean = mean(richness_est)) # 127.2075

allInsects_totsample %>% 
  dplyr::select(hab50, richness_est) %>% 
  dplyr::filter(hab50 == "Forest50") %>% 
  dplyr::summarise(mean = mean(richness_est)) # 135.755

allInsects_totsample %>% 
  dplyr::select(hab50, richness_est) %>% 
  dplyr::filter(hab50 == "Urban50") %>% 
  dplyr::summarise(mean = mean(richness_est)) # 63.50796

allInsects_totsample %>% 
  dplyr::select(hab50, richness_est) %>% 
  dplyr::filter(hab50 == "Mix") %>% 
  dplyr::summarise(mean = mean(richness_est)) # 106.989

##### Table: multcomp landcovers - simple model ##########################

###### richness ####

full_model <- lmer(richness_est ~ 
                     Urban_1000 +
                     Agriculture_1000 + 
                     Open.uncultivated.land_500+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_1000 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) # data=subset(allInsects, c(Open.uncultivated.land_1000 < 20, Wetland_50 < 40))

# pairwise comparison to urban
pair.ht <- glht(full_model, linfct = c("Urban_1000 - Agriculture_1000 = 0", "Urban_1000 - Forest_1000 = 0", "Urban_1000 - Wetland_50 = 0", "Urban_1000 - Open.uncultivated.land_500 = 0"))
summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to farmland
pair.ht <- glht(full_model, linfct = c("Agriculture_1000 - Forest_1000 = 0", "Agriculture_1000 - Wetland_50  = 0", "Agriculture_1000 - Open.uncultivated.land_500 = 0"))
summary(pair.ht) # semi-natural covers have higher biomass than farmland, but it is only significant for grassland, urban has significantly lower biomass
confint(pair.ht)

# pairwise comparison to grassland
pair.ht <- glht(full_model, linfct = c("Open.uncultivated.land_500 - Wetland_50 = 0", "Open.uncultivated.land_500 - Forest_1000 = 0"))
summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to wetland
pair.ht <- glht(full_model, linfct = c("Wetland_50 - Forest_1000 = 0"))
summary(pair.ht) 
confint(pair.ht)

###### Plot: coefficient plot ######################################################
#point estimate
myConfint <- confint(full_model)

coefDF <- data.frame(Landcover = names(fixef(full_model)),
                     estimate = as.numeric(fixef(full_model)),
                     lowerCI = as.numeric(myConfint[-c(1:3),1]),
                     upperCI = as.numeric(myConfint[-c(1:3),2]))

#subset to land covers
coefDF <- subset(coefDF, Landcover %in% c("Agriculture_1000","Urban_1000","Open.uncultivated.land_500","Wetland_50","Forest_1000"))

p <-
  coefDF %>% mutate(
    Landcover = fct_relevel(
      Landcover,
      "Urban_1000",
      "Agriculture_1000",
      "Open.uncultivated.land_500",
      "Wetland_50",
      "Forest_1000"
    )) %>% ggplot(aes(Landcover, estimate))

c<-p+scale_size_area(max_size = 1.5)

#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)

coefPlot <-
  d + geom_pointrange(aes(ymin = lowerCI, ymax = upperCI, colour = Landcover),
                      size = 1.5) + theme_minimal_grid() + theme(
                        legend.title = element_blank(),
                        legend.key = element_rect(size = 0.1),
                        legend.key.size = unit(1, 'cm')
                      ) + labs(x = "", y = "Effect on estimated flying arthropod\n richness & 95% CIs") + scale_x_discrete(
                        labels = c(
                          "Urban_1000" = "Urban",
                          "Agriculture_1000" = "Farmland",
                          "Open.uncultivated.land_500" = "Grassland",
                          "Wetland_50" = "Wetland",
                          "Forest_1000" = "Forest"
                        )
                      ) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.3,hjust=1), legend.position = "none", axis.title.y = element_text(size = 10)) + scale_colour_manual(values = landuseCols) + coord_fixed(0.03)

cowplot::save_plot("plots/DK_coefficient_richness_plot.png", coefPlot, base_width = 10, base_height = 5, dpi = 300)

###### biomass ####

full_model <- lmer(log(totalBiomass_mg+1) ~ 
                     Urban_1000 +
                     Agriculture_50 + 
                     Open.uncultivated.land_1000+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_500 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) # data=subset(allInsects, c(Open.uncultivated.land_1000 < 20, Wetland_50 < 40))

# pairwise comparison to urban
pair.ht <- glht(full_model, linfct = c("Urban_1000 - Agriculture_50 = 0", "Urban_1000 - Forest_500 = 0", "Urban_1000 - Wetland_50 = 0", "Urban_1000 - Open.uncultivated.land_1000 = 0"))
summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to farmland
pair.ht <- glht(full_model, linfct = c("Agriculture_50 - Forest_500 = 0", "Agriculture_50 - Wetland_50  = 0", "Agriculture_50 - Open.uncultivated.land_1000 = 0"))
summary(pair.ht) # semi-natural covers have higher biomass than farmland, but it is only significant for grassland, urban has significantly lower biomass
confint(pair.ht)

# pairwise comparison to grassland
pair.ht <- glht(full_model, linfct = c("Open.uncultivated.land_1000 - Wetland_50 = 0", "Open.uncultivated.land_1000 - Forest_500 = 0"))
summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to wetland
pair.ht <- glht(full_model, linfct = c("Wetland_50 - Forest_500 = 0"))
summary(pair.ht) 
confint(pair.ht)

##### # Plot: coefficient plot ######################################################
#point estimate
myConfint <- confint(full_model)

coefDF <- data.frame(Landcover = names(fixef(full_model)),
                     estimate = as.numeric(fixef(full_model)),
                     lowerCI = as.numeric(myConfint[-c(1:3),1]),
                     upperCI = as.numeric(myConfint[-c(1:3),2]))

#subset to land covers
coefDF <- subset(coefDF, Landcover %in% c("Agriculture_50","Urban_1000","Open.uncultivated.land_1000","Wetland_50","Forest_500"))

p <-
  coefDF %>% mutate(
    Landcover = fct_relevel(
      Landcover,
      "Urban_1000",
      "Agriculture_50",
      "Open.uncultivated.land_1000",
      "Wetland_50",
      "Forest_500"
    )) %>% ggplot(aes(Landcover, estimate))

c<-p+scale_size_area(max_size = 1.5)

#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)

coefPlot <-
  d + geom_pointrange(aes(ymin = lowerCI, ymax = upperCI, colour = Landcover),
                      size = 1.5) + theme_minimal_grid() + theme(
                        legend.title = element_blank(),
                        legend.key = element_rect(size = 0.1),
                        legend.key.size = unit(1, 'cm')
                      ) + labs(x = "", y = "Effect on estimated flying arthropod\n biomass (log+1) & 95% CIs") + scale_x_discrete(
                        labels = c(
                          "Urban_1000" = "Urban",
                          "Agriculture_50" = "Farmland",
                          "Open.uncultivated.land_1000" = "Grassland",
                          "Wetland_50" = "Wetland",
                          "Forest_500" = "Forest"
                        )
                      ) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.3,hjust=1), legend.position = "none", axis.title.y = element_text(size = 10)) + scale_colour_manual(values = landuseCols) + coord_fixed(1)

cowplot::save_plot("plots/DK_coefficient_biomass_plot.png", coefPlot, base_width = 10, base_height = 5, dpi = 300)

##### Table: exclude urban (more than 5% urban) ########################################

###### richness ####
summary(allInsects$Urban_1000) # NB! transformed land covers

#include only sites with less than 5% urban
allInsects_lowUrban <- subset(allInsects, Urban_1000 <5)
summary(allInsects_lowUrban$Urban_1000)

nrow(allInsects_lowUrban)
length(unique(allInsects_lowUrban$SampleID))
length(unique(allInsects_lowUrban$RouteID_JB))

lme1000 <- lmer(richness_est ~ 
                  Agriculture_1000 + 
                  Open.uncultivated.land_500+
                  Wetland_50 +
                  Forest_1000 +
                  Year +
                  Time_band + 
                  Time_band:cnumberTime + 
                  cTL + 
                  cyDay + 
                  (1|RouteID_JB) + (1|PID), data=allInsects_lowUrban)

summary(lme1000)
car::vif(lme1000)

tab_model(lme1000, digits = 2, collapse.ci = F, show.intercept = F, pred.labels = c("Farmland", "Grassland", "Wetland", "Forest", "Year", "Evening vs. midday", "Potential stops", "Day of year", "Time during midday", "Time during evening"), title = "Urban cover > 5% excluded")

# pairwise comparison to farmland
pair.ht <- glht(lme1000, linfct = c("Agriculture_1000 - Forest_250 = 0", "Agriculture_1000 - Wetland_50 = 0", "Agriculture_1000 - Open.uncultivated.land_1000 = 0"))
summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to grassland
pair.ht <- glht(lme1000, linfct = c(
  "Open.uncultivated.land_1000 - Wetland_50 = 0",
  "Open.uncultivated.land_1000 - Forest_250  = 0"))

summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to wetland
pair.ht <- glht(lme1000, linfct = c(
  "Wetland_50 - Forest_250  = 0"))

summary(pair.ht) 
confint(pair.ht)

###### biomass ####

lme1000 <- lmer(log(totalBiomass_mg+1) ~ 
                  Agriculture_1000 + 
                  Open.uncultivated.land_500+
                  Wetland_50 +
                  Forest_1000 +
                  Year +
                  Time_band + 
                  Time_band:cnumberTime + 
                  cTL + 
                  cyDay + 
                  (1|RouteID_JB) + (1|PID), data=allInsects_lowUrban)

summary(lme1000)
car::vif(lme1000)

tab_model(lme1000, digits = 2, collapse.ci = F, show.intercept = F, pred.labels = c("Farmland", "Grassland", "Wetland", "Forest", "Year", "Evening vs. midday", "Potential stops", "Day of year", "Time during midday", "Time during evening"), title = "Urban cover > 5% excluded")

# pairwise comparison to farmland
pair.ht <- glht(lme1000, linfct = c("Agriculture_1000 - Forest_250 = 0", "Agriculture_1000 - Wetland_50 = 0", "Agriculture_1000 - Open.uncultivated.land_1000 = 0"))
summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to grassland
pair.ht <- glht(lme1000, linfct = c(
  "Open.uncultivated.land_1000 - Wetland_50 = 0",
  "Open.uncultivated.land_1000 - Forest_250  = 0"))

summary(pair.ht) 
confint(pair.ht)

# pairwise comparison to wetland
pair.ht <- glht(lme1000, linfct = c(
  "Wetland_50 - Forest_250  = 0"))

summary(pair.ht) 
confint(pair.ht)

### richness AIC check ##############################################
#full and final model - the effect of each land cover on estimated flying insect richnes
full_model <- lmer(richness_est ~ 
                     Urban_1000 +
                     Agriculture_1000 + 
                     Open.uncultivated.land_500+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_1000 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects)

options(na.action = "na.fail")
dd <- MuMIn::dredge(full_model)
subset(dd, delta < 2)

# Visualize the model selection table:
par(mar = c(3,5,6,4))
plot(dd, labAsExpr = TRUE)

# Model average models with delta AICc < 4
#model.avg(dd, subset = delta < 2)

# best model with lowest AIC 
best_model <- lmer(richness_est ~ 
                     Agriculture_1000 + 
                     #Urban_1000 +
                     Open.uncultivated.land_500 + # test if grassland outlier modifies results
                     Wetland_50 +
                     Forest_1000 +
                     Time_band + cTL + Year + cnumberTime +
                     (1|RouteID_JB) + (1|PID), data= allInsects)
# data=subset(allInsects, Open.uncultivated.land_1000 < 0.2)
summary(best_model)
AICc(best_model)
tab_model(best_model, show.intercept = F, pred.labels = c("Farmland (1000 m)", "Grassland (500 m)", "Wetland (50 m)", "Forest (1000 m)", "Time band: midday vs evening", "Potential stops", "Year", "Time within a day"))
r.squaredGLMM(best_model)

# summary for both the full and the best fit model
tab_model(full_model, best_model, title = "Estimated richness", show.intercept = F, digits = 2, dv.labels = c("Full model", "Best fit model"), pred.labels = c("Urban (1000 m)", "Farmland (1000 m)", "Grassland (500 m)", "Wetland (50 m)", "Forest (1000 m)", "Time band: midday vs. evening", "Year: 2019 vs. 2018", "Potential stops", "Day of the year", "Time during evening", "Time during midday", "Time during the day"))

### biomass AIC check ##############################################

#full and final model - the effect of each land cover on estimated flying insect richness
full_model <- lmer(log(totalBiomass_mg+1) ~ 
                     Urban_1000 +
                     Agriculture_50 + 
                     Open.uncultivated.land_1000+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_500 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects)

options(na.action = "na.fail")
dd <- MuMIn::dredge(full_model)
subset(dd, delta < 2)

# Visualize the model selection table:
par(mar = c(3,5,6,4))
plot(dd, labAsExpr = TRUE)

# Model average models with delta AICc < 4
#model.avg(dd, subset = delta < 2)

# best model with lowest AIC 
best_model <- lmer(log(totalBiomass_mg+1) ~ 
                     #Agriculture_1000 + 
                     Urban_1000 +
                     #Open.uncultivated.land_500 + # test if grassland outlier modifies results
                     #Wetland_50 +
                     #Forest_1000 +
                     Time_band + 
                     #cTL + 
                     Year + cnumberTime +
                     (1|RouteID_JB) + (1|PID), data= allInsects)
# data=subset(allInsects, Open.uncultivated.land_1000 < 0.2)
summary(best_model)
AICc(best_model)
tab_model(best_model, show.intercept = F, pred.labels = c("Urban (1000 m)", "Time band: midday vs evening", "Year", "Time within a day"))
 nr.squaredGLMM(best_model)

# summary for both the full and the best fit model
tab_model(full_model, best_model, title = "Estimated biomass (log +1)", show.intercept = F, digits = 2, dv.labels = c("Full model", "Best fit model"), pred.labels = c("Urban (1000 m)", "Farmland (50 m)", "Grassland (1000 m)", "Wetland (50 m)", "Forest (500 m)", "Time band: midday vs. evening", "Year: 2019 vs. 2018", "Potential stops", "Day of the year", "Time during evening", "Time during midday", "Time during the day"))

##### Figure 4: effect plot ##########################

###### richness ##########
full_model <- lmer(richness_est ~ 
                     Urban_1000 +
                     Agriculture_1000 + 
                     Open.uncultivated.land_500+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_1000 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects)

# extract effects
gls1.alleffects <- allEffects(full_model)
plot(gls1.alleffects, 'Year', ylab="richness_est")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(full_model)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(full_model, ~ Agriculture_1000 + Urban_1000 + Open.uncultivated.land_500 + Forest_1000 + Wetland_50 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

# make data frames for each land cover
temp <- effectdata$Agriculture_1000
temp$landcover <- "Agriculture_1000"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Open.uncultivated.land
temp <- effectdata$Open.uncultivated.land_500
temp$landcover <- "Grassland_500"
grass <- temp %>% 
  dplyr::rename(
    propcover = Open.uncultivated.land_500
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Wetland
temp <- effectdata$Wetland_50
temp$landcover <- "Wetland_50"
wet <- temp %>% 
  dplyr::rename(
    propcover = Wetland_50
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Forest
temp <- effectdata$Forest_1000
temp$landcover <- "Forest_1000"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
temp <- effectdata$`Time_band:cnumberTime`
#temp$landcover <- "Time_band:cnumberTime"
time <- temp

# combine the data frames
test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_1000",
    "Grassland_500",
    "Wetland_50",
    "Forest_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover",
      "Farmland cover",
      "Grassland cover",
      "Wetland cover",
      "Forest cover"
    )
  ) +guides(colour=guide_legend(ncol=1, byrow = T)) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =14),
    legend.position = "",
    legend.spacing.x = unit(1.0, 'cm'),
    legend.key.height = unit(1, 'cm')
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "Land cover extent",
        y = "Estimated richness",
        subtitle = "A",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

cowplot::save_plot("plots/Fig_DK_effect_richness_landcover.tiff", effectplot, base_width = 8, base_height = 5, dpi = 800)

cowplot::save_plot("plots/Fig_DK_effect_richness_landcover.png", effectplot, base_width = 8, base_height = 5, dpi = 800)

###### biomass ##########
full_model <- lmer(log(totalBiomass_mg+1) ~ 
                     Urban_1000 +
                     Agriculture_50 + 
                     Open.uncultivated.land_1000+ # test if outlier drives the pattern
                     Wetland_50 +
                     Forest_500 +
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects)

# extract effects
gls1.alleffects <- allEffects(full_model)
plot(gls1.alleffects, 'Year', ylab="totalBiomass_mg (log+1)")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(full_model)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(full_model, ~ Agriculture_50 + Urban_1000 + Open.uncultivated.land_1000 + Forest_500 + Wetland_50 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

# make data frames for each land cover
temp <- effectdata$Agriculture_50
temp$landcover <- "Agriculture_50"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_50
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Open.uncultivated.land
temp <- effectdata$Open.uncultivated.land_1000
temp$landcover <- "Grassland_1000"
grass <- temp %>% 
  dplyr::rename(
    propcover = Open.uncultivated.land_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Wetland
temp <- effectdata$Wetland_50
temp$landcover <- "Wetland_50"
wet <- temp %>% 
  dplyr::rename(
    propcover = Wetland_50
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Forest
temp <- effectdata$Forest_500
temp$landcover <- "Forest_500"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_500
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
temp <- effectdata$`Time_band:cnumberTime`
#temp$landcover <- "Time_band:cnumberTime"
time <- temp

# combine the data frames
test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_50",
    "Grassland_1000",
    "Wetland_50",
    "Forest_500"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover",
      "Farmland cover",
      "Grassland cover",
      "Wetland cover",
      "Forest cover"
    )
  ) +guides(colour=guide_legend(nrow=1, byrow = T)) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =10),
    legend.position = "bottom",
    legend.spacing.x = unit(0.5, 'cm'),
    legend.key.height = unit(0.5, 'cm')
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "Land cover extent",
        y = "Estimated biomass (log+1)",
        subtitle = "B",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

cowplot::save_plot("plots/Fig_DK_effect_biomass_landcover.tiff", effectplot, base_width = 8, base_height = 5, dpi = 800)

cowplot::save_plot("plots/Fig_DK_effect_biomass_landcover.png", effectplot, base_width = 8, base_height = 5, dpi = 800)

##### Figure 5: time band ###############################
maxs <- c("Urban_1000", "Agriculture_1000", "Forest_1000")
facet_labs <- c("Midday", "Evening")
names(facet_labs) <- c("midday", "evening")

test <- allInsects_trans_landcovers[allInsects_trans_landcovers$Time_band=="midday",]
min(test$cnumberTime)
max(test$cnumberTime)

midday_plot <- allInsects_trans_landcovers[allInsects_trans_landcovers$Time_band=="midday",] %>% filter(maxLand_use %in% maxs) %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Urban_1000",
    "Agriculture_1000",
    "Forest_1000"
  )
) %>% ggplot(aes((cnumberTime), log(Biomass+1), colour = maxLand_use)) + geom_point(size = 3) + geom_smooth(method=lm, alpha = 0.3, size =1.5, show.legend = F)+ scale_colour_manual(values = landuseCols[c(1,2,5)], labels = c(
  "Urban",
  "Farmland",
  "Forest"
)) + facet_grid(.~Time_band, labeller = labeller(Time_band = facet_labs)) + scale_fill_manual(values = c("darkgrey", "darkgrey")) + labs(x = "", y= "log(biomass +1) (mg)", colour = "Sampling time") + theme_minimal() + theme(axis.text = element_text(size = 12), strip.text.x = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = "bottom", axis.title.y = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=8)))

midday_plot <- midday_plot + scale_x_continuous(breaks = c(-77.5, 10, 87.5), labels = c("12.00", "13.30", "15.00")) + ylim(0,8)

test <- allInsects_trans_landcovers[allInsects_trans_landcovers$Time_band=="evening",]
min(test$cnumberTime)
max(test$cnumberTime)

evening_plot <- allInsects_trans_landcovers[allInsects$Time_band=="evening",] %>% filter(maxLand_use %in% maxs) %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Urban_1000",
    "Agriculture_1000",
    "Forest_1000"
  )
) %>% ggplot(aes((cnumberTime), log(Biomass+1), colour = maxLand_use)) + geom_point(size = 3) + geom_smooth(method=lm, alpha = 0.3, size =1.5, show.legend = F)+ scale_colour_manual(values = landuseCols[c(1,2,5)], labels = c(
  "Urban",
  "Farmland",
  "Forest"
)) + facet_grid(.~Time_band, labeller = labeller(Time_band = facet_labs)) + scale_fill_manual(values = c("darkgrey", "darkgrey")) + labs(x = "", y= "", colour = "Sampling time") + theme_minimal() + theme(axis.text = element_text(size = 12), strip.text.x = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = "bottom", axis.title.y = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=8)))

evening_plot <- evening_plot + scale_x_continuous(breaks = c(-74, 21, 95), labels = c("17.00", "18.30", "20.00"))+ ylim(0,8)

plot_row <- cowplot::plot_grid(midday_plot, evening_plot)

# now add the title
title <- ggdraw() + 
  draw_label(
    "A: Denmark",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

sampling_time <- cowplot::plot_grid(
  title, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

cowplot::save_plot("plots/Fig5_DK_sampling_time_maxcover.tiff", sampling_time, base_width = 10, base_height = 6, dpi = 800)

cowplot::save_plot("plots/Fig5_DK_sampling_time_maxcover.png", sampling_time, base_width = 10, base_height = 6, dpi = 800)

# biomass difference between midday and evening
allInsects.long <- allInsects_trans_landcovers %>% 
  dplyr::select(Biomass, Time_band, cnumberTime, Urban_1000, Agriculture_1000, Open.uncultivated.land_1000, Wetland_50, Forest_250) %>% pivot_longer(-c(Biomass, cnumberTime, Time_band), names_to = "landcover", values_to = "cover")

head(allInsects.long)

#colorset = c('Urban_1000'=landuseCols[1],'Agriculture_1000'=landuseCols[2],'Open.uncultivated.land_1000'=landuseCols[3],'Wetland_50'=landuseCols[4], 'Forest_250' =landuseCols[5])

sampling_time <- allInsects.long  %>% ggplot(aes(cnumberTime, log(Biomass+1), colour = Time_band)) + geom_point() + scale_colour_manual(values = c("lightgrey", "darkgrey"), labels = c(
  "Midday",
  "Evening"
)) + geom_smooth(method=lm, aes(fill = Time_band), alpha = 0.1, size =1.5, show.legend = F) + scale_fill_manual(values = c("lightgrey", "darkgrey"))  + labs(x = "Standardised time", y= "log(biomass +1) (mg)", colour = "Land cover", subtitle = "A: Denmark") + theme(axis.text = element_blank(), plot.subtitle = element_text(size = 20, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 8), legend.position = "bottom")


##### Test of land cover diffs##############################

Ztest <- function(beta1,se1,beta2,se2){
  myZ <- (beta1 - beta2)/sqrt(beta1^2 + beta2^2)
  pvalue = 2*pnorm(abs(myZ), lower.tail = F)
  return(pvalue)
}

mySummary <-  summary(full_model)$coefficients

#Difference between farmland and open semi-natural habitats
# we predicted (1) to find lower biomass in agricultural areas compared to open semi-natural habitats (wetland and grassland) due to agricultural practices such as pesticide use, lower habitat complexity and increased human disturbance
Ztest(mySummary["Agriculture_1000","Estimate"],mySummary["Agriculture_1000","Std. Error"],
      mySummary["Open.uncultivated.land_1000","Estimate"],mySummary["Open.uncultivated.land_1000","Std. Error"])

Ztest(mySummary["Agriculture_1000","Estimate"],mySummary["Agriculture_1000","Std. Error"],
      mySummary["Wetland_50","Estimate"],mySummary["Wetland_50","Std. Error"])

#Difference between Urban and Open uncultivated
# we predicted (2) that urban cover would have the lowest biomass among all land covers due to the high proportion of impervious surfaces and low proportion of blue and green space.
Ztest(mySummary["Urban_1000","Estimate"],mySummary["Urban_1000","Std. Error"],
      mySummary["Open.uncultivated.land_1000","Estimate"],mySummary["Open.uncultivated.land_1000","Std. Error"])

Ztest(mySummary["Urban_1000","Estimate"],mySummary["Urban_1000","Std. Error"],
      mySummary["Wetland_50","Estimate"],mySummary["Wetland_50","Std. Error"])

Ztest(mySummary["Urban_1000","Estimate"],mySummary["Urban_1000","Std. Error"],
      mySummary["Forest_250","Estimate"],mySummary["Forest_250","Std. Error"])

Ztest(mySummary["Urban_1000","Estimate"],mySummary["Urban_1000","Std. Error"],
      mySummary["Agriculture_1000","Estimate"],mySummary["Agriculture_1000","Std. Error"])

##### DK biomass predictions ##############################

library(lme4)
library(lmerTest)

# urban
lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000 + 
                  Time_band + 
                  Time_band:cnumberTime +
                  log(cStops+1) + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), 
                data=allInsects_trans_landcovers)
summary(lme1000)

newData = data.frame(Urban_1000=0.5,
                     cStops=0,
                     cyDay = 0,
                     Time_band = "midday",
                     cnumberTime = 0)

predFun <- function(fit) {
  predict(fit,newData,re.form=NA)
}

#make predictions
Urban1 <- t(as_tibble(exp(predict(lme1000,newdata=newData,re.form=NA))))

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
urb <- bb[["data"]]
Urban2 <- t(as_tibble(exp(quantile(bb$t,c(0.025,0.975)))))

Urban <- cbind(Urban1, Urban2)
Urban <- as.data.frame(Urban)
colnames(Urban)
names(Urban)[1] <- "predBiomass"
names(Urban)[2] <- "lowCI"
names(Urban)[3] <- "highCI"
row.names(Urban) <- "Urban"

# farmland
lme1000 <- lmer(log(Biomass+1) ~ 
                  Agriculture_1000 + 
                  Time_band + 
                  Time_band:cnumberTime +
                  log(cStops+1) + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), 
                data=allInsects_trans_landcovers)
summary(lme1000)
newData = data.frame(Agriculture_1000=0.5,
                     cStops=0,
                     cyDay = 0,
                     Time_band = "midday",
                     cnumberTime = 0)


#make predictions
Farmland1 <- t(as_tibble(exp(predict(lme1000,newdata=newData,re.form=NA))))

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
farm <- bb[["data"]]
Farmland2 <- t(as_tibble(exp(quantile(bb$t,c(0.025,0.975)))))

Farmland <- cbind(Farmland1, Farmland2)
Farmland <- as.data.frame(Farmland)
colnames(Farmland)
names(Farmland)[1] <- "predBiomass"
names(Farmland)[2] <- "lowCI"
names(Farmland)[3] <- "highCI"
row.names(Farmland) <- "Farmland"

# Grassland
lme1000 <- lmer(log(Biomass+1) ~ 
                  Open.uncultivated.land_1000 + 
                  Time_band + 
                  Time_band:cnumberTime +
                  log(cStops+1) + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), 
                data=allInsects_trans_landcovers)
summary(lme1000)
newData = data.frame(Open.uncultivated.land_1000=0.5,
                     cStops=0,
                     cyDay = 0,
                     Time_band = "midday",
                     cnumberTime = 0)

#make predictions
Grassland1 <- t(as_tibble(exp(predict(lme1000,newdata=newData,re.form=NA))))

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
grass <- bb[["data"]]
Grassland2 <- t(as_tibble(exp(quantile(bb$t,c(0.025,0.975)))))

Grassland <- cbind(Grassland1, Grassland2)
Grassland <- as.data.frame(Grassland)
colnames(Grassland)
names(Grassland)[1] <- "predBiomass"
names(Grassland)[2] <- "lowCI"
names(Grassland)[3] <- "highCI"
row.names(Grassland) <- "Grassland"

# wetland
lme1000 <- lmer(log(Biomass+1) ~ 
                  Wetland_1000 + 
                  Time_band + 
                  Time_band:cnumberTime +
                  log(cStops+1) + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), 
                data=allInsects_trans_landcovers)
summary(lme1000)

newData = data.frame(Wetland_1000=0.5,
                     cStops=0,
                     cyDay = 0,
                     Time_band = "midday",
                     cnumberTime = 0)

#make predictions
wetland1 <- t(as_tibble(exp(predict(lme1000,newdata=newData,re.form=NA))))

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
wet <- bb[["data"]]
wetland2 <- t(as_tibble(exp(quantile(bb$t,c(0.025,0.975)))))

Wetland <- cbind(wetland1, wetland2)
Wetland <- as.data.frame(Wetland)
colnames(Wetland)
names(Wetland)[1] <- "predBiomass"
names(Wetland)[2] <- "lowCI"
names(Wetland)[3] <- "highCI"
row.names(Wetland) <- "Wetland"

# Forest
lme1000 <- lmer(log(Biomass+1) ~ 
                  Forest_1000 + 
                  Time_band + 
                  Time_band:cnumberTime +
                  log(cStops+1) + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), 
                data=allInsects_trans_landcovers)
summary(lme1000)

newData = data.frame(Forest_1000=0.5,
                     cStops=0,
                     cyDay = 0,
                     Time_band = "midday",
                     cnumberTime = 0)


#make predictions
Forest1 <- t(as_tibble(exp(predict(lme1000,newdata=newData,re.form=NA))))

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
fors <- bb[["data"]]
Forest2 <- t(as_tibble(exp(quantile(bb$t,c(0.025,0.975)))))

Forest <- cbind(Forest1, Forest2)
Forest <- as.data.frame(Forest)
colnames(Forest)
names(Forest)[1] <- "predBiomass"
names(Forest)[2] <- "lowCI"
names(Forest)[3] <- "highCI"
row.names(Forest) <- "Forest"

predConfData <- rbind(Urban, Farmland)
predConfData <- rbind(predConfData, Grassland)
predConfData <- rbind(predConfData, Wetland)
predConfData <- rbind(predConfData, Forest)
#predConfData <- rbind(predConfData, Unspecified)

predConfData <- rownames_to_column(predConfData, var = "landcover")

# plot
p <- predConfData %>%
  mutate(landcover = fct_relevel(
    landcover,
    "Urban",
    "Farmland",
    "Grassland",
    "Wetland",
    "Forest"
  )) %>% ggplot(aes(landcover, predBiomass, colour = landcover))

finalplot <-
  p + geom_pointrange(aes(ymin = lowCI, ymax = highCI),
                      size = 1.5,
                      show.legend = F) + scale_colour_manual(values = landuseCols) + theme_minimal_grid() + theme(legend.title = element_blank(),
                                                                                                                  legend.key = element_rect(size = 0.1),
                                                                                                                  legend.key.size = unit(1, 'cm')
                      ) + labs(x = "\nLand cover", y = "Predicted biomass (mg) and 95% CIs\n", subtitle = "B") + theme(plot.subtitle = element_text(size = 20, face = "bold")) + scale_y_log10()

cowplot::save_plot("plots/DK_predicted_biomass.png", finalplot, base_width = 8, base_height = 5)

##### combining the predicted data to re-run model and calculate effects ####
predeffect <- merge(urb, farm)
predeffect <- merge(predeffect, grass)
predeffect <- merge(predeffect, wet)
predeffect <- merge(predeffect, fors)

predeffect <- predeffect %>% rename(Biomass = `log(Biomass + 1)`) # be mindful that biomass is +1 and logtransformed here, the same for stops
predeffect <- predeffect %>% rename(cStops = `log(cStops + 1)`) 

# run model
library(nlme)
gls1 <- lme(Biomass ~ Agriculture_1000 + 
              Urban_1000 +
              Open.uncultivated.land_1000 +
              Forest_1000 +
              Wetland_1000 +
              Time_band +
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            data=predeffect)

summary(gls1) # urban non-significant
exp(1.867667-1) # farmland: 2.381349
exp(3.345569-1) # grassland: 10.43921
exp(2.202820-1) # forest: 3.329493
exp(2.369767-1) # wetland: 3.934434
exp(0.319483-1) # evening: 0.5063551

library(effects)
gls1.alleffects <- allEffects(gls1)
#plot(gls1.alleffects, 'Urban_1000', ylab="Biomass")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(gls1)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(gls1, ~ Agriculture_1000 + Urban_1000 + Open.uncultivated.land_1000 + Forest_1000 + Wetland_1000 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

##### predicted effect plot ####
temp <- effectdata$Agriculture_1000
temp$landcover <- "Agriculture_1000"
farm <- temp %>% 
  rename(
    propcover = Agriculture_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Open.uncultivated.land
temp <- effectdata$Open.uncultivated.land_1000
temp$landcover <- "Grassland_1000"
grass <- temp %>% 
  rename(
    propcover = Open.uncultivated.land_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Wetland
temp <- effectdata$Wetland_1000
temp$landcover <- "Wetland_1000"
wet <- temp %>% 
  rename(
    propcover = Wetland_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Forest
temp <- effectdata$Forest_1000
temp$landcover <- "Forest_1000"
forest <- temp %>% 
  rename(
    propcover = Forest_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_1000",
    "Grassland_1000",
    "Wetland_1000",
    "Forest_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover",
      "Farmland cover",
      "Grassland cover",
      "Wetland cover",
      "Forest cover"
    )
  ) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom"
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x, "%"))  + scale_y_continuous(
        limits = c(2.5, 7.5),
        labels = function(x)
          paste0(((x-1)) * 1, "%")) + geom_ribbon( # -1 since +1 is used in biomass response, but exp() also? - then it introduced a lot of decimals, but gives the right effect size
            aes(
              ymin = fit-se,
              ymax = fit+se,
              group = landcover
            ),
            linetype = 2,
            alpha = 0.2,
            show.legend = F
          ) + labs(
            x = "Land cover",
            y = "Predicted effect change on biomass",
            subtitle = "B",
            colour = "Land cover type"
          ) + scale_fill_manual(values = landuseCols)

cowplot::save_plot("plots/DK_predictedeffect_landcover.tiff", effectplot, base_width = 10, base_height = 6, dpi = 1200)

##### Test timeband interactions  ################################
# farmland
lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000 +
                  Time_band*Agriculture_1000 + 
                  Open.uncultivated.land_1000 +
                  Forest_250 +
                  Wetland_50 +
                  Time_band + 
                  cyDay + 
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_trans_landcovers)
summary(lme1000)

# grassland
lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000 +
                  Agriculture_1000 + 
                  Open.uncultivated.land_1000*Time_band +
                  Forest_250 +
                  Wetland_50 +
                  Time_band + 
                  cyDay + 
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_trans_landcovers)
summary(lme1000)

# urban
lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000*Time_band +
                  Agriculture_1000 + 
                  Open.uncultivated.land_1000 +
                  Forest_250 +
                  Wetland_50 +
                  Time_band + 
                  cyDay + 
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_trans_landcovers)
summary(lme1000)

# forest
lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000 +
                  Agriculture_1000 + 
                  Open.uncultivated.land_1000 +
                  Forest_250*Time_band +
                  Wetland_50 +
                  Time_band + 
                  cyDay + 
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_trans_landcovers)
summary(lme1000)

# wetland
lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000 +
                  Agriculture_1000 + 
                  Open.uncultivated.land_1000 +
                  Forest_250 +
                  Wetland_50*Time_band +
                  Time_band + 
                  cyDay + 
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_trans_landcovers)
summary(lme1000)

### Germany ##############################################

##### Figure 3: the data ##########################

qU <- ggplot(allInsects,aes(x=Urban_1000,y=(Biomass+1)))+
  geom_point(col=landuseCols[1])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70") + scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Urban cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qF <- ggplot(allInsects,aes(x=Agriculture_1000,y=(Biomass+1)))+
  geom_point(col=landuseCols[2])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Farmland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qD <- ggplot(allInsects,aes(x=Open.uncultivated_1000,y=(Biomass+1)))+
  geom_point(col=landuseCols[3])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    #limits = c(0, 0.16),
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("Biomass (mg)") + labs(subtitle = "Grassland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 24))

qW <- ggplot(allInsects,aes(x=Wetland_1000,y=(Biomass+1)))+
  geom_point(col=landuseCols[4])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Wetland cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

qFo <- ggplot(allInsects,aes(x=Forest_1000,y=(Biomass+1)))+
  geom_point(col=landuseCols[5])+
  scale_y_log10() +
  theme_bw() +
  geom_smooth(method="lm",color="grey70")+scale_x_continuous(
    labels = function(x)
      paste0(x, "%")) +
  xlab("") +ylab("") + labs(subtitle = "Forest cover") + theme(plot.subtitle = element_text(size = 12, face = "bold"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text = element_text(size = 16), axis.title.y = element_text(size = 16))

fig3 <- cowplot::plot_grid(qU,qF,qD,qW,qFo,ncol=1)

ggsave("plots/Landcover_percent.png",width=4,height=12)
cowplot::save_plot("plots/DE_Landcover_percent.tiff", fig3, base_width = 4, base_height = 12, dpi = 1200)


##### Figure 3: DE pie chart #####################################

library(dplyr)
library(tidyr)

routeMeans <- allInsects %>% 
  group_by(RouteID) %>%
  summarise(meanBiomass = mean(Biomass))

allInsects <- inner_join(allInsects,routeMeans,by="RouteID")

#remove duplicates
allInsects <- allInsects %>%
  select(RouteID,meanBiomass,Agriculture_1000,
         Forest_1000,Open.uncultivated_1000,
         Urban_1000,Wetland_1000) %>%
  distinct()

#fill in missing column
allInsects$totalLand <- apply(allInsects[,3:7],1,sum)
allInsects$Other_1000 <- 100-allInsects$totalLand

#divide up biomass into quantiles
allInsects$BiomassCats <- cut_number(allInsects$meanBiomass,n=5)

#mean land cover per biomass cats
allInsects_cat <- allInsects %>%
  group_by(BiomassCats) %>%
  summarise(Agriculture_1000 = mean(Agriculture_1000),
            Forest_1000 = mean(Forest_1000),
            Open.uncultivated_1000 = mean(Open.uncultivated_1000),
            Urban_1000 = mean(Urban_1000),
            Wetland_1000 = mean(Wetland_1000),
            Other_1000 = mean(Other_1000))

#plot data by biomass categories
allInsects_melt <- gather(allInsects_cat, Land_cover, value, -BiomassCats)

#plot
library(forcats)
allInsects_melt <- allInsects_melt %>% mutate(
  Land_cover = fct_relevel(
    Land_cover,
    "Urban_1000",
    "Agriculture_1000",
    "Open.uncultivated_1000",
    "Wetland_1000",
    "Forest_1000",
    "Other_1000"))  

levels(allInsects_melt$BiomassCats)

biomass.labs <- c("[0,46]"=" 0-46 mg", "(46,115]"="486-115 mg", "(115,209]"="115-209 mg", 
                  "(209,502]"="209-502 mg", "(502,1.38e+03]"="502-1380 mg")

forplot <- allInsects_melt %>% filter(BiomassCats %in% c("[0,46]", "(502,1.38e+03]"))

fig3_2 <- ggplot(forplot,aes(x="",y=value,fill=Land_cover, order = Land_cover))+
  geom_bar(stat="identity")+
  facet_wrap(~BiomassCats, labeller = labeller(BiomassCats=biomass.labs[c(1,5)]), ncol = 1)+
  coord_polar("y")+
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position="left", strip.text = element_text(size = 16), legend.text = element_text(size = 16), legend.key.size = unit(1.5, 'cm')) + scale_fill_manual(values = landuseCols, name = "", labels = c("Urban", "Farmland", "Grassland", "Wetland", "Forest", "Other")) 

cowplot::save_plot("plots/DE_landcover_biomass_proportions.tiff", fig3_2, base_width = 12, base_height = 6, dpi = 800)

##### Table 1: Linear Mixed Effects Model - Land covers #################

#full model and final
lme1000 <- full_model <- lmer(log(Biomass+1) ~ 
                                Agriculture_1000 + 
                                Urban_1000 +
                                Open.uncultivated_1000+
                                Wetland_1000 +
                                Forest_250 +
                                Time_band + 
                                Time_band:cnumberTime + 
                                cStops + 
                                cyDay + 
                                (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)
vif(lme1000)
tab_model(lme1000)
r.squaredGLMM(lme1000)
#R2m       R2c
#[1,] 0.3593968 0.8313989
confint(lme1000)

##### sqrt land covers ######
allInsects$sqrt_Agriculture_1000 <- sqrt(allInsects$Agriculture_1000)
allInsects$sqrt_Urban_1000 <- sqrt(allInsects$Urban_1000)
allInsects$sqrt_Open.uncultivated_1000 <- sqrt(allInsects$Open.uncultivated_1000)
allInsects$sqrt_Wetland_1000 <- sqrt(allInsects$Wetland_1000)
allInsects$sqrt_Forest_250 <- sqrt(allInsects$Forest_250)

lme1000 <- lmer(log(Biomass+1) ~ 
                  sqrt_Agriculture_1000 + 
                  sqrt_Urban_1000 +
                  sqrt_Open.uncultivated_1000 +
                  sqrt_Wetland_1000 +
                  sqrt_Forest_250 +
                  Time_band + 
                  Time_band:cnumberTime + 
                  cStops + 
                  cyDay + 
                  (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)
vif(lme1000)
#now agriculture is most important... but urban comes out on top with AIC
r.squaredGLMM(lme1000)
#R2m       R2c
#[1,] 0.3582978 0.8328536

#with model simplification
lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000 +
                  Time_band + 
                  Time_band:cnumberTime + 
                  (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)

#as gls
lme1000 <- lme(log(Biomass+1) ~ Agriculture_1000 + 
                 Urban_1000 +
                 Open.uncultivated_1000+
                 Wetland_1000 +
                 Forest_250 +
                 Time_band + 
                 Time_band:cnumberTime + 
                 cStops + 
                 cyDay,
               random=~1|PilotID,
               correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
               data=allInsects)
summary(lme1000)
vif(lme1000)
r.squaredGLMM(lme1000)

##### Table 1: multicomp land covers - simple model ##########################################

# pairwise comparison to farmland
pair.ht <- glht(lme1000, linfct = c(" Agriculture_1000 - Forest_250 = 0", 
                                    "Agriculture_1000 - Wetland_1000 = 0", 
                                    "Agriculture_1000 - Open.uncultivated_1000 = 0"))
summary(pair.ht) 

confint(pair.ht)

# pairwise comparison to urban
pair.ht <- glht(lme1000, linfct = c(
  "Urban_1000 - Forest_250  = 0",
  "Urban_1000 - Wetland_1000 = 0", 
  "Urban_1000 - Open.uncultivated_1000 = 0",
  "Urban_1000 - Agriculture_1000 = 0"))

summary(pair.ht) 

confint(pair.ht)

#compared to grassland
pair.ht <- glht(lme1000, linfct = c(
  "Open.uncultivated_1000 - Wetland_1000 = 0",
  "Open.uncultivated_1000 - Forest_250  = 0"))

summary(pair.ht) 

confint(pair.ht)


#compared to grassland
pair.ht <- glht(lme1000, linfct = c(
  "Wetland_1000 - Forest_250  = 0"))

summary(pair.ht) 

confint(pair.ht)

#squareroot
pair.ht <- glht(lme1000, linfct = c("sqrt_Forest_250 - sqrt_Agriculture_1000 = 0", "sqrt_Wetland_1000 - sqrt_Agriculture_1000 = 0", "sqrt_Open.uncultivated_1000 - sqrt_Agriculture_1000 = 0"))
confint(pair.ht)
pair.ht <- glht(lme1000, linfct = c("sqrt_Forest_250 - sqrt_Urban_1000 = 0","sqrt_Wetland_1000 - sqrt_Urban_1000 = 0", "sqrt_Open.uncultivated_1000 - sqrt_Urban_1000 = 0","sqrt_Agriculture_1000  - sqrt_Urban_1000  = 0"))
confint(pair.ht)
pair.ht <- glht(lme1000, linfct = c("sqrt_Wetland_1000 - sqrt_Open.uncultivated_1000 = 0","sqrt_Forest_250 - sqrt_Open.uncultivated_1000  = 0"))
confint(pair.ht)

##### Table 1: coefficient plot ######################################################
#point estimate

myConfint <- confint(full_model)

coefDF <- data.frame(Landcover = names(fixef(full_model)),
                     estimate = as.numeric(fixef(full_model)),
                     lowerCI = as.numeric(myConfint[-c(1:3),1]),
                     upperCI = as.numeric(myConfint[-c(1:3),2]))

#subset to land covers
coefDF <- subset(coefDF, Landcover %in% c("Agriculture_1000","Urban_1000","Open.uncultivated_1000","Wetland_1000","Forest_250"))

p <-coefDF %>% 
  mutate(Landcover = fct_relevel(
    Landcover,
    "Urban_1000",
    "Agriculture_1000",
    "Open.uncultivated_1000",
    "Wetland_1000",
    "Forest_250"
  )) %>% 
  ggplot(aes(Landcover, estimate))

c<-p+scale_size_area(max_size = 1.5)

#this puts in a dotted line at the point of group difference
d<-c+geom_hline(yintercept = 0, lty=2,size=1)

coefPlot <-
  d + geom_pointrange(aes(ymin = lowerCI, ymax = upperCI, colour = Landcover),
                      size = 1.5) + theme_minimal_grid() + theme(
                        legend.title = element_blank(),
                        legend.key = element_rect(size = 0.1),
                        legend.key.size = unit(1, 'cm')
                      ) + labs(x = "", y = "Effect on flying insect biomass & 95% CIs") + scale_x_discrete(
                        labels = c(
                          "Urban_1000" = "Urban",
                          "Agriculture_1000" = "Farmland",
                          "Open.uncultivated_1000" = "Grassland",
                          "Wetland_1000" = "Wetland",
                          "Forest_250" = "Forest"
                        )
                      ) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text.x = element_text(size = 10, angle = 90, vjust=0.3,hjust=1), legend.position = "none", axis.title.y = element_text(size = 10)) + scale_colour_manual(values = landuseCols) + coord_fixed(50)

cowplot::save_plot("plots/DE_coefficient_plot.tiff", coefPlot, base_width = 4, base_height = 5, dpi = 300)

##### ilr trans of land covers ################################

#Approach 1 using the complmrob package
library(complmrob)
#The variables on the right-hand-side of the 
#formula are transformed with the isometric log-ratio transformation (isomLR) and a robust linear regression model is fit to those transformed variables.

#The log-ratio methodology studies this relative relations using the quotient between parts, to be more precise, using the logarithm between ratios of parts

#need to add small values to zeros in the covariate data to get the function to work
myvars <-allInsects[,c("Agriculture_1000","Urban_1000","Open.uncultivated_1000","Wetland_1000","Forest_250")]
myvars<- data.frame(apply(myvars,c(1,2),
                          function(x)ifelse(x==0,0.001,x)))
names(myvars) <- c("cAgriculture_1000","cUrban_1000","cOpen.uncultivated_1000","cWetland_1000","cForest_250")
allInsects <- cbind(allInsects,myvars)

#original function from complmrob package (no random effects)
complm1 <- complmrob(log(Biomass+1) ~ cAgriculture_1000 + cUrban_1000 + cOpen.uncultivated_1000 + cWetland_1000 + cForest_250, data = allInsects)
summary(complm1)

#fit model as random effects model - setting urban as first composition - other land uses relative to this
landUses <- c("cAgriculture_1000","cUrban_1000", "cOpen.uncultivated_1000","cWetland_1000","cForest_250")
tmpPred <- data.frame(isomLR(as.matrix(allInsects[,landUses]),2))
tmpNames <- paste(names(tmpPred),collapse="+")
tmpCovariates <- allInsects[,c("Time_band","cnumberTime","cStops","cyDay","RouteID","PilotID")]
tmp <- cbind(tmpPred,tmpCovariates)

lme1000 <- lmer(as.formula(paste("log(Biomass+1) ~ Time_band + Time_band:cnumberTime + cStops + cyDay + 
                  (1|RouteID) + (1|PilotID) +",tmpNames)), data=allInsects)
summary(lme1000)
car::vif(lme1000)

#write function to do the same for all land uses
landUses <- c("cAgriculture_1000","cUrban_1000", "cOpen.uncultivated_1000","cWetland_1000","cForest_250")

fitISOLRmodel <- function(component=1){
  #get isoLR for specified component
  tmpPred <- data.frame(isomLR(as.matrix(allInsects[,landUses]),
                               component))
  #paste into a model formula
  tmpNames <- paste(names(tmpPred),collapse="+")
  #get other variables that we will need
  tmpCovariates <- allInsects[,c("Time_band","cnumberTime","cStops",
                                 "cyDay","RouteID","PilotID")]
  tmp <- cbind(tmpPred,tmpCovariates)
  #fit model
  require(lme4)
  require(lmerTest)
  lme1000 <- lmer(as.formula(paste("log(Biomass+1) ~ Time_band + 
                                   Time_band:cnumberTime + cStops + cyDay +
                                   (1|RouteID) + (1|PilotID) +",
                                   tmpNames)), data=allInsects)
  #extract data for land use variable that we are interesed in
  temp <- data.frame(summary(lme1000)$coefficients)
  temp <- temp[row.names(temp)==landUses[component],]
  return(temp)
  
}

#fit function to each land use
1:5 %>%
  purrr::map_df(~fitISOLRmodel(.))
#                             Estimate  Std..Error       df    t.value    Pr...t..
#cAgriculture_1000        0.007445503 0.005590892 57.13133  1.3317200 0.188241143
#cUrban_1000             -0.047679661 0.015818573 57.43560 -3.0141568 0.003831099
#cOpen.uncultivated_1000 -0.002310477 0.010498392 59.91710 -0.2200792 0.826557638
#cWetland_1000           -0.015208667 0.040177716 55.48551 -0.3785349 0.706477826
#cForest_250              0.013697219 0.011585161 46.42869  1.1823072 0.243103590

##### ilr comparison ##################################################

library(robCompositions)
?lmCoDaX

#transformation function
?pivotCoord
transVars <- pivotCoord(myvars,2)
head(transVars)

lmCoDaX(log(allInsects$Biomass+1), myvars, method="classical")
lmCoDaX(log(allInsects$Biomass+1), myvars, method="robust")

#compare isomLE with pivotCoord
landUses <- c("cAgriculture_1000","cUrban_1000", "cOpen.uncultivated_1000","cWetland_1000","cForest_250")
tmpPred <- data.frame(isomLR(as.matrix(allInsects[,landUses]),1))
tmpPred_pc <- pivotCoord(allInsects[,landUses],1)
qplot(tmpPred$cAgriculture_1000,tmpPred_pc$`cAgriculture_1000_cU-cO-cW-cF`)
qplot(tmpPred$cUrban_1000,tmpPred_pc$`cUrban_1000_cO-cW-cF`)
#they are the same...

##### AIC check ##############################################

library(MuMIn)
options(na.action = "na.fail")
dd <- dredge(lme1000)
subset(dd, delta < 2)

# Visualize the model selection table:
par(mar = c(3,5,6,4))
plot(dd, labAsExpr = TRUE)

# Model average models with delta AICc < 4
model.avg(dd, subset = delta < 2)

#only urban remains in best set

lme1000 <- lmer(log(Biomass+1) ~ 
                  Urban_1000 +
                  Time_band + 
                  (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)

##### Figure 4: effect plot ##########################

predEffects <- predictorEffects(lme1000)

getEffects <- function(predEffects,var){
  temp <- predEffects[var]
  data.frame(fit = temp[[1]]$fit,
             se = temp[[1]]$se,
             lower = temp[[1]]$lower,
             upper = temp[[1]]$upper,
             propcover = temp[[1]]$x[,1],
             landcover = rep(var,length(temp[[1]]$fit)))
}

#apply to each land cover
urb <- getEffects(predEffects,var="Urban_1000")
farm <- getEffects(predEffects,var="Agriculture_1000")
grass <- getEffects(predEffects,var="Open.uncultivated_1000")
wet <- getEffects(predEffects,var="Wetland_1000")
forest <- getEffects(predEffects,var="Forest_250")
test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot <- test %>% mutate(
  landcover = fct_relevel(landcover,
                          "Urban_1000",
                          "Agriculture_1000",
                          "Open.uncultivated_1000",
                          "Wetland_1000",
                          "Forest_250")) %>% 
  ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover",
      "Farmland cover",
      "Grassland cover",
      "Wetland cover",
      "Forest cover")) +
  guides(colour=guide_legend(ncol=1, byrow = T)) + 
  theme_minimal_grid() + 
  theme(plot.subtitle = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size =14),
        legend.position = "right",
        legend.spacing.x = unit(1.0, 'cm'),
        legend.key.height = unit(1, 'cm')) + 
  scale_x_continuous(limits = c(0, 100),
                     labels = function(x)
                       paste0(x, "%")) + scale_y_continuous(limits = c(1.5, 7)) + geom_ribbon(
                         aes(
                           ymin = fit-se,
                           ymax = fit+se,
                           group = landcover
                         ),
                         linetype = 2,
                         alpha = 0.2,
                         show.legend = F) + 
  labs(
    x = "Land cover extent",
    y = "log Predicted biomass (mg)",
    subtitle = "B: Germany",
    colour = "Land cover type") + 
  scale_fill_manual(values = landuseCols)

#cowplot::save_plot("plots/Fig4_DE_effect_landcover.tiff", effectplot, base_width = 12, base_height = 8, dpi = 800)

ggsave("plots/Fig4_DE_effect_landcover.png", width = 8, height = 5)

##### Test of land cover diffs##############################

Ztest <- function(beta1,se1,beta2,se2){
  myZ <- (beta1 - beta2)/sqrt(beta1^2 + beta2^2)
  pvalue = 2*pnorm(abs(myZ), lower.tail = F)
  return(pvalue)
}

mySummary <-  summary(lme1000)$coefficients

#Difference between Agriculture and Open uncultivated
Ztest(mySummary["Agriculture_1000","Estimate"],mySummary["Agriculture_1000","Std. Error"],
      mySummary["Open.uncultivated_1000","Estimate"],mySummary["Open.uncultivated_1000","Std. Error"])

#Difference between Urban and Open uncultivated
Ztest(mySummary["Urban_1000","Estimate"],mySummary["Urban_1000","Std. Error"],
      mySummary["Open.uncultivated_1000","Estimate"],mySummary["Open.uncultivated_1000","Std. Error"])

##### PCA axes as variables##################################

#run script in script 07 to get PCA axes variables

lme1000 <- lmer(log(Biomass+1) ~ 
                  Urbanization_gradient +
                  Forest_gradient +
                  Time_band + 
                  Time_band:cnumberTime + 
                  cStops + 
                  cyDay + 
                  (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)
vif(lme1000)
r.squaredGLMM(lme1000)

#Fixed effects:
#  Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)                   4.754e+00  2.604e-01  2.657e+01  18.259  < 2e-16 ***
#  Urbanization_gradient        -4.042e-01  1.881e-01  5.966e+01  -2.150  0.03566 *  
#  Forest_gradient              -2.180e-01  1.672e-01  5.724e+01  -1.304  0.19753    
#Time_bandevening              3.808e-01  1.152e-01  6.398e+01   3.305  0.00156 ** 
#cStops                       -3.212e-01  1.905e-01  5.009e+01  -1.686  0.09807 .  
#cyDay                        -6.489e-02  4.313e-02  5.681e+01  -1.504  0.13800    
#Time_bandmidday:cnumberTime  -4.606e-05  1.805e-03  8.419e+01  -0.026  0.97970    
#Time_bandevening:cnumberTime  5.106e-03  2.062e-03  8.362e+01   2.476  0.01531 *

lme1000 <- lmer(log(Biomass+1) ~ 
                  Urbanization_gradient * Time_band +
                  Forest_gradient * Time_band +
                  Time_band + 
                  Time_band:cnumberTime + 
                  cStops + 
                  cyDay + 
                  (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)
#Urbanization_gradient:Time_bandevening  0.2234471  0.1122381 61.6377367   1.991  0.05094


##### DE biomass predictions% ##############################

library(lme4)
library(lmerTest)

lme1000 <- lmer(log(Biomass+1) ~ 
                  sqrt(Forest_1000) + 
                  Time_band + 
                  Time_band:cnumberTime +
                  log(tr_signals+1) + cyDay + 
                  (1|RouteID) + (1|PilotID), 
                data=allInsects)
summary(lme1000)
newData = data.frame(Forest_1000=0.5,
                     tr_signals=0,
                     cyDay = 0,
                     Time_band = "midday",
                     cnumberTime = 0)


#make predictions
exp(predict(lme1000,newdata=newData,re.form=NA))

predFun <- function(fit) {
  predict(fit,newData,re.form=NA)
}

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
exp(quantile(bb$t,c(0.025,0.975)))

##### Test timeband interactions ################################

lme1000 <- lmer(log(Biomass+1) ~ 
                  Time_band*Agriculture_1000 + 
                  Time_band + 
                  Time_band:cnumberTime + cStops + cyDay + 
                  (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)

lme1000 <- lmer(log(Biomass+1) ~ 
                  Time_band*Urban_1000 + 
                  Time_band + 
                  Time_band:cnumberTime + cStops + cyDay + 
                  (1|RouteID) + (1|PilotID), data=allInsects)
summary(lme1000)


ggplot(allInsects,aes(x=as.numeric(Time_band), y=log(Biomass+1)))+
  geom_jitter(aes(colour=Land_use))+
  stat_smooth(aes(colour=Land_use),method="lm",se=FALSE)


##### Figure 5: time band #####################################################

#identify maxLanduse at each route
allInsects$maxLand_use <- apply(allInsects,1,function(x)which.max(x[c("Agriculture_1000","Forest_1000",
                                                                      "Open.uncultivated_1000","Wetland_1000",
                                                                      "Urban_1000")]))
allInsects$maxLand_use[allInsects$maxLand_use==1]<- "Agriculture_1000"
allInsects$maxLand_use[allInsects$maxLand_use==2]<- "Forest_1000"
allInsects$maxLand_use[allInsects$maxLand_use==3]<- "Open.uncultivated_1000"
allInsects$maxLand_use[allInsects$maxLand_use==4]<- "Wetland_1000"
allInsects$maxLand_use[allInsects$maxLand_use==5]<- "Urban_1000"


maxs <- c("Urban_1000", "Agriculture_1000", "Forest_1000")
facet_labs <- c("Midday", "Evening")
names(facet_labs) <- c("midday", "evening")

test <- allInsects[allInsects$Time_band=="midday",]
min(test$cnumberTime)
max(test$cnumberTime)

temp <- unique(cbind(test[,c("StartTime","cnumberTime")]))
temp[order(temp$cnumberTime),]

midday_plot <- allInsects[allInsects$Time_band=="midday",] %>% filter(maxLand_use %in% maxs) %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Urban_1000",
    "Agriculture_1000",
    "Forest_1000"
  )
) %>% ggplot(aes((cnumberTime), log(Biomass+1), colour = maxLand_use)) + geom_point(size = 3) + geom_smooth(method=lm, alpha = 0.3, size =1.5, show.legend = F)+ scale_colour_manual(values = landuseCols[c(1,2,5)], labels = c(
  "Urban",
  "Farmland",
  "Forest"
)) + facet_grid(.~Time_band, labeller = labeller(Time_band = facet_labs)) + scale_fill_manual(values = c("darkgrey", "darkgrey")) + labs(x = "", y= "log(biomass +1) (mg)", colour = "Sampling time") + theme_minimal() + theme(axis.text = element_text(size = 12), strip.text.x = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = "bottom", axis.title.y = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=8)))

midday_plot <- midday_plot + 
  scale_x_continuous(breaks = c(-90, 0, 90), 
                     labels = c("12.00", "13.30", "15.00")) + ylim(0,8)

test <- allInsects[allInsects$Time_band=="evening",]
min(test$cnumberTime)
max(test$cnumberTime)

temp <- unique(cbind(test[,c("StartTime","cnumberTime")]))
temp[order(temp$cnumberTime),]

evening_plot <- allInsects[allInsects$Time_band=="evening",] %>% filter(maxLand_use %in% maxs) %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Urban_1000",
    "Agriculture_1000",
    "Forest_1000"
  )
) %>% ggplot(aes((cnumberTime), log(Biomass+1), colour = maxLand_use)) + geom_point(size = 3) + geom_smooth(method=lm, alpha = 0.3, size =1.5, show.legend = F)+ scale_colour_manual(values = landuseCols[c(1,2,5)], labels = c(
  "Urban",
  "Farmland",
  "Forest"
)) + facet_grid(.~Time_band, labeller = labeller(Time_band = facet_labs)) + scale_fill_manual(values = c("darkgrey", "darkgrey")) + labs(x = "", y= "", colour = "Sampling time") + theme_minimal() + theme(axis.text = element_text(size = 12), strip.text.x = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 14), legend.position = "bottom", axis.title.y = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size=8)))

evening_plot <- evening_plot + 
  scale_x_continuous(breaks = c(-75, 15, 105), 
                     labels = c("17.00", "18.30", "20.00"))+ 
  ylim(0,8)

plot_row <- cowplot::plot_grid(midday_plot, 
                               evening_plot,ncol=2,nrow=1)

# now add the title
title <- ggdraw() + 
  draw_label(
    "B: Germany",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

sampling_time <- cowplot::plot_grid(
  title, plot_row,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave("plots/DE_sampling_time_maxcover.png", width = 10, height = 6)

##### Table S4.2: exclude urban cover ######################################################

summary(allInsects$Urban_1000)
#include only sites with less than 5% urban

allInsects_lowUrban <- subset(allInsects, Urban_1000 <5)

nrow(allInsects_lowUrban)

lme1000 <- lmer(log(Biomass+1) ~ 
                  Agriculture_1000 + 
                  Open.uncultivated_1000+
                  Wetland_1000 +
                  Forest_250 +
                  Time_band + 
                  Time_band:cnumberTime + 
                  cStops + 
                  cyDay + 
                  (1|RouteID) + (1|PilotID), data=allInsects_lowUrban)

# pairwise comparison to farmland
pair.ht <- glht(lme1000, linfct = c("Forest_250 - Agriculture_1000 = 0", "Wetland_1000 - Agriculture_1000 = 0", "Open.uncultivated_1000 - Agriculture_1000 = 0"))

summary(pair.ht) 

confint(pair.ht)

pair.ht <- glht(lme1000, linfct = c(
  "Wetland_1000 - Open.uncultivated_1000 = 0",
  "Forest_250 - Open.uncultivated_1000  = 0"))

summary(pair.ht) 

confint(pair.ht)

#simple regression models
lme1000 <- lmer(log(Biomass+1) ~ 
                  Agriculture_1000 +
                  (1|RouteID) + (1|PilotID), data=allInsects_lowUrban)
summary(lme1000)
lme1000 <- lmer(log(Biomass+1) ~ 
                  Wetland_1000 +
                  (1|RouteID) + (1|PilotID), data=allInsects_lowUrban)
summary(lme1000)
lme1000 <- lmer(log(Biomass+1) ~ 
                  Open.uncultivated_1000 +
                  (1|RouteID) + (1|PilotID), data=allInsects_lowUrban)
summary(lme1000)
lme1000 <- lmer(log(Biomass+1) ~ 
                  Forest_250 +
                  (1|RouteID) + (1|PilotID), data=allInsects_lowUrban)
summary(lme1000)