
#### load libraries ####
library(ggplot2)
library(dplyr)
library(webr)
library(ggiraphExtra)
library(ggpubr)
library(tidyverse)

library(car) # Companion to Applied Regression
library(lme4) # Fit linear and generalized linear mixed-effects models
library(lmerTest) # Provides p-values in type I, II or III anova and summary tables for lmer model fits (cf. lme4) 
library(nlme) # Fit and compare Gaussian linear and nonlinear mixed-effects models
library(effects) # Graphical and tabular effect displays, e.g., of interactions, for various statistical models with linear predictors
library(MuMIn)#AIC, R2

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified") # if including unspecified/other category

#### load data ####
# load Danish data
#allInsects <- read_rds("data/cleaned_data/data_richness_BINs.RDS")
taxonomy <- readRDS("data/cleaned_data/taxonomy_BIN_combo.rds")
allInsects <- readRDS("data/cleaned_data/data_richness_BINs.rds")

### plot estimated richness for each buffer size ###########

# change the five land covers to be 0-100 instead of 0-1
allInsects_untrans_landcovers <- allInsects # save the untransformed land covers
allInsects[,33:60] <- allInsects[,33:60]*100

#urban
u50 <- ggplot(allInsects)+
  geom_point(aes(x=Urban_50,y=(richness_est)),
             col=landuseCols[1])+
  #scale_y_log10() + 
  theme_bw() +
  xlab("Urban cover at 50m") +ylab("Estimated richness")

u250 <- ggplot(allInsects)+
  geom_point(aes(x=Urban_250,y=(richness_est)),
             col=landuseCols[1])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Urban cover at 250m") +ylab("Estimated richness")

u500 <- ggplot(allInsects)+
  geom_point(aes(x=Urban_500,y=(richness_est)),
             col=landuseCols[1])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Urban cover at 500m") +ylab("Estimated richness")

u1000 <- ggplot(allInsects, aes(x=Urban_1000,y=(richness_est)))+
  geom_point(col=landuseCols[1])+
  #geom_smooth(method = "lm", formula = y ~ x, colour = "darkgrey", fill = "lightgrey") +
  #scale_y_log10() +
  theme_bw() +
  xlab("Urban cover at 1000m") +ylab("Estimated richness")

cowplot::plot_grid(u50,u250,u500,u1000,ncol=1)

#agriculture
a50 <- ggplot(allInsects)+
  geom_point(aes(x=Agriculture_50,y=(richness_est)),
             col=landuseCols[2])+
  #scale_y_log10() + 
  theme_bw() +
  xlab("Farmland cover at 50m") +ylab("Estimated richness")

a250 <- ggplot(allInsects)+
  geom_point(aes(x=Agriculture_250,y=(richness_est)),
             col=landuseCols[2])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Farmland cover at 250m") +ylab("Estimated richness")

a500 <- ggplot(allInsects)+
  geom_point(aes(x=Agriculture_500,y=(richness_est)),
             col=landuseCols[2])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Farmland cover at 500m") +ylab("Estimated richness")

a1000 <- ggplot(allInsects, aes(x=Agriculture_1000,y=(richness_est)))+
  geom_point(col=landuseCols[2])+
  #geom_smooth(method = "lm", formula = y ~ x, colour = "darkgrey", fill = "lightgrey") +
  #scale_y_log10() +
  theme_bw() +
  xlab("Farmland cover at 1000m") +ylab("Estimated richness")

cowplot::plot_grid(a50,a250,a500,a1000,ncol=1)

# Grassland
g50 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Open.uncultivated.land_50,y=(richness_est)),
             col=landuseCols[3])+
  #scale_y_log10() + 
  theme_bw() +
  xlab("Grassland cover at 50m") +ylab("Estimated richness")

g250 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Open.uncultivated.land_250,y=(richness_est)),
             col=landuseCols[3])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Grassland cover at 250m") +ylab("Estimated richness")

g500 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Open.uncultivated.land_500,y=(richness_est)),
             col=landuseCols[3])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Grassland cover at 500m") +ylab("richness_est")

g1000 <- ggplot(allInsects_trans_landcovers, aes(x=Open.uncultivated.land_1000,y=(richness_est)))+
  geom_point(col=landuseCols[3])+
  #geom_smooth(method = "lm", formula = y ~ x, colour = "darkgrey", fill = "lightgrey") +
  #scale_y_log10() +
  theme_bw() +
  xlab("Grassland cover at 1000m") +ylab("Estimated richness")

cowplot::plot_grid(g50,g250,g500,g1000,ncol=1)

# wetland
w50 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Wetland_50,y=(richness_est)),
             col=landuseCols[4])+
  #scale_y_log10() + 
  theme_bw() +
  xlab("Wetland cover at 50m") +ylab("Estimated richness")

w250 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Wetland_250,y=(richness_est)),
             col=landuseCols[4])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Wetland cover at 250m") +ylab("Estimated richness")

w500 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Wetland_500,y=(richness_est)),
             col=landuseCols[4])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Wetland cover at 500m") +ylab("Estimated richness")

w1000 <- ggplot(allInsects_trans_landcovers, aes(x=Wetland_1000,y=(richness_est))) +
  geom_point(col=landuseCols[4])+
  #geom_smooth(method = "lm", formula = y ~ x, colour = "darkgrey", fill = "lightgrey") +
  #scale_y_log10() +
  theme_bw() +
  xlab("Wetland cover at 1000m") +ylab("Estimated richness")

cowplot::plot_grid(w50,w250,w500,w1000,ncol=1)

#forest
f50 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Forest_50,y=(richness_est)),
             col=landuseCols[5])+
  #scale_y_log10() + 
  theme_bw() +
  xlab("Forest cover at 50m") +ylab("Estimated richness")

f250 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Forest_250,y=(richness_est)),
             col=landuseCols[5])+
  #scale_y_log10() +
  theme_bw() +
  xlab("Forest cover at 250m") +ylab("Estimated richness")

f500 <- ggplot(allInsects_trans_landcovers)+
  geom_point(aes(x=Forest_500,y=(richness_est)),
             col=landuseCols[5])+
  ##scale_y_log10() +
  theme_bw() +
  xlab("Forest cover at 500m") +ylab("Estimated richness")

f1000 <- ggplot(allInsects_trans_landcovers, aes(x=Forest_1000,y=(richness_est)))+
  geom_point(col=landuseCols[5])+
  #geom_smooth(method = "lm", formula = y ~ x, colour = "darkgrey", fill = "lightgrey") +
  #scale_y_continuous() +
  theme_bw() +
  xlab("Forest cover at 1000m") +ylab("Estimated richness")

cowplot::plot_grid(f50,f250,f500,f1000,ncol=1)


cowplot::plot_grid(u1000,a1000,g1000,w1000,f1000,ncol=1, labels = c("A", "B", "C", "D", "E"))
ggsave("plots/Landcover_percent_with_regline.png",width=5,height=15)

### Buffer effect plots ######################
str(allInsects)
hist(allInsects$richness_est)
hist(log(allInsects$richness_est))
hist(sqrt(allInsects$richness_est)) # DK: sqrt works for normalisation!

#function to extract summary components of simpel buffer effect models
getEffect <- function(model){
  coefs <- summary(model)$coef[2,]
  temp <- confint(model)
  cis <- temp[5,]
  data.frame(t(coefs),t(cis))
}

##### Estimated richness #####

# we will not normalise the response to ease interpretation, although square root transformation normalises the distribution

#agriculture
hist(allInsects$Agriculture_1000)
hist(sqrt(allInsects$Agriculture_1000)) # does not help

lme50 <- lmer(richness_est ~ (Agriculture_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(richness_est ~ (Agriculture_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(richness_est ~ (Agriculture_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(richness_est ~ (Agriculture_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outAgri <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outAgri <- as.data.frame(outAgri)
outAgri$Buffer <- c(50,250,500,1000)
outAgri$Land_use <- "Farmland"

#urban
hist(allInsects$Urban_1000)#should we sqrt it?
hist(sqrt(allInsects$Urban_1000)) #sqrt is between

lme50 <- lmer(richness_est ~ (Urban_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(richness_est ~ (Urban_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(richness_est ~ (Urban_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(richness_est ~ (Urban_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outUrban <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUrban <- as.data.frame(outUrban)
outUrban$Buffer <- c(50,250,500,1000)
outUrban$Land_use <- "Urban"

#Open.uncultivated
hist(allInsects$Open.uncultivated.land_250)

lme50 <- lmer(richness_est ~ (Open.uncultivated.land_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(richness_est ~ (Open.uncultivated.land_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(richness_est ~ (Open.uncultivated.land_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(richness_est ~ (Open.uncultivated.land_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outOpen.uncultivated <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outOpen.uncultivated <- as.data.frame(outOpen.uncultivated)
outOpen.uncultivated$Buffer <- c(50,250,500,1000)
outOpen.uncultivated$Land_use <- "Grassland"

#forest
hist(allInsects$Forest_250)

lme50 <- lmer(richness_est ~ (Forest_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(richness_est ~ (Forest_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(richness_est ~ (Forest_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(richness_est ~ (Forest_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outForest<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outForest <- as.data.frame(outForest)
outForest$Buffer <- c(50,250,500,1000)
outForest$Land_use <- "Forest"

#wetland
hist(allInsects$Wetland_1000)

lme50 <- lmer(richness_est ~ (Wetland_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay +  
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(richness_est ~ (Wetland_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(richness_est ~ (Wetland_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(richness_est ~ (Wetland_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outWetland<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outWetland <- as.data.frame(outWetland)
outWetland$Buffer <- c(50,250,500,1000)
outWetland$Land_use <- "Wetland"

#unspecified
hist(allInsects$Unspecified.land.cover_1000)
hist(sqrt(allInsects$Unspecified.land.cover_1000))

lme50 <- lmer(richness_est ~ Unspecified.land.cover_50 + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay +  
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(richness_est ~ Unspecified.land.cover_250 + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(richness_est ~ Unspecified.land.cover_500 + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(richness_est ~ Unspecified.land.cover_1000 + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outUnspecified<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUnspecified <- as.data.frame(outUnspecified)
outUnspecified$Buffer <- c(50,250,500,1000)
outUnspecified$Land_use <- "Unspecified"

#combine all effects
outAll <- rbind(outForest,outAgri,outUrban,outWetland,outOpen.uncultivated)
outAll$Land_use <- factor(outAll$Land_use,levels=landuseOrder)

ggplot(outAll)+
  geom_crossbar(aes(x=factor(Buffer),y=Estimate,
                    ymin=X2.5..,ymax=X97.5..,
                    fill=Land_use),
                width=0.5)+
  facet_wrap(~Land_use,scales="free",ncol=1)+
  scale_fill_manual(values=landuseCols)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")+
  geom_hline(yintercept=0,colour="black",linetype="dashed")+
  xlab("Buffer size (m)") + ylab("Effect of land cover on\n estimated richness") + labs(subtitle = "A") + theme(plot.subtitle=element_text(size=18, face="bold", color="black"))

ggsave("plots/landcover_estrichness_buffer_effect_size.png",width=4,height=8)

##### Biomass #####

#agriculture
lme50 <- lmer(totalBiomass_mg ~ (Agriculture_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(totalBiomass_mg ~ (Agriculture_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(totalBiomass_mg ~ (Agriculture_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(totalBiomass_mg ~ (Agriculture_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outAgri <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outAgri <- as.data.frame(outAgri)
outAgri$Buffer <- c(50,250,500,1000)
outAgri$Land_use <- "Farmland"

#urban
# singularity an issue for all buffer zone lower than 1000
lme50 <- lmer(totalBiomass_mg ~ (Urban_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(totalBiomass_mg ~ (Urban_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(totalBiomass_mg ~ (Urban_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(totalBiomass_mg ~ (Urban_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outUrban <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUrban <- as.data.frame(outUrban)
outUrban$Buffer <- c(50,250,500,1000)
outUrban$Land_use <- "Urban"

#Open.uncultivated
lme50 <- lmer(totalBiomass_mg ~ (Open.uncultivated.land_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(totalBiomass_mg ~ (Open.uncultivated.land_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(totalBiomass_mg ~ (Open.uncultivated.land_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(totalBiomass_mg ~ (Open.uncultivated.land_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outOpen.uncultivated <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outOpen.uncultivated <- as.data.frame(outOpen.uncultivated)
outOpen.uncultivated$Buffer <- c(50,250,500,1000)
outOpen.uncultivated$Land_use <- "Grassland"

#forest
lme50 <- lmer(totalBiomass_mg ~ (Forest_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay + 
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(totalBiomass_mg ~ (Forest_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(totalBiomass_mg ~ (Forest_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(totalBiomass_mg ~ (Forest_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outForest<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outForest <- as.data.frame(outForest)
outForest$Buffer <- c(50,250,500,1000)
outForest$Land_use <- "Forest"

#wetland
# model failed to coverge for the 50 m buffer for DK
lme50 <- lmer(totalBiomass_mg ~ (Wetland_50) + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay +  
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(totalBiomass_mg ~ (Wetland_250) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(totalBiomass_mg ~ (Wetland_500) + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(totalBiomass_mg ~ (Wetland_1000) + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outWetland<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outWetland <- as.data.frame(outWetland)
outWetland$Buffer <- c(50,250,500,1000)
outWetland$Land_use <- "Wetland"

#unspecified
# singularity for the 500 m buffer
lme50 <- lmer(totalBiomass_mg ~ Unspecified.land.cover_50 + Time_band + Year + 
                Time_band:cnumberTime + cTL + cyDay +  
                (1|RouteID_JB) + (1|PID), data=allInsects)
lme250 <- lmer(totalBiomass_mg ~ Unspecified.land.cover_250 + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme500 <- lmer(totalBiomass_mg ~ Unspecified.land.cover_500 + Time_band + Year + 
                 Time_band:cnumberTime + cTL + cyDay +  
                 (1|RouteID_JB) + (1|PID), data=allInsects)
lme1000 <- lmer(totalBiomass_mg ~ Unspecified.land.cover_1000 + Time_band + Year + 
                  Time_band:cnumberTime + cTL + cyDay +  
                  (1|RouteID_JB) + (1|PID), data=allInsects)
outUnspecified<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUnspecified <- as.data.frame(outUnspecified)
outUnspecified$Buffer <- c(50,250,500,1000)
outUnspecified$Land_use <- "Unspecified"

#combine all effects
outAll <- rbind(outForest,outAgri,outUrban,outWetland,outOpen.uncultivated)
outAll$Land_use <- factor(outAll$Land_use,levels=landuseOrder)

ggplot(outAll)+
  geom_crossbar(aes(x=factor(Buffer),y=Estimate,
                    ymin=X2.5..,ymax=X97.5..,
                    fill=Land_use),
                width=0.5)+
  facet_wrap(~Land_use,scales="free",ncol=1)+
  scale_fill_manual(values=landuseCols)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")+
  geom_hline(yintercept=0,colour="black",linetype="dashed")+
  xlab("Buffer size (m)") + ylab("Effect of land cover on\n total sample biomass") + labs(subtitle = "B")+ theme(plot.subtitle=element_text(size=18, face="bold", color="black"))

ggsave("plots/landcover_biomass_buffer_effect_size.png",width=4,height=8)

### taxonomy pie-donut chart ####
tax_red <- taxonomy %>% select(order, family)

# get counts and frequencies of combination of order and family, and remove <5% frequencies
tax_red_pd <- tax_red %>%
  group_by(order, family) %>%
  drop_na(family) %>% 
  summarise(n = n()) #%>% mutate(freq = n / sum(n)) 

# look at the object and rename the low abundant orders
namestochange <- tax_red_pd %>%
  group_by(order, family) %>%
  arrange(desc(n))

# hard-coded - not ideal! we only keep in order >5%
pd_other <-
  tax_red_pd %>% dplyr::mutate(
    Order = dplyr::recode(
      order,
      "Lepidoptera" = "Other order",
      "Trombidiformes" = "Other order",
      "Psocodea" = "Other order",
      "Mesostigmata" = "Other order",
      "Araneae" = "Other order",
      "Sarcoptiformes"= "Other order",
      "Trichoptera"= "Other order",
      "Odonata"= "Other order",
      "Ephemeroptera"= "Other order",
      "Neuroptera"= "Other order",
      "Orthoptera"= "Other order",
      "Plecoptera"= "Other order",
      "Thysanoptera"= "Other order",
      "Dermaptera"= "Other order",
      "Mecoptera"= "Other order",
      "Opiliones" = "Other order",
      "Pseudoscorpiones" = "Other order"
    )
  ) 

str(pd_other)
pd_other <- as_tibble(pd_other)

# work around -pick top four of each of the main orders and write below
pd_other <- pd_other %>% select(Order, family, n)

# work on each order level - remove <5% families within each order

# Diptera
pd_other %>% 
  filter(Order == "Diptera") %>%  
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) 

pd_diptera <-
  pd_other %>% 
  filter(Order == "Diptera") %>% 
  mutate(freq = n / sum(n)) %>% 
  filter(!freq < 0.05)

pd_diptera %>% top_n(4, n)
  
pd_diptera <-
  pd_diptera %>% 
  group_by(Order, family) %>% 
  summarise(freq = sum(n))

pd_diptera <- as_tibble(pd_diptera)

# Hymenoptera
pd_other %>% 
  filter(Order == "Hymenoptera") %>%  
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) 

pd_hymenoptera <-
  pd_other %>% 
  filter(Order == "Hymenoptera") %>% 
  mutate(freq = n / sum(n)) %>% 
  filter(!freq < 0.05)

pd_hymenoptera %>% top_n(4, n)

#pd_hymenoptera$final_family = replace(x = pd_hymenoptera$family, list =  !pd_hymenoptera$family %in% c('Braconidae','Ichneumonidae','Megaspilidae','Tenthredinidae'),values =  'Other Hymenoptera')

pd_hymenoptera <- pd_hymenoptera %>% 
  group_by(Order, family) %>% 
  summarise(freq = sum(n))

pd_hymenoptera <- as_tibble(pd_hymenoptera)

# Coleoptera
pd_other %>% 
  filter(Order == "Coleoptera") %>%  
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) 

pd_coleoptera <-
  pd_other %>% 
  filter(Order == "Coleoptera") %>% 
  mutate(freq = n / sum(n)) %>% 
  filter(!freq < 0.05)

pd_coleoptera %>% top_n(4, n)

pd_coleoptera <- pd_coleoptera %>% 
  group_by(Order, family) %>% 
  summarise(freq = sum(n))

pd_coleoptera <- as_tibble(pd_coleoptera)

# Hemiptera
pd_other %>% 
  filter(Order == "Hemiptera") %>%  
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) 

pd_hemiptera <-
  pd_other %>% 
  filter(Order == "Hemiptera") %>% 
  mutate(freq = n / sum(n)) %>% 
  filter(!freq < 0.05)

pd_hemiptera %>% top_n(4, n)

pd_hemiptera <- pd_hemiptera %>% 
  group_by(Order, family) %>% 
  summarise(freq = sum(n))

pd_hemiptera <- as_tibble(pd_hemiptera)

# other orders - only keep two
pd_other %>% 
  filter(Order == "Other order") %>%  
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) 

dp_other <-
  pd_other %>% 
  filter(Order == "Other order") %>% 
  mutate(freq = n / sum(n)) %>% 
  filter(!freq < 0.05)

dp_other %>% top_n(4, n)

dp_other <- dp_other %>% 
  group_by(Order, family) %>% 
  summarise(freq = sum(n))

dp_other <- as_tibble(dp_other)

dp <- rbind(pd_diptera, pd_hymenoptera)
dp <- rbind(dp, pd_coleoptera)
dp <- rbind(dp, pd_hemiptera)
dp <- rbind(dp, dp_other)

str(dp)

# Pie-Donut chart with the webr package
PieDonut(dp, aes(Order, family, count = freq), r0 = 0.25, maxx = 1.3, r1 = 0.8, r2 = 1, labelposition = 1, pieLabelSize = 4, donutLabelSize = 3, showRatioThreshold = 0.01, ratioByGroup = F) 

tiff(file="plots/piedonut_taxonomy.tiff",
     width=18, height=12, units="in", res=400)
PieDonut(dp, aes(Order, family, count = freq), r0 = 0.25, maxx = 1.3, r1 = 0.8, r2 = 1, labelposition = 1, pieLabelSize = 4, donutLabelSize = 3, showRatioThreshold = 0.01) 
dev.off()

# using the ggiraphExtra package
ggPieDonut(data = dp, mapping = aes(pies = Order, donuts = family, count = freq), interactive = F, labelposition = 1) 