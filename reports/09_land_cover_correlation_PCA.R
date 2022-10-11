### Load required libraries ###########################

#library(ggbiplot)
library(corrplot)
library(Hmisc)
library(ggfortify) # to plot PCA
library(MuMIn)

### Load data ###################################################
# NB! choose which country to run the analysis on and load the corresponding data

# load Danish data
allInsects <- read_rds("data/cleaned_data/data_richness_BINs.RDS")

### Correlation plot ##################

# prepare plotting
#par(mfrow = c(2, 2))

#names(allInsects)

# select variables for PCA - we will only use land cover at 1000 m and cStops
#biomass.pca <- prcomp(allInsects[,c(12,142,44:49)], center = TRUE,scale. = TRUE) #choose biomass, stops and land covers
#summary(biomass.pca)
#str(biomass.pca)

#ggbiplot(biomass.pca)

# correlation plot for 1000 m buffer
someInsects <- allInsects[,c(85,54:55,57,59:60)]
colnames(someInsects)
colnames(someInsects) <- c("Stops", "Farmland", "Forest", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

png(file = "plots/correlation_buffers.png", height = 1400, width = 1400, res = 200)
#tiff("plots/correlation_buffers.tiff", units="in", width=7, height=7, res=300)

# prepare plotting
par(mfrow = c(2, 2))

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color",
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 1000 m", mar=c(0,0,1,0))

# correlation plot for 500 m buffer
someInsects <- allInsects[,c(85,47:48,50,52:53)]
colnames(someInsects)
colnames(someInsects) <- c("Stops", "Farmland", "Forest", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color",
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 500 m", mar=c(0,0,1,0))

# correlation plot for 250 m buffer
someInsects <- allInsects[,c(85,40:41,43,45:46)]
colnames(someInsects)
colnames(someInsects) <- c("Stops", "Farmland", "Forest", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color",
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 250 m", mar=c(0,0,1,0))

# correlation plot for 50 m buffer
someInsects <- allInsects[,c(85,33:34,36,38:39)]
colnames(someInsects)
colnames(someInsects) <- c("Stops", "Farmland", "Forest", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color",
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 50 m", mar=c(0,0,1,0))

dev.off()

### PCA #########

# subset data for each land cover buffer prior to analysis
names(allInsects)
cor50 <- allInsects[,c(1,33:34,36,38:39)]
cor250 <- allInsects[,c(1,40:41,43,45:46)]
cor500 <- allInsects[,c(1,47:48,50,52:53)]
cor1000 <- allInsects[,c(1,54:55,57,59:60)]

# 50 m buffer correlation of land cover variables
fit <- princomp(cor50, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# 250 m buffer correlation of land cover variables
fit <- princomp(cor250, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# 250 m buffer correlation of land cover variables
fit <- princomp(cor500, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# 250 m buffer correlation of land cover variables
fit <- princomp(cor1000, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

mydata <- allInsects[,c("cTL",names(allInsects)[grepl("_1000",names(allInsects))])]
names(mydata)
mydata <- mydata[,c(2:3,5,7:8)]
names(mydata)[names(mydata)=="Open.uncultivated.land_1000"] <- "Grassland_1000"
names(mydata)[names(mydata)=="Agriculture_1000"] <- "Farmland_1000"
names(mydata) <- gsub("_1000","",names(mydata))

allInsects$Land_use <- as.character(allInsects$Land_use)
allInsects$Land_use[allInsects$Land_use=="Dryland"] <- "Grassland"
allInsects$Land_use[allInsects$Land_use=="Open uncultivated land"] <- "Grassland"
landuseOrder
allInsects$Land_use <- factor(allInsects$Land_use, levels=landuseOrder)

fit <- princomp(mydata, cor=TRUE)

#pca with rotation
library(psych)
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
pca_rotated <- psych::principal(mydata, rotate="varimax", nfactors=2, scores=TRUE)
print(pca_rotated)

png(file = "plots/rotated_pca_landcovers.png", height = 1800, width = 1800, res = 200)
biplot(pca_rotated, main = "")
dev.off()

#add PCA axes scores to the dataset
allInsects$Urbanization_gradient <- pca_rotated$scores[,1]
allInsects$Forest_gradient <- pca_rotated$scores[,2] # really a forest/farmland gradient

#with ggplot
autoplot(fit)
dk_autoplot <- autoplot(fit, data = allInsects, 
                        loadings = TRUE, 
                        loadings.colour = 'black',
                        loadings.label = TRUE, 
                        loadings.label.size = 5) + 
  #scale_colour_manual(values = landuseCols[1:6])+
  theme_bw()

ggsave("plots/pca_1000_DK.png")
dk_autoplot <- plot_grid(dk_autoplot, labels = "AUTO")
save_plot("plots/pca_1000_DK_numbered.png", dk_autoplot, base_height = 8, base_width = 12)

### model w. gradients ####
#### biomass ####
lme1000 <- lmer(log(totalBiomass_mg+1) ~ 
                  Urbanization_gradient + 
                  Forest_gradient +
                  Time_band + 
                  Time_band:cnumberTime + cTL + cyDay + Year +
                  (1|RouteID_JB) + (1|PID), data=allInsects)
summary(lme1000)
r.squaredGLMM(lme1000)

# first save table to html file
tab_model(lme1000, collapse.ci = F, show.intercept = F, pred.labels = c("Urbanization gradient", "Farmland/Forest gradient", "Time band: evening vs midday", "Potential stops", "Day of year", "Year","Time within midday", "Time within evening"), digits = 2, file = "plots/biomass_PCA_gradients.docx")

#### estimated richness ####
lme1000 <- lmer(richness_est ~ 
                  Urbanization_gradient + 
                  Forest_gradient +
                  Time_band + 
                  Time_band:cnumberTime + cTL + cyDay + Year +
                  (1|RouteID_JB) + (1|PID), data=allInsects)
summary(lme1000)
r.squaredGLMM(lme1000)

tab_model(lme1000, collapse.ci = F, show.intercept = F, pred.labels = c("Urbanization gradient", "Farmland/Forest gradient", "Time band: evening vs midday", "Potential stops", "Day of year", "Year","Time within midday", "Time within evening"), digits = 2, file = "plots/estrichness_PCA_gradients.html")

###pca analysis GOT TO HERE###############################################
#taken from the Quick R website

mydata <- allInsects[,c("cStops",names(allInsects)[grepl("_1000",names(allInsects))])]
names(mydata)
mydata <- mydata[,2:6]
names(mydata)[names(mydata)=="Open.uncultivated_1000"] <- "Grassland_1000"
names(mydata)[names(mydata)=="Agriculture_1000"] <- "Farmland_1000"
names(mydata) <- gsub("_1000","",names(mydata))
allInsects$Land_use <- as.character(allInsects$Land_use)
allInsects$Land_use[allInsects$Land_use=="Dryland"] <- "Grassland"
allInsects$Land_use[allInsects$Land_use=="Open uncultivated"] <- "Grassland"
landuseOrder
allInsects$Land_use <- factor(allInsects$Land_use, levels=landuseOrder)

fit <- princomp(mydata, cor=TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="lines") # scree plot
fit$scores # the principal components
biplot(fit)

#Maximum Likelihood Factor Analysis - from Quick R
# with varimax rotation
fit <- factanal(mydata, 2, rotation="varimax")
print(fit, digits=2, sort=TRUE)
load <- fit$loadings[,1:2]
plot(load,type="n")
text(load,labels=names(mydata),cex=.7)  

#pca with rotation
library(psych)
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
pca_rotated <- psych::principal(mydata, rotate="varimax", nfactors=2, scores=TRUE)
biplot(pca_rotated,main="B: Germany")
print(pca_rotated)
#                       RC1  RC2
#SS loadings           1.49 1.46
#Proportion Var        0.30 0.29


ggsave("plots/pca_with_rotation_1000_DE.png",width=8,height=12)

#add PCA axes scores to the dataset
allInsects$Urbanization_gradient <- pca_rotated$scores[,1]
allInsects$Forest_gradient <- pca_rotated$scores[,2]

#with ggplot
autoplot(fit)
autoplot(fit, data = allInsects, colour = 'Land_use',
         loadings = TRUE, 
         loadings.colour = 'black',
         loadings.label = TRUE, 
         loadings.label.size = 2.5) + 
  scale_colour_manual(values = landuseCols[1:5])+
  theme_bw()

ggsave("plots/pca_1000_DE.png")
Footer
© 2022 GitHub, Inc.
Footer navigation

Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About

InsectMobile_Biomass/07_correlation_PCA.R at master · CecSve/InsectMobile_Biomass