library(tidyverse)
library(ggcorrplot)
library(GGally)

resFolders <- "data/cleaned_data"

#read in biodiversity data
df <- readRDS(paste(resFolders,"data_richness.RDS", sep="/"))


#subset to response metrics
df <- df %>%
        dplyr::select(contains("richness"),
                      contains("shannon"),
                      contains("Biomass"),
                      contains("n_reads")) %>%
        select(-est_richness_model, -biomassUncertainty)

#sep for Year and Time_band eventually
ggpairs(df)

#get and plot correlation matrix
corr <- round(cor(df, use="pairwise.complete.obs"), 2)
ggcorrplot(corr, hc.order = TRUE, type = "lower",lab = TRUE,
           outline.col = "white",
           ggtheme = ggplot2::theme_minimal,
           colors = c("#6D9EC1", "white", "#E46726"))
