###############################################
#### Code for: Sinno et al.                ####
#### Submission to Ecological Applications ####                       
#### Trait metrics and mixed models        ####
#### December 2023                         #### 
###############################################

# Load packages
library(mgcv) 
library(dplyr)
library(tidyr)
library(purrr)
library(GGally)
#library(ggplot2)

# Load data frames
df.bees <- read.csv('input/Sinno_Bees_Summarized_2020.csv') # Read in "Sinno_Bees_Summarized_2020"
df.floral <- read.csv('input/Sinno_FloralTraits_2020.csv') # Read in "Sinno_FloralTraits_2020"
head(df.bees)
head(df.floral)

# Calculate Metrics for Floral Traits ####

# create function for FDvar (for corolla and nectar)
FDvar.C <- function(df.floral){
  wi <- df.floral$Abundance/sum(df.floral$Abundance)
  lnxi <- log(df.floral$Mean_Corolla)
  lnx.bar <- sum(wi * lnxi)
  V <- sum(wi*(lnxi - lnx.bar)^2, na.rm = TRUE)
  FDvar.C <- (2/pi)*atan(5*V)
  return(FDvar.C)
}

FDvar.N <- function(df.floral){
  wi <- df.floral$Abundance/sum(df.floral$Abundance)
  lnxi <- log(df.floral$Mean_Nectar)
  lnx.bar <- sum(wi * lnxi)
  V <- sum(wi*(lnxi - lnx.bar)^2, na.rm = TRUE)
  FDvar.N <- (2/pi)*atan(5*V)
  return(FDvar.N)
}

PlotMetrics<- df.floral %>%
  group_by(Plot_ID) %>%
  dplyr::summarise(
    Floral_density = sum(Abundance),
    Floral_richness = n(),
    Corolla_CWM = weighted.mean(Mean_Corolla, Abundance, na.rm=T),
    Nectar_CWM = weighted.mean(Mean_Nectar, Abundance, na.rm=T)
  )

Corolla_CWV <- df.floral %>% 
  group_by(Plot_ID) %>%
  do(data.frame(Corolla_CWV=FDvar.C(.)))

df.floral2<-df.floral %>% drop_na(Mean_Nectar) #Remove rows with NA first
Nectar_CWV <- df.floral2 %>% 
  group_by(Plot_ID) %>%
  do(data.frame(Nectar_CWV=FDvar.N(.)))

# merge dataframes to create one data frame for analysis and plotting
df <- list(df.bees, PlotMetrics, Corolla_CWV, Nectar_CWV) # create list of data frames
df <- df %>% reduce(full_join, by='Plot_ID') #merge all data frames in list
head(df)

# Make sampling round and site factors
df$Sampling_round <- factor(df$Sampling_round) # make factor
class(df$Sampling_round) # confirm
df$Site_label <- factor(df$Site_label) # make factor
class(df$Site_label) # confirm

# save
saveRDS(df, 'output/CleanedData.rds')

# Running generalized linear mixed models (GLMMs), with a negative binomial distribution (for count data)
# Using package mgcv 

# Full model (including ALL predictor variables)

MFull.Corolla <- gam(Bee_richness ~ Floral_richness + Floral_density + Corolla_CWV + Corolla_CWM + 
                       Sampling_round  + s(Site_label,bs="re"), data = df, family=nb)
summary(MFull.Corolla)

# Look at variable correlations

ggpairs(df[,c(6:12)])

# Corolla CWV is strongly correlated with floral richness (>0.6), as we would expect (both measures of "diversity"). Note: Floral richness & density also somewhat correlated (>0.5)
# Best not to include both CWV and floral richness in the same model
# To decide which of the diversity metrics to include in the model, drop each in turn from the Full model (above)

# Drop CWV (leaving Floral Richness)
M2.Corolla <- gam(Bee_richness ~ Floral_density + Floral_richness + Corolla_CWM + 
                    Sampling_round  + s(Site_label,bs="re"), data = df, family=nb)
summary(M2.Corolla) # Floral richness "marginally" significant
gam.check(M2.Corolla) # Residuals look good

# Drop floral richness (this model tests Trait Diversity + Optimal Trait hypotheses in one model)
M3.Corolla <- gam(Bee_richness ~ Floral_density + Corolla_CWV + Corolla_CWM + 
                    Sampling_round  + s(Site_label,bs="re"), data = df, family=nb)
summary(M3.Corolla) # Model explains similar amount of variation, both CWV and CWM are marginally significant predictors of Bee richness here
gam.check(M3.Corolla) # Residuals look good


#### Nectar traits ####

# Full model (including all variables)
MFull.Nectar <- gam(Bee_richness ~ Floral_density + Floral_richness + Nectar_CWV + Nectar_CWM + 
                      Sampling_round  + s(Site_label,bs="re"), data = df, family=nb)
summary(MFull.Nectar) # only floral density significant, neither floral richness nor nectar traits are significant here
gam.check(MFull.Nectar) # Residuals look good


#### Nectar + Corolla interaction ####

Nectar.Corolla <- gam(Bee_richness ~ Floral_density + Floral_richness + Corolla_CWM*Nectar_CWM + 
                        Sampling_round + s(Site_label,bs="re"), data = df, family=nb)
summary(Nectar.Corolla)
gam.check(Nectar.Corolla)
