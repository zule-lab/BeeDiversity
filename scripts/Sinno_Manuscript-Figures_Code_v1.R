###############################################
#### Code for: Sinno et al.                ####
#### Submission to Ecological Applications ####                       
#### Manuscript Figures                    ####
#### December 2023                         #### 
###############################################

library(mgcv) 
library(ggplot2)
library(gridExtra)

# Run Sinno_Manuscript_Code_v1 script to generate data frame
df.fig <- readRDS('output/CleanedData.rds')

# Plot bee richness ~ floral density
Density_plot <- ggplot(df.fig, aes(x = Floral_density, y = Bee_richness)) +
  geom_point(aes(fill=Sampling_round), shape = 21, size = 1) + theme_bw() +
  scale_fill_manual(values=c("black", "white")) + 
  geom_smooth(method="gam", method.args=list(family="nb"), formula=y~x, alpha=0.3, colour="black") +
  theme(aspect.ratio = 1, axis.title.x = element_text(size = 12, face = "bold", colour = "black"),    
        axis.title.y = element_text(size = 12, face ="bold", colour = "black"), legend.position = "none") + 
  labs(x = "Floral density (%)", y = "Wild bee species richness") 
Density_plot
ggsave('graphics/density_plot.png', width = 10, height = 12, units = 'in')

#### Residual Plots ####

# Extract residuals for CWV
M.resid <- gam(Bee_richness ~ Floral_density + Corolla_CWM + Sampling_round +s(Site_label,bs="re"), data = df.fig, family=nb)
resid<-residuals.gam(M.resid, type="deviance") # List of residuals

# Add residuals to data frame
df.fig2<-df.fig[complete.cases(df.fig[,7]) ,] # Remove NAs for floral density, corolla traits
df.fig2 <- cbind(df.fig2,resid) 

# Plot residuals of bee richness and floral density ~ corolla length CWV
CWV_residuals <- ggplot(df.fig2, aes(x = Corolla_CWV, y = resid)) +
  geom_point(aes(fill=Sampling_round), shape = 21, size = 1) + theme_bw() +
  scale_fill_manual(values=c("black", "white")) + 
  scale_y_continuous(limits= c(-3.2,3.2)) +
  geom_smooth(method="lm", alpha=0.3, colour = "black") +
  theme(aspect.ratio = 1, axis.title.x = element_text(size = 12, face = "bold", colour = "black"),    
        axis.title.y = element_text(size = 12, face ="bold", colour = "black"), legend.position = "none") +
  labs(x = "Corolla length CWV", y = "Model residuals") 
CWV_residuals
ggsave('graphics/CWV_residuals.png', width = 10, height = 12, units = 'in')

# Residuals for CWM
M.resid2 <- gam(Bee_richness ~ Floral_density + Corolla_CWV + Sampling_round +s(Site_label,bs="re"), data = df.fig, family=nb)
resid2<-residuals.gam(M.resid2, type="deviance") # List of residuals
df.fig2 <- cbind(df.fig2,resid2) 

CWM_residuals <- ggplot(df.fig2, aes(x = Corolla_CWM, y = resid2)) +
  geom_point(aes(fill=Sampling_round), shape = 21, size = 1) + theme_bw() +
  scale_fill_manual(values=c("black", "white")) + 
  geom_smooth(method="lm", alpha=0.3, colour = "black", linetype="dashed") +
  scale_y_continuous(limits= c(-3.15,3.15)) +
  theme(aspect.ratio = 1, axis.title.x = element_text(size = 12, face = "bold", colour = "black"),    
        axis.title.y = element_text(size = 12, face ="bold", colour = "black"), legend.position = "none") +
  labs(x = "Corolla length CWM", y = "Model residuals") 
CWM_residuals
ggsave('graphics/CWM_residuals.png', width = 10, height = 12, units = 'in')


# Combine plots
p <- grid.arrange(Density_plot, CWV_residuals, CWM_residuals, ncol=3)
ggsave('graphics/all_plots.png', p, width = 14, height = 10, units = 'in')
