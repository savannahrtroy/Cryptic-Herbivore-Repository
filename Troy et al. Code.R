#########################################################
#### Savannah Troy et al. Cryptic Herbivore Analysis ####
#########################################################

#### Library ####

# Packages

library(tidyverse)
library(emmeans)
library(blme)
library(glmmTMB)
library(car)

# Files

seed <- readRDS("ProcB_CodePackage/bothseason_seed.rds")
germ <- readRDS("ProcB_CodePackage/bothseason_germinant.rds")

#seed$obs <- 1:nrow(seed)
#germ$obs <- 1:nrow(germ)

seed$TrialNo <- as.factor(seed$TrialNo)
germ$TrialNo <- as.factor(germ$TrialNo)

#### Main Analysis ####

# Response variable is proportion of unconsumed seeds out of total deployed

#Use glmmTMB template model builder to determine dimensions
both_failed_seed <- glmmTMB(FailSeed ~ PlotType + Treatment + season + offset(log(exposure)) +
                              
                              PlotType:Treatment + PlotType:season + season:Treatment +
                              
                              (1|PlotID/Cluster) + (1|TrialNo),
                            
                            family = poisson, data = seed)

dim(model.matrix(both_failed_seed)) #number of columns specifies the number of parameters- here, it's 18

# model response is ratio of remaining seeds to consumed seeds

m1both_seed <- bglmer(cbind(D5SW,FailSeed) ~ PlotType + Treatment + season + offset(log(exposure)) +
                        PlotType:Treatment + PlotType:season + season:Treatment + (1|PlotID/Cluster) + (1|TrialNo), 
                      family = binomial,
                      fixef.prior = normal(cov=diag(9,18)), 
                      control = glmerControl(optimizer="bobyqa"), 
                      data = seed)
m1both_germ <- bglmer(cbind(D5GW,FailGerm) ~ PlotType + Treatment + season + offset(log(exposure)) +
                        PlotType:Treatment + PlotType:season + season:Treatment + (1|PlotID/Cluster) + (1|TrialNo), 
                      family = binomial,
                      fixef.prior = normal(cov=diag(9,18)), 
                      control = glmerControl(optimizer="bobyqa"), 
                      data = germ)

# goodness of fit

library(DHARMa)
sim_res <- simulateResiduals(fittedModel = m1both_germ, n=1000)
sim_res_seed <- simulateResiduals(fittedModel = m1both_seed, n=1000)
plot(sim_res)
plot(sim_res_seed)

#emmeans
emm_int_sb <- emmeans(m1both_seed, ~ Treatment | PlotType | season, offset = log(1))
emm_int_sb
emm_int_gb <- emmeans(m1both_germ, ~ Treatment | PlotType | season, offset = log(1))
emm_int_gb

# General effect: C, B, F, RB
# Contrast for each taxon plus a contrast for pooled effect of all cryptic taxa
eff_bs <- contrast(emm_int_sb, list(invertebrate = c(0, 0, 1, -1),
                                    bird = c(-1, 1, 0, 0),
                                    rodent = c(0, -1, 0, 1),
                                    all=c(-1, 0, 1, 0)))

output_sb <- confint(eff_bs,level=0.95) 

eff_bg <- contrast(emm_int_gb, list(invertebrate = c(0, 0, 1, -1),
                                    bird = c(-1, 1, 0, 0),
                                    rodent = c(0, -1, 0, 1),
                                    all=c(-1, 0, 1, 0)))
output_gb <- confint(eff_bg,level=0.95)

# Large mammal contrast

lgtest <- emmeans(m1both_seed, ~ Treatment * PlotType | season, offset = log(1))
lgtest.eff <- contrast(lgtest, list(large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
lgmam.seed <- confint(lgtest.eff)

lg.germ <- emmeans(m1both_germ, ~ Treatment * PlotType | season, offset = log(1))
lg.germ.eff <- contrast(lg.germ, list(large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
lg.germ.output <- confint(lg.germ.eff)

#add these in to the output objects so it will be easier to make figures
lgmam.seed$PlotType <- NA
output_sb <- rbind(output_sb, lgmam.seed)
lg.germ.output$PlotType <- NA
output_gb <- rbind(output_gb, lg.germ.output)

Anova(m1both_germ)
Anova(m1both_seed)

#### Analyses in Supplement ####

### Wet season only analysis: day 5 vs day 10

wseed <- readRDS("ProcB_CodePackage/wetseason_seed.rds")
wgerm <- readRDS("ProcB_CodePackage/wetseason_germinant.rds")

# Day 5 models
# all the trials were same length during wet season, so no offset of exposure term in these formulas

ws5 <- bglmer(cbind(D5SW,failure5) ~ PlotType*Treatment  + (1|PlotID) + (1|TrialNo), 
              family = binomial,
              fixef.prior = normal(cov = diag(9,12)),  
              control = glmerControl(optimizer="bobyqa"), 
              data = wseed)
wg5 <- bglmer(cbind(D5GW,failure5) ~ PlotType*Treatment  + (1|PlotID) + (1|TrialNo), 
              family = binomial,
              fixef.prior = normal(cov = diag(9,12)), 
              control = glmerControl(optimizer="bobyqa"), 
              data = wgerm)
#emmeans
emm_int_sw5 <- emmeans(ws5, ~ Treatment | PlotType, offset = log(1))
emm_int_sw5
emm_int_gw5 <- emmeans(wg5, ~ Treatment | PlotType, offset = log(1))
emm_int_gw5
#General effect: C, B, F, RB
eff_ws5 <- contrast(emm_int_sw5, list(invertebrate = c(0, 0, 1, -1),
                                      bird = c(-1, 1, 0, 0),
                                      rodent = c(0, -1, 0, 1),
                                      all=c(-1, 0, 1, 0)))
output_sw5 <- confint(eff_ws5,level=0.95) 
output_sw5
eff_wg5 <- contrast(emm_int_gw5, list(invertebrate = c(0, 0, 1, -1),
                                      bird = c(-1, 1, 0, 0),
                                      rodent = c(0, -1, 0, 1),
                                      all=c(-1, 0, 1, 0)))
output_gw5 <- confint(eff_wg5,level=0.95)
output_gw5

#Add large mammals
w5s.lg <- emmeans(ws5, ~ Treatment * PlotType, offset = log(1))
w5slg.eff <- contrast(w5s.lg, list(large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w5slg.output <- confint(w5slg.eff,level=0.95)

w5g.lg <- emmeans(wg5, ~ Treatment * PlotType, offset = log(1))
w5glg.eff <- contrast(w5g.lg, list(large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w5glg.output <- confint(w5glg.eff,level=0.95)

#add these in to the output objects so it will be easier to make figures
w5slg.output$PlotType <- "Large Mammals"
output_sw5 <- rbind(output_sw5, w5slg.output)
w5glg.output$PlotType <- "Large Mammals"
output_gw5 <- rbind(output_gw5, w5glg.output)

# Day 10 models

mod_wsbayes <- bglmer(cbind(D10SW,failure10) ~ PlotType*Treatment  + (1|PlotID) + (1|TrialNo), 
                      family = binomial,
                      fixef.prior = normal(cov = diag(9,12)), 
                      control = glmerControl(optimizer="bobyqa"), 
                      data = wseed)
mod_wgbayes <- bglmer(cbind(D10GW,failure10) ~ PlotType*Treatment  + (1|PlotID) + (1|TrialNo), 
                      family = binomial,
                      fixef.prior = normal(cov = diag(9,12)), 
                      control = glmerControl(optimizer="bobyqa"), 
                      data = wgerm)
#emmeans
emm_int_sw10 <- emmeans(mod_wsbayes, ~ Treatment | PlotType, offset = log(1))
emm_int_sw10
emm_int_gw10 <- emmeans(mod_wgbayes, ~ Treatment | PlotType, offset = log(1))
emm_int_gw10
#General effect: C, B, F, RB
eff_ws10 <- contrast(emm_int_sw10, list(invertebrate = c(0, 0, 1, -1),
                                        bird = c(-1, 1, 0, 0),
                                        rodent = c(0, -1, 0, 1),
                                        all=c(-1, 0, 1, 0)))

output_sw10 <- confint(eff_ws10,level=0.95) 
output_sw10
eff_wg10 <- contrast(emm_int_gw10, list(invertebrate = c(0, 0, 1, -1),
                                        bird = c(-1, 1, 0, 0),
                                        rodent = c(0, -1, 0, 1),
                                        all=c(-1, 0, 1, 0)))
output_gw10 <- confint(eff_wg10,level=0.95)
output_gw10

#Add large mammals
w10s.lg <- emmeans(mod_wsbayes, ~ Treatment * PlotType, offset = log(1))
w10slg.eff <- contrast(w10s.lg, list(large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w10slg.output <- confint(w10slg.eff,,level=0.95)

w10g.lg <- emmeans(mod_wgbayes, ~ Treatment * PlotType, offset = log(1))
w10glg.eff <- contrast(w10g.lg, list(large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w10glg.output <- confint(w10glg.eff,,level=0.95)

#add these in to the output objects so it will be easier to make figures
w10slg.output$PlotType <- "Large Mammals"
output_sw10 <- rbind(output_sw10, w10slg.output)
w10glg.output$PlotType <- "Large Mammals"
output_gw10 <- rbind(output_gw10, w10glg.output)

# Analysis of deviance results for SI table 1
Anova(mod_wsbayes)
Anova(mod_wgbayes)


### Rodent Minimum Number Alive

rodent.all <- read.csv("ProcB_CodePackage/rodentsurvey.csv")
rodent.mna <- rodent.all %>% group_by(Plot,season,treatment) %>% summarise(MNA = n_distinct(UniqueID))
rodent.treatgroup <- rodent.mna %>% group_by(season,treatment) %>% summarise(mean.mna = mean(MNA),
                                                                             se.mna = sd(MNA,na.rm=T)/sqrt(3))
# Poisson mna ~ treatment

rodent.mna$season <- as.factor(rodent.mna$season)
rodent.mna$treatment <- as.factor(rodent.mna$treatment)

rod.mod <- glm(MNA ~ treatment + season + treatment*season,
               data=rodent.mna,
               family="poisson")

#anoda
rod.test <- anova(rod.mod)
rod.test

#correlation with consumption effects
rod.eff.g <- output_gb %>% filter(contrast == "rodent")
rod.eff.g$response <- "Germinant"
rod.eff.s <- output_sb %>% filter(contrast == "rodent")
rod.eff.s$response <- "Seed"
rod.eff <- rbind(rod.eff.g, rod.eff.s)

names(rodent.treatgroup)[2] <- "PlotType"
# match case
rodent.treatgroup$PlotType <- case_when(rodent.treatgroup$PlotType=="Control" ~ "control",
                                        rodent.treatgroup$PlotType=="Fenced" ~ "fence",
                                        rodent.treatgroup$PlotType=="Habitat-Modified" ~ "habitat")

rod.eff <- merge(rod.eff,rodent.treatgroup,by=c("season","PlotType"),all.x=T)
rod.eff$logmna <- log(rod.eff$mean.mna+1)
cor(rod.eff$estimate[rod.eff$response=="Germinant"],rod.eff$logmna[rod.eff$response=="Germinant"])
cor(rod.eff$estimate[rod.eff$response=="Seed"],rod.eff$logmna[rod.eff$response=="Seed"])



### Vegetation Structure

canopy <- read.csv("ProcB_CodePackage/canopycover.csv")
shrub <- read.csv("ProcB_CodePackage/shrubcover.csv")

# get a key for plot number - treatment type

plot.treat <- germ %>% select(PlotID,PlotType) %>% distinct()

sh.sel <- shrub %>% dplyr::select(plot, shrub) %>% group_by(plot) %>% summarise(avg=mean(shrub),sd=sd(shrub)) %>%
  filter(plot %in% c(10,11,18,19,8,6,5,3,2))
ca.sel <- canopy %>% dplyr::select(plot,cover) %>% group_by(plot) %>% summarise(avg=mean(cover),sd=sd(cover)) %>%
  filter(plot %in% c(10,11,18,19,8,6,5,3,2))
ca.sel$length <- 22
sh.sel$length <- 10
ca.sel$SE <- ca.sel$sd / sqrt(ca.sel$length)
sh.sel$SE <- sh.sel$sd / sqrt(sh.sel$length)
sh.sel$ID <- "Shrub Density"
ca.sel$ID <- "Canopy Closure"

allveg <- rbind(ca.sel,sh.sel)
allveg <- merge(allveg, plot.treat, by.x="plot", by.y="PlotID")

#ANOVA
canopy.sub <- canopy %>% dplyr::select(plot,cover) %>% filter(plot %in% c(10,11,18,19,8,6,5,3,2)) %>% merge(plot.treat,by.x="plot",by.y="PlotID") 
ca.aov <- aov(cover ~ PlotType, data=canopy.sub)
summary(ca.aov)
ca.tu <- TukeyHSD(ca.aov)
plot(ca.tu)

shrub.sub <- shrub %>% dplyr::select(plot,shrub) %>% filter(plot %in% c(10,11,18,19,8,6,5,3,2)) %>% merge(plot.treat,by.x="plot",by.y="PlotID") 
sh.aov <- aov(shrub ~ PlotType, data=shrub.sub)
summary(sh.aov)
sh.tu <- TukeyHSD(sh.aov)
plot(sh.tu)


