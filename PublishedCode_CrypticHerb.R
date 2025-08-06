# CRYPTIC HERBIVORES DATA ANALYSIS
# Savannah Troy - Both Seasons

# LOADING PACKAGES
library(tidyverse)
library(emmeans)
library(ggplot2)
library(glmmTMB)
library(lubridate)
library(car)
library(blme)
library(cowplot)
library(gridExtra)
library(knitr)
library(data.table)
library(ggh4x)
library(lme4)
library(ggpattern)

#Load data
herb.wet <- read.csv("data/wetseason.csv",stringsAsFactors=FALSE)
herb.dry <- read.csv("data/dryseason.csv", stringsAsFactors=FALSE)

#### Dry only cleaning and analysis ####
#This section recreates the analyses done in the Cryptic Analysis_rf file but with a formatting
#that more easily matches wet season dataset

herb.dry <- herb.dry %>% filter(trial != "Pilot") %>% rename(PlotType=exclosure)
herb.dry$PlotType[herb.dry$PlotType=="fenced"] <- "fence"
herb.dry$PlotType[herb.dry$PlotType=="Habitat-Disturbed"] <- "habitat"

hd.seed <- herb.dry %>% dplyr::select(plot,PlotType,cluster,treatment,trial,trialday,
                                      initial.seeds,seeds.intact,seeds.partial,exposure)
hd.seed <- pivot_wider(hd.seed, names_from=trialday, values_from = c(seeds.intact,seeds.partial))
hd.germ <- herb.dry %>% dplyr::select(plot,PlotType,cluster,treatment,trial,trialday,
                                      initial.germ,no.germ.intact,germ.partial,exposure)
hd.germ <- pivot_wider(hd.germ, names_from=trialday, values_from = c(no.germ.intact,germ.partial))

table(hd.seed$trial,hd.seed$plot)
table(hd.germ$trial, hd.germ$plot)

#Make failure columns
hd.seed$failure <- hd.seed$initial.seeds - hd.seed$seeds.intact_5
hd.germ$failuregerm <- hd.germ$initial.germ - hd.germ$no.germ.intact_5

#Some failure values are negative because more seeds sprouted during the trial in some trays.
#Change negative values to zeros for now
hd.germ$failuregerm[hd.germ$failuregerm<0]<-0 

#Relevel reference categories
hd.germ$treatment <- as.factor(hd.germ$treatment)
hd.germ$PlotType <- as.factor(hd.germ$PlotType)
hd.germ$treatment <- relevel(hd.germ$treatment, ref = "Control")
hd.germ$PlotType <- relevel(hd.germ$PlotType, ref = "control")

hd.seed$treatment <- as.factor(hd.seed$treatment)
hd.seed$PlotType <- as.factor(hd.seed$PlotType)
hd.seed$treatment <- relevel(hd.seed$treatment, ref = "Control")
hd.seed$PlotType <- relevel(hd.seed$PlotType, ref = "control")

#### Wet only cleaning and analysis ####
#Gives comparison between 5 and 10 day duration
#All new data, with a longer trial time (up to 10 days)
#Prep data types

#Add ratio columns
herb.wet$seed.ratio5 <- herb.wet$D5SW/herb.wet$InitSeed
herb.wet$germ.ratio5 <- herb.wet$D5GW/herb.wet$InitGerm

#Rename (doesn't work later, weird factor/character issue)
herb.wet$Treatment[herb.wet$Treatment=="C"] <- "Control"
herb.wet$Treatment[herb.wet$Treatment=="F"] <- "Full"
herb.wet$Treatment[herb.wet$Treatment=="RB"] <- "RodentBird"
herb.wet$Treatment[herb.wet$Treatment=="B"] <- "Bird"

#Relevel reference categories
herb.wet$PlotType <- as.factor(herb.wet$PlotType)
herb.wet$Cluster <- as.factor(herb.wet$Cluster)
herb.wet$Treatment <- as.factor(herb.wet$Treatment)
herb.wet$Treatment <- relevel(herb.wet$Treatment, ref = "Control")
herb.wet$PlotType <- relevel(herb.wet$PlotType, ref = "control")

#Barplots
ggplot(herb.wet, aes(x = Treatment, y = seed.ratio5)) + geom_bar(stat = "identity") + xlab("Treatment") + ylab("Seeds Intact")+ theme_bw()
ggplot(herb.wet, aes(x = PlotType, y = seed.ratio5)) + geom_bar(stat = "identity") + xlab("Exclosure") + ylab("Seeds Intact")+ theme_bw()

ggplot(herb.wet, aes(x = Treatment, y = germ.ratio5)) + geom_bar(stat = "identity") + xlab("Treatment") + ylab("Germinants Intact")+ theme_bw()
ggplot(herb.wet, aes(x = PlotType, y = germ.ratio5)) + geom_bar(stat = "identity") + xlab("Exclosure") + ylab("Germinants Intact")+ theme_bw()

#Germinants, d10 and 5 as model options
wgerm <- herb.wet %>% dplyr::select(PlotID,InitDate,PlotType,Cluster,Treatment,TrialNo,
                                    InitGerm,D1GW,D2GW,D3GW,D4GW,D5GW,D8GW,D10GW)
wgerm <- wgerm %>% filter(TrialNo != 1)
wgerm$failure5 <- wgerm$InitGerm-wgerm$D5GW
wgerm$failure10 <- wgerm$InitGerm-wgerm$D10GW

#Some failure values are negative because more seeds sprouted during the trial in some trays.
#Change negative values to zeros for now
wgerm$failure5[wgerm$failure5<0]<-0 
wgerm$failure10[wgerm$failure10<0] <- 0

#seeds day 5 and 10 as model options
wseed <- herb.wet %>% dplyr::select(PlotID,InitDate,PlotType,Cluster,Treatment,TrialNo,
                                    InitSeed,D1SW,D2SW,D3SW,D4SW,D5SW,D8SW,D10SW)
wseed <- wseed %>% filter(TrialNo != 4)
wseed$failure5 <- wseed$InitSeed-wseed$D5SW
wseed$failure5[wseed$failure5 < 0] <- 0
wseed$failure10 <- wseed$InitSeed-wseed$D10SW
wseed$failure10[wseed$failure10 < 0] <- 0

#5 day model
#Need to create an exposure column (5 for all)
wseed$exposure5 <- 5
wgerm$exposure <- 5

ws5 <- bglmer(cbind(D5SW,failure5) ~ PlotType*Treatment  + (1|PlotID), 
              family = binomial,
              fixef.prior = normal(cov = diag(9,12)), #var = 9, sd = 3
              control = glmerControl(optimizer="bobyqa"), 
              data = wseed)
wg5 <- bglmer(cbind(D5GW,failure5) ~ PlotType*Treatment  + (1|PlotID), 
              family = binomial,
              fixef.prior = normal(cov = diag(9,12)), #var = 9, sd = 3
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
                                      rodent = c(0, -1, 0, 1)))
output_sw5 <- confint(eff_ws5,level=0.95) 
output_sw5
eff_wg5 <- contrast(emm_int_gw5, list(invertebrate = c(0, 0, 1, -1),
                                      bird = c(-1, 1, 0, 0),
                                      rodent = c(0, -1, 0, 1)))
output_gw5 <- confint(eff_wg5,level=0.95)
output_gw5

#Add large mammals
w5s.lg <- emmeans(ws5, ~ Treatment * PlotType, offset = log(1))
w5slg.eff <- contrast(w5s.lg, list(#mam=c(1,0,0,0,-1,0,0,0,0,0,0,0),
  large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w5slg.output <- confint(w5slg.eff,level=0.95)
#use mam2, averages the effects from habitat and control 
w5g.lg <- emmeans(wg5, ~ Treatment * PlotType, offset = log(1))
w5glg.eff <- contrast(w5g.lg, list(#mam=c(1,0,0,0,-1,0,0,0,0,0,0,0),
  large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w5glg.output <- confint(w5glg.eff,level=0.95)

#add these in to the output objects so it will be easier to make figures
w5slg.output$PlotType <- "Large Mammals"
output_sw5 <- rbind(output_sw5, w5slg.output)
w5glg.output$PlotType <- "Large Mammals"
output_gw5 <- rbind(output_gw5, w5glg.output)

#10 day model
#Need to create an exposure column (10 for all)
wseed$exposure10 <- 10
wgerm$exposure10 <- 10

mod_wsbayes <- bglmer(cbind(D10SW,failure10) ~ PlotType*Treatment  + (1|PlotID), 
                      family = binomial,
                      fixef.prior = normal(cov = diag(9,12)), #var = 9, sd = 3
                      control = glmerControl(optimizer="bobyqa"), 
                      data = wseed)
mod_wgbayes <- bglmer(cbind(D10GW,failure10) ~ PlotType*Treatment  + (1|PlotID), 
                      family = binomial,
                      fixef.prior = normal(cov = diag(9,12)), #var = 9, sd = 3
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
                                        rodent = c(0, -1, 0, 1)))

output_sw10 <- confint(eff_ws10,level=0.95) 
output_sw10
eff_wg10 <- contrast(emm_int_gw10, list(invertebrate = c(0, 0, 1, -1),
                                        bird = c(-1, 1, 0, 0),
                                        rodent = c(0, -1, 0, 1)))
output_gw10 <- confint(eff_wg10,level=0.95)
output_gw10

#Add large mammals
w10s.lg <- emmeans(mod_wsbayes, ~ Treatment * PlotType, offset = log(1))
w10slg.eff <- contrast(w10s.lg, list(#mam=c(1,0,0,0,-1,0,0,0,0,0,0,0),
  large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w10slg.output <- confint(w10slg.eff,,level=0.95)
#use mam2, averages the effects from habitat and control 
w10g.lg <- emmeans(mod_wgbayes, ~ Treatment * PlotType, offset = log(1))
w10glg.eff <- contrast(w10g.lg, list(#mam=c(1,0,0,0,-1,0,0,0,0,0,0,0),
  large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
w10glg.output <- confint(w10glg.eff,,level=0.95)

#add these in to the output objects so it will be easier to make figures
w10slg.output$PlotType <- "Large Mammals"
output_sw10 <- rbind(output_sw10, w10slg.output)
w10glg.output$PlotType <- "Large Mammals"
output_gw10 <- rbind(output_gw10, w10glg.output)


#### Both season analysis ####
#Make the column names the same and combine seasons into one df

herb.wet$seed.ratio5 <- herb.wet$D5SW/herb.wet$InitSeed
herb.wet$germ.ratio5 <- herb.wet$D5GW/herb.wet$InitGerm

names(herb.dry)
names(herb.wet)

dry.match <- herb.dry %>% dplyr::select(-c(X,obsid,id,date,observer,seeds.remaining,seeds.partial,
                                           germ.consumed,notes,scribe,entry,maxday,max,germ.partial)) %>%
  rename(PlotID=plot,Cluster=cluster,Treatment=treatment,TrialNo=trial,
         D5SW=seeds.intact,D5GW=no.germ.intact,
         seed.ratio5=seed.ratio,germ.ratio5=germ.ratio,FailSeed=failure,FailGerm=failuregerm,
         InitSeed=initial.seeds, InitGerm=initial.germ) 

dry.match$season <- "DRY"
dry.match$Cluster[dry.match$Cluster=="A"] <- "OPEN"
dry.match$Cluster[dry.match$Cluster=="B"] <- "CLOSED"
dry.match <- dry.match[,-6] #redundant to exposure column

wet.match <- herb.wet %>% dplyr::select(-c(InitDate,D1SW,D1SP,D1GW,D1GP,D2SW,D2SP,D2GW,D2GP,D3SW,D3SP,
                                           D3GW,D3GP,D4GW,D4GP,D4SW,D4GP,D8SW,D8SP,D8GW,D8GP,D10SW,D10SP,
                                           D10GW,D10GP,D4SP,D5SP,D5GP))
wet.match$season <- as.character("WET")
wet.match$Cluster <- as.character(wet.match$Cluster)
wet.match$Cluster[wet.match$Cluster=="O"] <- "OPEN"
wet.match$Cluster[wet.match$Cluster=="C"] <- "CLOSED"

wet.match$FailSeed <- wet.match$InitSeed-wet.match$D5SW
wet.match$FailGerm <- wet.match$InitGerm-wet.match$D5GW
wet.match$exposure <- 5

both.match <- rbind(wet.match,dry.match)
unique(both.match$PlotType)

both.match$PlotID <- as.factor(both.match$PlotID)
both.match$PlotType <- as.factor(both.match$PlotType)
both.match$Cluster <- as.factor(both.match$Cluster)
both.match$Treatment <- as.factor(both.match$Treatment)
both.match$Treatment <- relevel(both.match$Treatment, ref = "Control")
both.match$PlotType <- relevel(both.match$PlotType, ref = "control")
both.match$season <- as.factor(both.match$season)

both.seed <- both.match %>% dplyr::select(-c(D5GW,germ.ratio5,FailGerm,InitGerm)) %>%
  filter(TrialNo != 4)
both.germ <- both.match %>% dplyr::select(-c(D5SW,seed.ratio5,FailSeed,InitSeed)) %>% 
  filter(!is.na(germ.ratio5))
both.germ$FailGerm[both.germ$FailGerm<0]<-0 
both.seed$FailSeed[both.seed$FailSeed<0]<-0 

#Bayesian GLMER with the addition of season as a covariate
#First use template model builder to figure out dimensions
both_failed_seed <- glmmTMB(FailSeed ~ PlotType + Treatment + season + offset(log(exposure)) +
                              
                              PlotType:Treatment + PlotType:season + season:Treatment +
                              
                              (1|PlotID),
                            
                            family = poisson, data = both.seed)

dim(model.matrix(both_failed_seed)) #number of columns specifies the number of parameters


m1both_seed <- bglmer(cbind(D5SW,FailSeed) ~ PlotType + Treatment + season + offset(log(exposure)) +
                        PlotType:Treatment + PlotType:season + season:Treatment + (1|PlotID), 
                      family = binomial,
                      fixef.prior = normal(cov=diag(9,18)), 
                      control = glmerControl(optimizer="bobyqa"), 
                      data = both.seed)
m1both_germ <- bglmer(cbind(D5GW,FailGerm) ~ PlotType + Treatment + season + offset(log(exposure)) +
                        PlotType:Treatment + PlotType:season + season:Treatment + (1|PlotID), 
                      family = binomial,
                      fixef.prior = normal(cov=diag(9,18)), 
                      control = glmerControl(optimizer="bobyqa"), 
                      data = both.germ)

#emmeans
emm_int_sb <- emmeans(m1both_seed, ~ Treatment | PlotType | season, offset = log(1))
emm_int_sb
emm_int_gb <- emmeans(m1both_germ, ~ Treatment | PlotType | season, offset = log(1))
emm_int_gb

#General effect: C, B, F, RB
#Want to do absent - present
#July 2025: added a contrast for pooled effect of all cryptic taxa
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

#large mammal contrast should be absent - present; fence control - control control
#would need 12 numbers, all 0 except control control and fence control
#probably need to make a new object (eff)
#need to set new equation with emmeans function so you can make a contrast across both factors
#not sure what symbol to use instead of |
lgtest <- emmeans(m1both_seed, ~ Treatment * PlotType | season, offset = log(1))
lgtest.eff <- contrast(lgtest, list(#mam=c(1,0,0,0,-1,0,0,0,0,0,0,0),
  large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
lgmam.seed <- confint(lgtest.eff)
#use mam2, averages the effects from habitat and control 
lg.germ <- emmeans(m1both_germ, ~ Treatment * PlotType | season, offset = log(1))
lg.germ.eff <- contrast(lg.germ, list(#mam=c(1,0,0,0,-1,0,0,0,0,0,0,0),
  large_mammal =c(-0.5,0,0,0,1,0,0,0,-0.5,0,0,0)))
lg.germ.output <- confint(lg.germ.eff)

#add these in to the output objects so it will be easier to make figures
lgmam.seed$PlotType <- NA
output_sb <- rbind(output_sb, lgmam.seed)
lg.germ.output$PlotType <- NA
output_gb <- rbind(output_gb, lg.germ.output)


############################################################################################
### Analyses in Supplement ###

#### Rodent Minimum Number Alive ####

rodent.all <- read.csv("data/rodentsurvey.csv")
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
summary(rod.test)

#### Vegetation Structure ####

canopy <- read.csv("data/canopy dry 23.csv")
shrub <- read.csv("data/shrub dry 23.csv")

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
allveg <- full_join(allveg, plot.treat, by="plot")

#ANOVA
canopy.sub <- canopy %>% dplyr::select(plot,cover) %>% filter(plot %in% c(10,11,18,19,8,6,5,3,2)) %>% full_join(plot.treat,by="plot") 
ca.aov <- aov(cover ~ treatment, data=canopy.sub)
summary(ca.aov)
ca.tu <- TukeyHSD(ca.aov)
plot(ca.tu)

shrub.sub <- shrub %>% dplyr::select(plot,shrub) %>% filter(plot %in% c(10,11,18,19,8,6,5,3,2)) %>% full_join(plot.treat,by="plot") 
sh.aov <- aov(shrub ~ treatment, data=shrub.sub)
summary(sh.aov)
sh.tu <- TukeyHSD(sh.aov)
plot(sh.tu)

