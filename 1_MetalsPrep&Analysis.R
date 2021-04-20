# M Schrandt
# July 17, 2020

# Read in the metals data and analyze via ANCOVA
# Question: Does sex or collection area of Tripletail affect/relate
#           to different heavy metal concentrations in fish muscle?
#           Use TL of fish as a covariate in the model to account 
#            for fish size and ontogenetic dietary changes

library("tidyverse")
library("MASS") # MASS package for stepAIC function

# Read in data table
TTT <- read.csv("data/TTTmetals.csv")

summary(TTT)
str(TTT)
# Check minimum TL
min(TTT$TL)
# min is 409 mm TL; legal size is 457 mm TL; difference is 48 mm, or 1.89 inches
# for now let's keep all fish, regardless of size because we are taking size into account in the model
# the ecological questions are whether the amount of metals differ with area, sex, and fish size...it doesn't realy matter what size
# we can still discuss the consumption issue since we are looking at size and the vast majority of the samples are legal-size

#### Subsetting Data ####

subset.EMS = subset(TTT, Area == "E_MS_Sound")
subset.WMS = subset(TTT, Area == "W_MS_Sound")
subset.NMOB = subset(TTT, Area == "N_Mob_Bay")
subset.SMOB = subset(TTT, Area == "S_Mob_Bay")
subset.MS = subset(TTT, Area.Combined == "MissSound")
subset.AL = subset(TTT, Area.Combined == "MobBay")
subset.M = subset(TTT, Sex == "M")
subset.F = subset(TTT, Sex == "F")

#### Exploring Mean Values and SE ####
library("plyr")
library("ggplot2")
# Mercury Means
Hg_sum = ddply(TTT, c("Sex", "Area.Combined"), summarise,
               N    = length(TotalHg_ppm_ww),
               mean = mean(TotalHg_ppm_ww),
               sd   = sd(TotalHg_ppm_ww),
               se   = sd / sqrt(N)
)
Hg_sum

#Mercury Plot
# Error bars represent standard error of the mean
Hg_plot = ggplot(Hg_sum, aes(x=Area.Combined, y=mean, fill=Sex)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Total Mercury (mg/kg wet weight)") 
Hg_plot

# Arsenic Means
As_sum = ddply(TTT, c("Sex", "Area.Combined"), summarise,
               N    = length(As_ppm_ww),
               mean = mean(As_ppm_ww),
               sd   = sd(As_ppm_ww),
               se   = sd / sqrt(N)
)
As_sum
# Arsenic Plot
As_plot = ggplot(As_sum, aes(x=Area.Combined, y=mean, fill=Sex)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Arsenic (mg/kg wet weight)") 
As_plot

# Cadmium Means
Cd_sum = ddply(TTT, c("Sex", "Area.Combined"), summarise,
               N    = length(Cd_ppm_ww),
               mean = mean(Cd_ppm_ww),
               sd   = sd(Cd_ppm_ww),
               se   = sd / sqrt(N)
)
Cd_sum

# Cadmium Plot
Cd_plot = ggplot(Cd_sum, aes(x=Area.Combined, y=mean, fill=Sex)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Cadmium (mg/kg wet weight)") 
Cd_plot

# Lead Means
Pb_sum = ddply(TTT, c("Sex", "Area.Combined"), summarise,
               N    = length(Pb_ppm_ww),
               mean = mean(Pb_ppm_ww),
               sd   = sd(Pb_ppm_ww),
               se   = sd / sqrt(N)
)
Pb_sum

# Lead Plot
Pb_plot = ggplot(Pb_sum, aes(x=Area.Combined, y=mean, fill=Sex)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Lead (mg/kg wet weight)") 
Pb_plot

# Selenium Means
Se_sum = ddply(TTT, c("Sex", "Area.Combined"), summarise,
               N    = length(Se_ppm_ww),
               mean = mean(Se_ppm_ww),
               sd   = sd(Se_ppm_ww),
               se   = sd / sqrt(N)
)
Se_sum

# Selenium Plot
Se_plot = ggplot(Se_sum, aes(x=Area.Combined, y=mean, fill=Sex)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  ggtitle("Selenium (mg/kg wet weight)") 
Se_plot

#### Step-wise GLMs ####
# Metals are As, Cd, Hg, Pb, Se
# Covariate is TL - good practice to put control/covariate variables into equation before variables of interest
# Factors are Area, Sex

# Arsenic
# fit a regular glm first for all areas
As.glm = glm(As_ppm_ww ~ TL+Area+Sex, data=TTT)
summary(As.glm) #only TL comes up as significant
# use step-wise, backward selection to get rid of unnecessary factors in model
As.step=step(As.glm)
summary(As.step)
# step-wise procedure suggests only TL is significant, therefore fit that model below
As.glm2 = glm(As_ppm_ww ~ TL, data = TTT)
summary(As.glm2)
plot(As.glm2)
#get predicted values
pred.As=predict(As.glm2)
pred.As

# fit glm for combined areas
As.glm.comb = glm(As_ppm_ww ~ TL+Area.Combined+Sex, data=TTT) #TL is significant, area marginally
summary(As.glm.comb)
As.comb.step=step(As.glm.comb)
summary(As.comb.step)
# step-wise suggests including TL and area.combined so fit that model
As.glm.comb2 = glm(As_ppm_ww ~ TL+Area.Combined, data=TTT)
summary(As.glm.comb2)
plot(As.glm.comb2)
#get predicted values
pred.As.comb=predict(As.glm.comb2)
pred.As.comb

# Cadmium
# fit a regular glm first for all areas
Cd.glm = glm(Cd_ppm_ww ~ TL+Area+Sex, data=TTT)
summary(Cd.glm) #only areas are significant, only S Mob Bay
# use step-wise, backward selection to get rid of unnecessary factors in model
Cd.step=step(Cd.glm)
summary(Cd.step)
# step-wise procedure suggests keeping TL and area, therefore fit that model below
Cd.glm2 = glm(Cd_ppm_ww ~ TL+Area, data = TTT)
summary(Cd.glm2)
plot(Cd.glm2)
#get predicted values
pred.Cd=predict(Cd.glm2)
pred.Cd

# fit glm for combined areas
Cd.glm.comb = glm(Cd_ppm_ww ~ TL+Area.Combined+Sex, data=TTT) #area is significant, TL close
summary(Cd.glm.comb)
Cd.comb.step=step(Cd.glm.comb)
summary(Cd.comb.step)
# step-wise suggests including TL and area.combined so fit that model
Cd.glm.comb2 = glm(Cd_ppm_ww ~ TL+Area.Combined, data=TTT)
summary(Cd.glm.comb2)
plot(Cd.glm.comb2)
#get predicted values
pred.Cd.comb=predict(Cd.glm.comb2)
pred.Cd.comb

# Total Mercury
# fit a regular glm first for all areas
Hg.glm = glm(TotalHg_ppm_ww ~ TL+Area+Sex, data=TTT)
summary(Hg.glm) #only TL significant
# use step-wise, backward selection to get rid of unnecessary factors in model
Hg.step=step(Hg.glm)
summary(Hg.step)
# step-wise procedure suggests just keeping TL, therefore fit that model below
Hg.glm2 = glm(TotalHg_ppm_ww ~ TL+Area, data = TTT)
summary(Hg.glm2)
plot(Hg.glm2)
#get predicted values
pred.Hg=predict(Hg.glm2)
pred.Hg

# fit glm for combined areas
Hg.glm.comb = glm(TotalHg_ppm_ww ~ TL+Area.Combined+Sex, data=TTT) #TL significant
summary(Hg.glm.comb)
Hg.comb.step=step(Hg.glm.comb)
summary(Hg.comb.step)
# step-wise suggests including TL, so fit that model
Hg.glm.comb2 = glm(TotalHg_ppm_ww ~ TL, data=TTT)
summary(Hg.glm.comb2) #TL significant
plot(Hg.glm.comb2)
#get predicted values
pred.Hg.comb=predict(Hg.glm.comb2)
pred.Hg.comb

# Lead
# fit a regular glm first for all areas
Pb.glm = glm(Pb_ppm_ww ~ TL+Area+Sex, data=TTT)
summary(Pb.glm) #nothing technically significant
# use step-wise, backward selection to get rid of unnecessary factors in model
Pb.step=step(Pb.glm)
summary(Pb.step)
# step-wise procedure suggests keeping TL, therefore fit that model below
Pb.glm2 = glm(Pb_ppm_ww ~ TL, data = TTT)
summary(Pb.glm2)
plot(Pb.glm2)
#get predicted values
pred.Pb=predict(Pb.glm2)
pred.Pb

# fit glm for combined areas
Pb.glm.comb = glm(Pb_ppm_ww ~ TL+Area.Combined+Sex, data=TTT) #nothing technically significant
summary(Pb.glm.comb)
Pb.comb.step=step(Pb.glm.comb)
summary(Pb.comb.step)
# step-wise suggests including TL and area.combined so fit that model
Pb.glm.comb2 = glm(Pb_ppm_ww ~ TL+Area.Combined, data=TTT)
summary(Pb.glm.comb2) #nothing significant
plot(Pb.glm.comb2)
#get predicted values
pred.Pb.comb=predict(Pb.glm.comb2)
pred.Pb.comb

# Selenium
# fit a regular glm first for all areas
Se.glm = glm(Se_ppm_ww ~ TL+Area+Sex, data=TTT)
summary(Se.glm) #nothing technically significant
# use step-wise, backward selection to get rid of unnecessary factors in model
Se.step=step(Se.glm)
summary(Se.step)
# step-wise procedure suggests keeping TL and sex, therefore fit that model below
Se.glm2 = glm(Se_ppm_ww ~ TL+Sex, data = TTT)
summary(Se.glm2)
plot(Se.glm2)
#get predicted values
pred.Se=predict(Se.glm2)
pred.Se

# fit glm for combined areas
Se.glm.comb = glm(Se_ppm_ww ~ TL+Area.Combined+Sex, data=TTT) #nothing technically significant
summary(Se.glm.comb)
Se.comb.step=step(Se.glm.comb)
summary(Se.comb.step)
# step-wise suggests including TL and Sex so fit that model
Se.glm.comb2 = glm(Se_ppm_ww ~ TL+Sex, data=TTT)
summary(Se.glm.comb2) #TL and sex significant
#get predicted values
pred.Se.comb=predict(Se.glm.comb2)
pred.Se.comb
plot(Se.glm.comb2)

#### Figures ####
# since we have a continuous covariate that is almost always sig., it's helpful to plot TL on x-axis
# plot concentration on y-axis
# then put different lines for each other significant factor (e.g., Area.combined for arsenic)
# the figures below are a start, but we need to look more into how to properly predict values and plot those

library(ggplot2)
library("gridExtra") #allows basically par mfrow function with ggplot2
# Arsenic
As.plot=ggplot(TTT, aes(y = As_ppm_ww, x = TL, group = Area.Combined, color=Area.Combined,
                        shape = Area.Combined)) + 
  geom_point() + geom_smooth(method = "glm")
As.plot

# Cadmium
Cd.plot=ggplot(TTT, aes(y = Cd_ppm_ww, x = TL, group = Area.Combined, color=Area.Combined,
                        shape = Area.Combined)) + 
  geom_point() + geom_smooth(method = "glm")
Cd.plot

# Total Mercury
Hg.plot=ggplot(TTT, aes(y = TotalHg_ppm_ww, x = TL, group = Area.Combined, color=Area.Combined,
                        shape = Area.Combined)) + 
  geom_point() + geom_smooth(method = "glm")
Hg.plot

# Lead
Pb.plot=ggplot(TTT, aes(y = Pb_ppm_ww, x = TL, group = Area.Combined, color=Area.Combined,
                        shape = Area.Combined)) + 
  geom_point() + geom_smooth(method = "glm")
Pb.plot

# Selenium
Se.plot=ggplot(TTT, aes(y = Se_ppm_ww, x = TL, group = Sex, color=Sex,
                        shape = Sex)) + 
  geom_point() + geom_smooth(method = "glm")
Se.plot

# make a panel graph essentially
grid.arrange(As.plot, Cd.plot, Hg.plot,
             Pb.plot, Se.plot, ncol=3)