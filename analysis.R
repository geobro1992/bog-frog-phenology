library(dplyr)
library(tidyverse)
library(lme4)
library(qpcR) # akaike.weights()
library(chron) # used for re-converting scaled dates and times
library(ordinal)
library(ggplot2)
library(ggeffects)
library(ggpubr)
library(readxl)
library(optimx)

dat <- read.csv("bf_dat.csv")

dat$y = as.numeric(dat$RANOKA)
dat$RANCLA = as.numeric(dat$RANCLA)
dat$RANCLA_Index = as.factor(dat$RANCLA_Index)
dat$RANOKA_Index = as.factor(dat$RANOKA_Index)
dat$WindScale <- as.factor(dat$WindScale)
dat$WindScale = recode(dat$WindScale, `3` = "2")
dat$LightConditions <- as.factor(dat$LightConditions)
dat$SkyCode <- as.factor(dat$SkyCode)
dat$SiteID <- as.factor(dat$SiteID)
dat$Rain = as.factor(dat$Rain)
dat$temp = scale(dat$Temp)
dat$RH = scale(dat$RH)

dat = dat %>%
  group_by(SiteID) %>%
  filter(sum(y) > 0) %>%
  filter(Year > 2006) %>%
  filter(Year < 2020) %>%
  droplevels() %>%
  ungroup()

length(unique(dat$SiteID))
table(dat$y)
table(dat$RANCLA)
table(dat$Year)

table(dat$y, dat$Year)
table(dat$RANCLA, dat$Year)


m1 = glmer(y ~ LightConditions + SkyCode + RH + Temp + time + date + I(date^2) + WindScale + (1|Year), 
            data = dat, family = binomial)

drop1(m1, test = "Chisq")


m1 = clm(RANOKA_Index ~ LightConditions + SkyCode + RH + Temp + time + date + I(date^2) + WindScale, 
           data = dat, family = binomial)

drop1(m1, test = "Chisq")


# RANCLA

m1 = glmer(RANCLA ~ LightConditions + SkyCode + RH + Temp + time + date + I(date^2) + WindScale + (1|Year), 
           data = dat, family = binomial)

drop1(m1, test = "Chisq")

m1 = clm(RANCLA_Index ~ LightConditions + SkyCode + RH + Temp + time + date + I(date^2) + WindScale, 
         data = dat, family = binomial)

drop1(m1, test = "Chisq")

##############
# species comparison

dat2 <- gather(dat,
               key = "taxa",
               value = "y",
               c(RANOKA, RANCLA))

dat2$taxa = as.factor(dat2$taxa)

m1 = glmer(y ~ LightConditions*taxa + SkyCode*taxa + RH*taxa + Temp*taxa + time*taxa + date*taxa + I(date^2)*taxa + WindScale*taxa + (1|Year), 
            data = dat2, family = binomial, control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))

summary(m1)
drop1(m1, test = "Chisq")

model0 = glm(y ~ x, family = binomial())
anova(model0,model, test = 'LRT')



##################
# Index

dat3 <- gather(dat,
                key = "taxa",
                value = "y",
                c(RANOKA_Index, RANCLA_Index))

dat3$taxa = as.factor(dat3$taxa)
levels(dat3$taxa) = c("RANCLA", "RANOKA")
dat3$y = as.factor(dat3$y)

m1 = clm(y ~ LightConditions*taxa + SkyCode*taxa + RH*taxa + Temp*taxa + time*taxa + date*taxa + I(date^2)*taxa + WindScale*taxa, 
         data = dat3, family = binomial)

drop1(m1, test = "Chisq")


