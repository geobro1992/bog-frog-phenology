# required libraries
library(tidyverse)
library(lme4)
library(qpcR)
library(ordinal)
library(ggplot2)
library(ggeffects)
library(ggpubr)
library(lubridate)

# read in data
dat <- read.csv("bf_dat.csv")

# format variables
dat$SiteID <- as.factor(dat$SiteID)
dat$RH = scale(dat$RH)

dat = dat %>%
  mutate(SurveyDate = as.Date(SurveyDate, format = "%m/%d/%Y")) %>%
  mutate(survdate = as.Date(survdate, format = "%m/%d/%Y")) %>%
  mutate(Year = year(SurveyDate))


#####################
# logistic regression
#####################

# gather presence/absence records for both species
dat2 <- gather(dat,
               key = "taxa",
               value = "y",
               c(RANOKA, RANCLA))

# format variables
dat2$taxa = as.factor(dat2$taxa)
dat2$Year = as.factor(dat2$Year)


# best suppported model (temperature and wind speed omitted)
m1.top = glmer(y ~ date + I(date^2) + time + taxa*RH + taxa*LightConditions + taxa*SkyCode + 
                 (1|SiteID) + (1|Year), data = dat2, family = binomial, 
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(m1.top)

# full model
m1.full = glmer(y ~ date + I(date^2) + taxa*time + taxa*RH + taxa*LightConditions + taxa*SkyCode + taxa*WindScale + taxa*Temp +
                  (1|SiteID) + (1|Year), data = dat2, family = binomial, 
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(m1.full)

# compare model support with AIC
AIC(m1.top, m1.full)


#######################
# predictions and plots

# Date
preds = ggpredict(m1.top, terms = c("date [all]", "taxa"), ci.lvl = 0.95, type = "fixed")

p1 = plot(preds) + 
  ylim(0,1) +
  theme_Publication() + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(breaks=c(-0.15, 0.15, 0.45, 0.75, 1.05), labels = c("May", "Jun", "Jul", "Aug", "Sep"))+
  theme(legend.title = element_blank()) +
  ggtitle("") + xlab("Date") + ylab("Detection Probability")


# Time
preds = ggpredict(m1.top, terms = c("time [all]", "taxa"), ci.lvl = 0.95, type = "fixed")

p2 = plot(preds) + 
  ylim(0,1) +
  theme_Publication() + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.title = element_blank()) +
  scale_x_continuous(breaks=c(-2.25, -0.75, 0.75), labels = c("21:00", "23:00", "01:00"))+
  ggtitle("") + xlab("Time") + ylab("Detection Probability")


# Humidity
preds = ggpredict(m1.top, terms = c("RH [all]", "taxa"), ci.lvl = 0.95, type = "fixed")

p3 = plot(preds) + 
  ylim(0,1) +
  theme_Publication() + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_x_continuous(breaks=c(-5.63, -3.54903074, -1.37860756, 0.79181562), labels = c("30", "50", "70", "90"))+
  ggtitle("") + xlab("Relative Humidity") + ylab("Detection Probability")


# Sky Code
preds = ggpredict(m1.top, terms = c("SkyCode [all]", "taxa"), ci.lvl = 0.95, type = "fixed")

p4 = plot(preds) + 
  ylim(0,1) +
  theme_Publication() + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_x_continuous(breaks=c(0, 1, 2), labels = c("0", "1", "2"))+
  ggtitle("") + xlab("Cloud Cover") + ylab("Detection Probability")


# Light Conditions
preds = ggpredict(m1.top, terms = c("LightConditions [all]", "taxa"), ci.lvl = 0.95, type = "fixed")

p5 = plot(preds) + 
  ylim(0,1) +
  theme_Publication() + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.title = element_blank(), legend.position = "none") +
  ggtitle("") + xlab("Lunar Illumination") + ylab("Detection Probability")


f1 = ggarrange(p1, p2, 
               nrow = 1, ncol = 2, 
               common.legend = T, legend = "top", 
               labels = c("A", "B")) +
  theme(plot.margin = margin(0.1,0.1,0.1,0.1, "mm"))
ggsave("Figure1.pdf", f1, device = cairo_pdf, width = 8, height = 4, dpi = 600)


pleg <- plot(preds)+ ggtitle("") +
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  16),
        legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=8)))

f2 = ggarrange(p3, p5, p4, pleg,
               nrow = 2, ncol = 2, 
               common.legend = F, 
               labels = c("A", "B", "C")) +
  theme(plot.margin = margin(0.1,0.1,0.1,0.1, "mm"))
ggsave("Figure2.pdf", f2, device = cairo_pdf, width = 8, height = 8, dpi = 600)

####################
# ordinal regression
####################

# gather calling index data for both species
dat.I <- gather(dat,
                key = "taxa",
                value = "y",
                c(RANOKA_Index, RANCLA_Index))

# format variables
dat.I$taxa = as.factor(dat.I$taxa)
dat.I$y = as.factor(dat.I$y)
dat.I$Year = as.factor(dat.I$Year)


# top model
m2.top = clm(y ~ date + I(date^2) + time + taxa*RH + taxa*LightConditions + taxa*SkyCode + taxa*WindScale + Temp, 
             data = dat.I, Hess = T)
summary(m2.top)

m2.full = clm(y ~ date + I(date^2) + taxa*time + taxa*RH + taxa*LightConditions + taxa*SkyCode + taxa*WindScale + taxa*Temp, 
              data = dat.I, Hess = T)
summary(m2.full)

AIC(m2.top, m2.full)

# plotting predictions
preds = ggpredict(m2.top, terms = c("date [all]", "taxa"), ci.lvl = 0.95, type = "fixed")
plot(preds)
preds = ggpredict(m2.top, terms = c("time [all]", "taxa"), ci.lvl = 0.95, type = "fixed")
plot(preds)
preds = ggpredict(m2.top, terms = c("RH [all]", "taxa"), ci.lvl = 0.95, type = "fixed")
plot(preds)
preds = ggpredict(m2.top, terms = c("LightConditions [all]", "taxa"), ci.lvl = 0.95, type = "fixed")
plot(preds)
preds = ggpredict(m2.top, terms = c("SkyCode [all]", "taxa"), ci.lvl = 0.95, type = "fixed")
plot(preds)
preds = ggpredict(m2.top, terms = c("WindScale [all]", "taxa"), ci.lvl = 0.95, type = "fixed")
plot(preds)
preds = ggpredict(m2.top, terms = c("Temp [all]", "taxa"), ci.lvl = 0.95, type = "fixed")
plot(preds)
