# Taylor van Doren
# May 6, 2024

# REVISED
# 1918 influenza pandemic BMI analyses

# load libraries
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(tidyverse)
library(reshape2)
library(lmtest)
library(pscl)

# setwd and load data
setwd("/Users/taylorvandoren/Dropbox/Papers/Body size papers/Data")
data <- read.csv("DeathRecords.csv")

data <- data %>% select(
  file_name,
  judicial_district,
  sex,
  race_cat,
  birth_season,
  death_season,
  death_year,
  age_group,
  reassign_age,
  birthplace_cat,
  primary_death_gencat,
  primary_death_duration_years,
  secondary_death_gencat,
  secondary_death_duration_years,
  height_in,
  weight_lb
) %>% rename(
  index = file_name,
  district = judicial_district,
  sex = sex,
  race_cat = race_cat,
  birth_season = birth_season,
  death_season = death_season,
  yod = death_year,
  age_group = age_group,
  age = reassign_age,
  bp_cat = birthplace_cat,
  cod = primary_death_gencat,
  cod_years = primary_death_duration_years,
  sec_cod = secondary_death_gencat,
  sec_cod_years = secondary_death_duration_years,
  height = height_in,
  weight = weight_lb
)


data$age <- as.numeric(data$age)
data$height <- as.numeric(data$height)
data$weight <- as.numeric(data$weight)
data$yod <- as.factor(data$yod)

# calculate BMI
meter <- 39.37
kg <- 2.20462

data$height.m <- data$height / meter
data$weight.kg <- data$weight / kg
data$bmi <- data$weight.kg / (data$height.m^2)

# categorize data
data$cod[which(data$cod == "pneumonia")] <- "p&i"
data$cod[which(data$cod == "influenza")] <- "p&i"
data$cod[which(data$cod == "tuberculosis")] <- "tb"
data$cod[which(data$cod == "other")] <- "other causes"
data$cod[which(data$cod == "nutritional disorder")] <- "nutritional disorders"

data <- data %>% filter(
  age >= 18,
  race_cat != "",
  !is.na(bmi)
)

tb <- subset(data, data$cod == "tb")
flu <- subset(data, data$cod == "p&i")
acc <- subset(data, data$cod == "violence, accidents")
circ <- subset(data, data$cod == "circulatory system")
nerv <- subset(data, data$cod == "nervous system")

data <- rbind(tb, flu, acc, circ, nerv)

# get ready to analyze
data$logcod <- ifelse(data$cod == "p&i", 1, 0)
data$sex <- as.factor(data$sex)
data$race_cat <- as.factor(data$race_cat)
data$district <- as.factor(data$district)
data$sex <- ifelse(data$sex == "f", "Female", "Male")

pand <- subset(data, data$yod == 1918 | data$yod == 1919 | data$yod == 1920)
nonpand <- subset(data, data$yod == 1917 | data$yod == 1921 | data$yod == 1922 | data$yod == 1923 | 
                    data$yod == 1924 | data$yod == 1925)


###################
# BASELINE MODELS #
###################

### LINEAR

## PANDEMIC
base.lin.pan <- glm(logcod ~ bmi, data = pand, family = "binomial") # fit model
exp(summary(base.lin.pan)$coefficients) # exponentiate coefficients
pR2(base.lin.pan) # McF = 0.0011

## NON-PANDEMIC
base.lin.nonpan <- glm(logcod ~ bmi, data = nonpand, family = "binomial")
exp(summary(base.lin.nonpan)$coefficients)
pR2(base.lin.nonpan) # McF = 0.00012

### QUADRATIC

## PANDEMIC
base.quad.pan <- glm(logcod ~ bmi + I(bmi^2), data = pand, family = "binomial")
exp(summary(base.quad.pan)$coefficients)
pR2(base.quad.pan) # McF = 0.0032

# optimize model
f <- function(x) { base.quad.pan$coefficients[1] + base.quad.pan$coefficients[2]*x + base.quad.pan$coefficients[3]*(x^2) }
optimize(f, interval = c(10,40), maximum = T) # vertex = 29.23128

## NONPANDEMIC
base.quad.nonpan <- glm(logcod ~ bmi + I(bmi^2), data = nonpand, family = "binomial")
exp(summary(base.quad.nonpan)$coefficients)
pR2(base.quad.nonpan) # McF = 0.0013
f <- function(x) { base.quad.nonpan$coefficients[1] + base.quad.nonpan$coefficients[2]*x + base.quad.pan$coefficients[3]*(x^2) }
optimize(f, interval = c(10,40), maximum = T) # vertex = 23.89542

### CUBIC

## PANDEMIC
base.cub.pan <- glm(logcod ~ bmi + I(bmi^2) + I(bmi^3), data = pand, family = "binomial")
exp(summary(base.cub.pan)$coefficients)
pR2(base.cub.pan) # McF = 0.0097
probs.cub.pan <- predict(base.cub.pan, pand, type = "response")
test <- data.frame(pand, probs.cub.pan)
ggplot(test, aes(bmi, probs.cub.pan)) +
  geom_point(alpha = 0.5)
f <- function(x) { 
  base.cub.pan$coefficients[1] + base.cub.pan$coefficients[2]*x +
    base.cub.pan$coefficients[3]*(x^2) + base.cub.pan$coefficients[4]*(x^3)
}
optimize(f, interval = c(0,30), maximum = TRUE) # vertex = 24.40079
optimize(f, interval = c(25, 45)) # vertex = 35.33805

## NONPANDEMIC
base.cub.nonpan <- glm(logcod ~ bmi + I(bmi^2) + I(bmi^3), data = nonpand, family = "binomial")
exp(summary(base.cub.nonpan)$coefficients)
pR2(base.cub.nonpan) # McF = 0.0019
probs.cub.nonpan <- predict(base.cub.nonpan, nonpand, type = "response")
test <- data.frame(nonpand, probs.cub.nonpan)
ggplot(test, aes(bmi, probs.cub.nonpan)) +
  geom_point(alpha = 0.5)
f <- function(x) {
  base.cub.nonpan$coefficients[1] + base.cub.nonpan$coefficients[2]*x +
    base.cub.nonpan$coefficients[3]*(x^2) + base.cub.nonpan$coefficients[4]*(x^3)
}
optimize(f, interval = c(10,30), maximum = TRUE) # vertex = 23.88139
optimize(f, interval = c(30,45)) # vertex = 37.6605

f# PANDEMIC
# A list of fitted coefficients and p-values of baseline models
list(summary(base.lin.pan)$coefficients, 
     summary(base.quad.pan)$coefficients, 
     summary(base.cub.pan)$coefficients)

# A list of exponentiated coefficients of baseline models
list(exp(summary(base.lin.pan)$coefficients), 
     exp(summary(base.quad.pan)$coefficients), 
     exp(summary(base.cub.pan)$coefficients))

# NONPANDEMIC
# A list of fitted coefficients and p-values of baseline models
list(summary(base.lin.nonpan)$coefficients, 
     summary(base.quad.nonpan)$coefficients, 
     summary(base.cub.nonpan)$coefficients)

# A list of exponentiated coefficients of baseline models
list(exp(summary(base.lin.nonpan)$coefficients), 
     exp(summary(base.quad.nonpan)$coefficients), 
     exp(summary(base.cub.nonpan)$coefficients))


#######################
# MULTIVARIATE MODELS #
#######################

# varibales to be fit: bmi + district + age + sex + race_cat
# we are varying the fit of BMI

### LINEAR

## PANDEMIC
set.seed(2) # keep things consistent
lin.pan.obs <- sample(1:nrow(pand), 0.8*nrow(pand), replace = FALSE) # split into train and test
lin.pan.train <- pand[lin.pan.obs,] # training data set for model fitting
lin.pan.test <- pand[-lin.pan.obs,] # testing data set for which to predict values
lin.pan.mod <- glm(logcod ~ bmi + district + age + sex + race_cat, data = lin.pan.train, family = "binomial") # model fit with training data
probs <- predict.glm(lin.pan.mod, lin.pan.test, type = "response") # predict probability of P&I death and add to testing dataset
lin.pan.test$probs <- probs
test.pred <- rep("Other", length(probs)) # set up default vector with test predictions
test.pred[probs > 0.5] <- "P&I" # every prediction > 0.5 will be coded as a probable P&I death
(table(test.pred, lin.pan.test$logcod)[1] + table(test.pred, lin.pan.test$logcod)[4])/sum(table(test.pred, lin.pan.test$logcod)) # predicts correctly 68.2% of the time
total.lin.pan <- glm(logcod ~ bmi + district + age + sex + race_cat, data = pand, family = "binomial") # predict values for all data
pand$probs <- predict(total.lin.pan, pand, type = "response") # add predicted values to pandemic data
ggplot(pand, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(shape = 1, alpha = 0.35) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.175) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1)

# get coefficients and exponentiate
lin.pan.mod <- glm(logcod ~ bmi + district + age + sex + race_cat, data = pand, family = "binomial")
summary(lin.pan.mod)
exp(summary(lin.pan.mod)$coefficients)
pR2(lin.pan.mod) # McF = 0.1933817

# optimize model
# we don't have to optimize a linear model :-)


## NONPANDEMIC
set.seed(2)
lin.nonpan.obs <- sample(1:nrow(nonpand), 0.8*nrow(nonpand), replace = FALSE)
lin.nonpan.train <- nonpand[lin.nonpan.obs,]
lin.nonpan.test <- nonpand[-lin.nonpan.obs,]
lin.nonpan.mod <- glm(logcod ~ bmi + district + age + sex + race_cat, data = lin.nonpan.train, family = "binomial")
probs <- predict.glm(lin.nonpan.mod, lin.nonpan.test, type = "response")
lin.nonpan.test$probs <- probs
test.pred <- rep("Other", length(probs))
test.pred[probs>0.5] <- "P&I"
table(test.pred, lin.nonpan.test$logcod) # this model will NEVER predict a flu death
(28+1)/sum(table(test.pred, lin.nonpan.test$logcod)) # 9.5%
total.lin.nonpan <- glm(logcod ~ bmi + district + age + sex + race_cat, data = nonpand, family = "binomial")
nonpand$probs <- predict(total.lin.nonpan, nonpand, type = "response")
ggplot(nonpand, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(shape = 1, alpha = 0.35) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.175) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1)

# get coefficients and exponentiate
lin.nonpan.mod <- glm(logcod ~ bmi + district + age + sex + race_cat, data = nonpand, family = "binomial")
summary(lin.nonpan.mod)
exp(summary(lin.nonpan.mod)$coefficients)
pR2(lin.nonpan.mod) # McF = 0.06804145

# optimize model
# we don't have to optimize a linear model :-)

### QUADRATIC

## PANDEMIC
set.seed(2) # keep things consistent
quad.pan.obs <- sample(1:nrow(pand), 0.8*nrow(pand), replace = FALSE) # split into train and test
quad.pan.train <- pand[quad.pan.obs,] # training data set for model fitting
quad.pan.test <- pand[-quad.pan.obs,] # testing data set for which to predict values
quad.pan.mod <- glm(logcod ~ bmi + I(bmi^2) + district + age + sex + race_cat, data = quad.pan.train, family = "binomial") # model fit with training data
probs <- predict.glm(quad.pan.mod, quad.pan.test, type = "response") # predict probability of P&I death and add to testing dataset
quad.pan.test$probs <- probs
quad.pred <- rep("Other", length(probs)) # set up default vector with test predictions
quad.pred[probs > 0.5] <- "P&I" # every prediction > 0.5 will be coded as a probable P&I death
(table(quad.pred, quad.pan.test$logcod)[1] + table(quad.pred, quad.pan.test$logcod)[4])/sum(table(quad.pred, quad.pan.test$logcod)) # predicts correctly 70.5% of the time
total.quad.pan <- glm(logcod ~ bmi + I(bmi^2) + district + age + sex + race_cat, data = pand, family = "binomial")
pand$probs <- predict(total.quad.pan, pand, type = "response")
ggplot(pand, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(shape = 1, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 0.5, alpha = 0.175) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1)

# get coefficients and exponentiate
quad.pan.mod <- glm(logcod ~ bmi + I(bmi^2) + district + age + sex + race_cat, data = pand, family = "binomial")
summary(quad.pan.mod)
exp(summary(quad.pan.mod)$coefficients)
pR2(quad.pan.mod) # McF = 0.1982829

# optimize model
quad.function <- function(x) {
  quad.pan.mod$coefficients[1] + quad.pan.mod$coefficients[2]*x + quad.pan.mod$coefficients[3]*(x^2) + 
    quad.pan.mod$coefficients[8]*x
}
optimize(quad.function, interval = c(20,40), maximum = T) # vertex = 30.51238

# hypothetical individual: everything the same except Indigeneity
# bmi = vertex, district = 2, age = 30, sex = F, race_cat = Indigenous
# versus
# bmi = vertex, district = 2, age = 30, sex = F, race_cat = non-Indigenous
exp(quad.pan.mod$coefficients[1] + quad.pan.mod$coefficients[2]*31.5 + quad.pan.mod$coefficients[3]*(31.5^2) +
  quad.pan.mod$coefficients[4] + quad.pan.mod$coefficients[8]*30) /
  (1 + exp(quad.pan.mod$coefficients[1] + quad.pan.mod$coefficients[2]*31.5 + quad.pan.mod$coefficients[3]*(31.5^2) +
             quad.pan.mod$coefficients[4] + quad.pan.mod$coefficients[8]*30)) # AKN probability of P&I = 89%

exp(quad.pan.mod$coefficients[1] + quad.pan.mod$coefficients[2]*31.5 + quad.pan.mod$coefficients[3]*(31.5^2) +
      quad.pan.mod$coefficients[4] + quad.pan.mod$coefficients[8]*30 + quad.pan.mod$coefficients[10]) /
  (1 + exp(quad.pan.mod$coefficients[1] + quad.pan.mod$coefficients[2]*31.5 + quad.pan.mod$coefficients[3]*(31.5^2) +
             quad.pan.mod$coefficients[4] + quad.pan.mod$coefficients[8]*30 + quad.pan.mod$coefficients[10]))


## NONPANDEMIC
set.seed(2)
quad.nonpan.obs <- sample(1:nrow(nonpand), 0.8*nrow(nonpand), replace = FALSE)
quad.nonpan.train <- nonpand[quad.nonpan.obs,]
quad.nonpan.test <- nonpand[-quad.nonpan.obs,]
quad.nonpan.mod <- glm(logcod ~ bmi + I(bmi^2) + district + age + sex + race_cat, data = quad.nonpan.train, family = "binomial")
probs <- predict.glm(quad.nonpan.mod, quad.nonpan.test, type = "response")
quad.nonpan.test$probs <- probs
quad.pred <- rep("Other", length(probs))
quad.pred[probs>0.5] <- "P&I"
table(quad.pred, quad.nonpan.test$logcod) 
(28+1)/sum(table(quad.pred, quad.nonpan.test$logcod)) # 9.5%
total.quad.nonpan <- glm(logcod ~ bmi + I(bmi^2) + district + age + sex + race_cat, data = nonpand, family = "binomial")
nonpand$probs <- predict(total.quad.nonpan, nonpand, type = "response")
ggplot(nonpand, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(shape = 1, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 0.5, alpha = 0.175) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1)

# get coefficients and exponentiate
quad.nonpan.mod <- glm(logcod ~ bmi + I(bmi^2) + district + age + sex + race_cat, data = nonpand, family = "binomial")
summary(quad.nonpan.mod)
exp(summary(quad.nonpan.mod)$coefficients)
pR2(quad.nonpan.mod) # McF = 0.06968274

# optimize model
f <- function(x) {
  quad.nonpan.mod$coefficients[1] + quad.nonpan.mod$coefficients[2]*x + quad.nonpan.mod$coefficients[3]*(x^2) +
    quad.nonpan.mod$coefficients[8]
}
optimize(f, interval = c(15,40), maximum = T) # vertex = 26.01624


### CUBIC

## PANDEMIC
set.seed(2) # keep things consistent
cub.pan.obs <- sample(1:nrow(pand), 0.8*nrow(pand), replace = FALSE) # split into train and test
cub.pan.train <- pand[cub.pan.obs,] # training data set for model fitting
cub.pan.test <- pand[-cub.pan.obs,] # testing data set for which to predict values
cub.pan.mod <- glm(logcod ~ bmi + I(bmi^2) + I(bmi^3) + district + age + sex + race_cat, data = cub.pan.train, family = "binomial") # model fit with training data
probs <- predict.glm(cub.pan.mod, cub.pan.test, type = "response") # predict probability of P&I death and add to testing dataset
cub.pan.test$probs <- probs
cub.pred <- rep("Other", length(probs)) # set up default vector with test predictions
cub.pred[probs > 0.5] <- "P&I" # every prediction > 0.5 will be coded as a probable P&I death
(table(cub.pred, cub.pan.test$logcod)[1] + table(cub.pred, cub.pan.test$logcod)[4])/sum(table(cub.pred, cub.pan.test$logcod)) # predicts correctly 70.5% of the time
total.cub.pan <- glm(logcod ~ bmi + I(bmi^2) + I(bmi^3) + district + age + sex + race_cat, data = pand, family = "binomial")
pand$probs <- predict(total.cub.pan, pand, type = "response")
ggplot(pand, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(shape = 1, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), size = 0.5, alpha = 0.175) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1)

# get coefficients and exponentiate
cub.pan.mod <- glm(logcod ~ bmi + I(bmi^2) +I(bmi^3) + district + age + sex + race_cat, data = pand, family = "binomial")
summary(cub.pan.mod)
exp(summary(cub.pan.mod)$coefficients)
pR2(cub.pan.mod) # McF = 0.2109958

# optimize
f <- function(x) {
  cub.pan.mod$coefficients[1] + cub.pan.mod$coefficients[2]*x + cub.pan.mod$coefficients[3]*(x^2) + cub.pan.mod$coefficients[4]*(x^3) +
    cub.pan.mod$coefficients[9]
}
optimize(f, interval = c(10,30), maximum = T) # vertex = 26.01623
optimize(f, interval = c(30,45)) # vertex = 35.25268

## NONPANDEMIC
set.seed(2)
cub.nonpan.obs <- sample(1:nrow(nonpand), 0.8*nrow(nonpand), replace = FALSE)
cub.nonpan.train <- nonpand[cub.nonpan.obs,]
cub.nonpan.test <- nonpand[-cub.nonpan.obs,]
cub.nonpan.mod <- glm(logcod ~ bmi + I(bmi^2) + I(bmi^3) + district + age + sex + race_cat, data = cub.nonpan.train, family = "binomial")
probs <- predict.glm(cub.nonpan.mod, cub.nonpan.test, type = "response")
cub.nonpan.test$probs <- probs
cub.pred <- rep("Other", length(probs))
cub.pred[probs>0.5] <- "P&I"
table(cub.pred, cub.nonpan.test$logcod) 
(28+1)/(sum(table(cub.pred, cub.nonpan.test$logcod))) # 9.5%
total.cub.nonpan <- glm(logcod ~ bmi + I(bmi^2) + I(bmi^3) + district + age + sex + race_cat, data = nonpand, family = "binomial")
nonpand$probs <- predict(total.cub.nonpan, nonpand, type = "response")
ggplot(nonpand, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(shape = 1, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), size = 0.5, alpha = 0.175) +
  scale_color_viridis_d(direction = -1) +
  scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1)

# get coefficients and exponentiate
cub.nonpan.mod <- glm(logcod ~ bmi + I(bmi^2) +I(bmi^3) + district + age + sex + race_cat, data = nonpand, family = "binomial")
summary(cub.nonpan.mod)
exp(summary(cub.nonpan.mod)$coefficients)
pR2(cub.nonpan.mod)

# optimize model
f <- function(x) {
  cub.nonpan.mod$coefficients[1] + cub.nonpan.mod$coefficients[2]*x + cub.nonpan.mod$coefficients[3]*(x^2) + cub.nonpan.mod$coefficients[4]*(x^3) +
    cub.nonpan.mod$coefficients[9]
}
optimize(f, interval = c(10,30), maximum = T) # vertex =25.2786
optimize(f, interval = c(30,45)) # vertex = 42.33692

######################################
# COMPARE MULTIVARIATES TO BASELINES #
######################################

# pandemic
lrtest(base.lin.pan, total.lin.pan) # p < 0.001
lrtest(base.quad.pan, total.quad.pan) # p < 0.001
lrtest(base.cub.pan, total.cub.pan) # p < 0.001

# nonpandemic
lrtest(base.lin.nonpan, total.lin.nonpan) # p < 0.001
lrtest(base.quad.nonpan, total.quad.nonpan) # p < 0.001
lrtest(base.cub.nonpan, total.cub.nonpan) # p < 0.001


#######################################
# MAKE FIGURE PUTTING IT ALL TOGETHER #
#######################################

# add color palette
colors <- c("#A3FAD0", "#43524E", "#FA815C", "#74A598", "#7A6E6A")

# PANDEMIC 

lin.pan.train$probs <- predict(total.lin.pan, lin.pan.train, type = "response")
linear <- lin.pan.train
linear$model <- "x"

quad.pan.train$probs <- predict(total.quad.pan, quad.pan.train, type = "response")
quadratic <- quad.pan.train
quadratic$model <- "x^2"

cub.pan.train$probs <- predict(total.cub.pan, cub.pan.train, type = "response")
cubic <- cub.pan.train
cubic$model <- "x^3"

agg <- rbind(linear, quadratic, cubic)
agg$model <- as.factor(agg$model)
agg$district <- ifelse(agg$district == "1", "Southeast", 
                       ifelse(agg$district == "2", "Western", 
                              ifelse(agg$district == "3", "South Central",
                                     ifelse(agg$district == "4", "Southwest",
                                            ifelse(agg$district == "5", "Interior", NA)))))
agg$race_cat <- ifelse(agg$race_cat == "indig", "Alaska Native", "Non-Alaska Native")
agg$district <- as.factor(agg$district)

tiff("PANDModelPredictions.tiff", height = 3.6, width = 7, units = "in", res = 600)
ggplot(agg, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(alpha = 0.15) +
  geom_smooth(data = subset(agg, model == "x"), method = "lm", size = 0.5, alpha = 0.175) +
  geom_smooth(data = subset(agg, model == "x^2"), method = "lm", formula = y~x+I(x^2), size = 0.5, alpha = 0.175) +
  geom_smooth(data = subset(agg, model == "x^3"), method = "lm", formula = y~x+I(x^2)+I(x^3), size = 0.5, alpha = 0.175) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  #scale_color_viridis_d(direction = -1) +
  #scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(model), rows = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1) +
  labs(
    x = "BMI",
    y = "Probability of P&I Death",
    color = "District",
    fill = "District"
  )
dev.off()


# NONPANDEMIC
lin.nonpan.train$probs <- predict(total.lin.nonpan, lin.nonpan.train, type = "response")
linear <- lin.nonpan.train
linear$model <- "x"

quad.nonpan.train$probs <- predict(total.quad.nonpan, quad.nonpan.train, type = "response")
quadratic <- quad.nonpan.train
quadratic$model <- "x^2"

cub.nonpan.train$probs <- predict(total.cub.nonpan, cub.nonpan.train, type = "response")
cubic <- cub.nonpan.train
cubic$model <- "x^3"

agg <- rbind(linear, quadratic, cubic)
agg$model <- as.factor(agg$model)
agg$district <- ifelse(agg$district == "1", "Southeast", 
                       ifelse(agg$district == "2", "Western", 
                              ifelse(agg$district == "3", "South Central",
                                     ifelse(agg$district == "4", "Southwest",
                                            ifelse(agg$district == "5", "Interior", NA)))))
agg$district <- as.factor(agg$district)
agg$race_cat <- ifelse(agg$race_cat == "indig", "Alaska Native", "Non-Alaska Native")

tiff("NONPANDModelPredictions.tiff", height = 3.6, width = 7, units = "in", res = 600)
ggplot(agg, aes(x = bmi, y = probs, color = district, fill = district)) +
  geom_point(alpha = 0.15) +
  geom_smooth(data = subset(agg, model == "x"), method = "lm", size = 0.5, alpha = 0.175) +
  geom_smooth(data = subset(agg, model == "x^2"), method = "lm", formula = y~x+I(x^2), size = 0.5, alpha = 0.175) +
  geom_smooth(data = subset(agg, model == "x^3"), method = "lm", formula = y~x+I(x^2)+I(x^3), size = 0.5, alpha = 0.175) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  #scale_color_viridis_d(direction = -1) +
  #scale_fill_viridis_d(direction = -1) +
  facet_grid(cols = vars(model), rows = vars(race_cat)) +
  theme_minimal() +
  ylim(0,1) +
  labs(
    x = "BMI",
    y = "Probability of P&I Death",
    color = "District",
    fill = "District",
  )
dev.off()


#################
# OTHER FIGURES #
#################

# boxplot of BMI of flu and other per year
bp.dat <- rbind(pand, nonpand)
bp.dat$logcod <- as.factor(bp.dat$logcod)
bp.dat$logcod <- ifelse(bp.dat$logcod == "1", "P&I", "Other")
tiff("Boxplots.tiff", height = 3.75, width = 10, units = "in", res = 600)
ggplot(bp.dat, aes(x = yod, y = bmi, col = cod, fill = cod)) +
  geom_boxplot(outlier.shape = 1, alpha = 0.7) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  #scale_color_viridis_d(option = "mako") +
  #scale_fill_viridis_d(option = "mako") +
  theme_minimal() +
  labs(
    x = "Year",
    y = "BMI",
    color = "Cause of Death",
    fill = "Cause of Death"
  )
dev.off()

# histograms of indigenous and non-indigenous BMIs stratified by cause of death
data$race_cat <- ifelse(data$race_cat == "indig", "Alaska Native", "Non-Alaska Native")
data <- data %>% filter(!is.na(bmi))

data2 <- rbind(pand, nonpand)
data2$race_cat <- ifelse(data2$race_cat == "indig", "Alaska Native", "Non-Alaska Native")

list(summary(subset(data2, data2$race_cat == "Alaska Native")[19]), summary(subset(data2, data2$race_cat == "Non-Alaska Native")[19]))
nat.bmi <- 23.3
non.bmi <- 25.2
bmis <- data.frame(as.factor(c("nat.bmi", "non.bmi")), c(nat.bmi, non.bmi), as.factor(c("μ(AK Nat)", "μ(Non-AK Nat)")))
colnames(bmis) <- c("group", "bmi", "line")

# tests for differences in mean between AK Native and non-AK Native subsets of the population
t.test(subset(data2, data2$race_cat == "Alaska Native")[19], subset(data2, data2$race_cat == "Non-Alaska Native")[19]) # p < 0.001***

tiff("BMIHistograms.tiff", height = 2.75, width = 10, units = "in", res = 600)
ggplot(data2, aes(x = bmi, y = ..density.., color = cod, fill = cod, group = cod)) +
  geom_density(alpha = 0.1) +
  facet_grid(col = vars(race_cat)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  geom_vline(data = bmis, aes(xintercept = bmi, linetype = line)) +
  annotate("text", x = 35.5, y = 0.14, label = "p < 0.001 ***") +
  ylim(0,0.18) +
  xlim(0,50) +
  xlab("BMI") +
  ylab("Density") +
  labs(
    fill = "Cause of Death",
    color = "Cause of Death"
  )
dev.off()

# calculate t.tests between flu and others per year
t.17 <- t.test(subset(data, data$yod == "1917" & data$logcod == 1)$bmi, subset(data, data$yod == "1917" & data$logcod == 0)$bmi)
t.18 <- t.test(subset(data, data$yod == "1918" & data$logcod == 1)$bmi, subset(data, data$yod == "1918" & data$logcod == 0)$bmi)
t.19 <- t.test(subset(data, data$yod == "1919" & data$logcod == 1)$bmi, subset(data, data$yod == "1919" & data$logcod == 0)$bmi)
t.20 <- t.test(subset(data, data$yod == "1920" & data$logcod == 1)$bmi, subset(data, data$yod == "1920" & data$logcod == 0)$bmi)
t.21 <- t.test(subset(data, data$yod == "1921" & data$logcod == 1)$bmi, subset(data, data$yod == "1921" & data$logcod == 0)$bmi)
t.22 <- t.test(subset(data, data$yod == "1922" & data$logcod == 1)$bmi, subset(data, data$yod == "1922" & data$logcod == 0)$bmi)
t.23 <- t.test(subset(data, data$yod == "1923" & data$logcod == 1)$bmi, subset(data, data$yod == "1923" & data$logcod == 0)$bmi)
t.24 <- t.test(subset(data, data$yod == "1924" & data$logcod == 1)$bmi, subset(data, data$yod == "1924" & data$logcod == 0)$bmi)
t.25 <- t.test(subset(data, data$yod == "1925" & data$logcod == 1)$bmi, subset(data, data$yod == "1925" & data$logcod == 0)$bmi)
p.val <- c(t.17$p.value, t.18$p.value, t.19$p.value, t.20$p.value, t.21$p.value, t.22$p.value, t.23$p.value, t.24$p.value, t.25$p.value)
flu.mean <- c(t.17$estimate, t.18$estimate, t.19$estimate, t.20$estimate, t.21$estimate, t.22$estimate, t.23$estimate, t.24$estimate, t.25$estimate)
other.mean <- c(t.17$estimate[2], t.18$estimate[2], t.19$estimate[2], t.20$estimate[2], t.21$estimate[2], t.22$estimate[2], t.23$estimate[2], t.24$estimate[2], t.25$estimate[2])
low.intervals <- c(t.17$conf.int, t.18$conf.int, t.19$conf.int, t.20$conf.int, t.21$conf.int, t.22$conf.int, t.23$conf.int, t.24$conf.int, t.25$conf.int)
high.intervals <- c(t.17$conf.int[2], t.18$conf.int[2], t.19$conf.int[2], t.20$conf.int[2], t.21$conf.int[2], t.22$conf.int[2], t.23$conf.int[2], t.24$conf.int[2], t.25$conf.int[2])
year <- 1917:1925
sig.test <- data.frame(year, flu.mean, other.mean, p.val, low.intervals, high.intervals)




# put together dataset
pand$time <- rep("pandemic", nrow(pand))
nonpand$time <- rep("non-pandemic", nrow(nonpand))
data <- rbind(pand, nonpand)
write_csv(data, "data.csv")
