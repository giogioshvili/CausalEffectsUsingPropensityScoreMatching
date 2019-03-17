#Load packages
library(corrplot) 
library(data.table)
library(ggplot2)
library(dplyr)
library(MatchIt)
library(stats)
library(ggplot2)
library(rbounds)
library(cobalt)


#Load the data
data <- fread('criteo-uplift.csv') #This step takes around 30 sec

#---------------------------------Analysis of RCT-----------------------------------

#Review & summarize data 
str(data)
summary(data)
data[, .N, treatment]
data[, .N, exposure]
data[, .N, conversion]



#Calculate the correlation matrix and draw plot
cormat=cor(data)
corrplot(cormat, method = "circle",tl.col = "black") #plot matrix


#Simple T test for Causal Effects 
with(data, t.test(conversion ~ treatment))

#calculate mean conversion rates for all 4 groups as a percentage
a <- data[treatment==1, mean(conversion)]*100
b <- data[treatment==0, mean(conversion)]*100
c <- data[exposure==1, mean(conversion)]*100     #conversion rate for exposed users
d <- data[exposure==0, mean(conversion)]*100     #conversion rate for unexposed users
e <- data[exposure==0 & treatment==0, mean(conversion)]*100
#Draw 1x4 matrix to represent average conversions among all groups 
matrix <- c(a,b,c,d,e)
average_conversion <- matrix(data=matrix, ncol=5, nrow=1, byrow=TRUE)
colnames(average_conversion) <- c("Average Conversion when Treatment is 1", "Average Conversion when Treatment is 0", "Average Conversion when Exposure is 1", "Average Conversion when Exposure is 0", "Average conversion when exposure & treatment is 0")

#Intend to treat (ITT) & lift as percentage 
ITT <- (a-b)
ITT.lift <- (a-b)/b*100

#Calculate exposure as percentage
exposure_rate <- mean(data$exposure)*100                     


#Average treatment on treated (ATT) & lift as percentage
ATT <- ITT/(exposure_rate/100)
ATT.lift <- ATT/(c-ATT)*100

#####Re-structure dataset to make RCT comparable to observational method and computationally inexpensive

#Set a random seed & choose sample size
set.seed(314)
sample_size = 50000

#Create sub-sample
#Since this size of the sample (50K) mimics the whole dataset and also gives similar intuition as well as similar causal inference, there is no need to increase the sample size
sub_sample <- data[sample(nrow(data), sample_size),]
summary(sub_sample)

#Simple test if there is a significant difference between sample and population (full dataset) exposures,
(mean(sub_sample$exposure) - mean(data$exposure))/sqrt(sd(sub_sample$exposure)*sd(data$exposure))
(mean(sub_sample$conversion) - mean(data$conversion))/sqrt(sd(sub_sample$conversion)*sd(data$conversion))
#=> no significant difference between distribution of sample and distribution of population



###Repeat some of the above steps to obtain causal effects on a new sub-sampled data

with(sub_sample, t.test(conversion ~ exposure))

a <- sub_sample[treatment==1, mean(conversion)]*100
b <- sub_sample[treatment==0, mean(conversion)]*100
c <- sub_sample[exposure==1, mean(conversion)]*100
d <- sub_sample[exposure==0, mean(conversion)]*100
e <- sub_sample[exposure==0 & treatment==0, mean(conversion)]*100

matrix <- c(a,b,c,d,e)
average_conversion <- matrix(data=matrix, ncol=5, nrow=1, byrow=TRUE)
colnames(average_conversion) <- c("Average Conversion when Treatment is 1", "Average Conversion when Treatment is 0", "Average Conversion when Exposure is 1", "Average Conversion when Exposure is 0","Average conversion when exposure & treatment is 0")

ITT <- (a-b)
ITT.lift <- (a-b)/b*100

exposure_rate <- mean(sub_sample$exposure)*100                     


ATT <- ITT/(exposure_rate/100)
ATT.lift <- ATT/(c-ATT)*100


#Causal effects are also represented by coefficient OLS of simple/full model.
#Simple model (exposure & conversion)
casual_effect_linear <- lm(conversion ~ exposure, sub_sample)
summary(casual_effect_linear)

#full model (between conversion & exposure + all covariates)
full_formula <-as.formula('conversion ~ exposure + f0 + f1 + f2 +f3 + f4 + f5 +f6 +f7 + f8 + f9 + f10 + f11')
casual_effect_linear_full <- lm(full_formula, sub_sample)
summary(casual_effect_linear_full)

###-----------------------------Observational Method Analysis-----------------------
##------------------------------Propensity Score Matching---------------------------



#PS matching model: Applied on the same data sample as the RCT("data1")
##Propensity score formula definition
PS_formula <- as.formula('exposure ~ f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11')

#Estimate PS scores and run the matching model 
#Add time tracker to check how much time does the algoritm need on the particular sample size
st <- Sys.time()
match_model = matchit(PS_formula,
                      data = sub_sample, 
                      method = "nearest",  
                      distance = "logit")
Sys.time() - st

summary(match_model)

##Match the data 
PS_matched <- match.data(match_model)
# convert into data.table
PS_matched <- as.data.table(PS_matched)
##Review the data 
dim(PS_matched)
summary(PS_matched)

#Visualize the performance of the matching model 
#QQplot
plot(match_model)
plot(match_model, type = 'jitter', interactive = FALSE)

##Review distribution feature by feature distribution.
bal.plot(match_model, 'f2')
bal.plot(match_model, var.name = "f1", which = "both",
         type = "histogram", mirror = TRUE)

bal.plot(match_model, "distance", which = "both")

##Review all features difference in adjusteed and unadjusted samples at once using love chart 
love.plot(bal.tab(match_model), threshold = .1)


#Obtaining Causal Effects (ATT)
##Simple T test
with(PS_matched, t.test(conversion ~ exposure))

f <- PS_matched[exposure==1, mean(conversion)]*100   #conversion rate for exposed users (it should be same as the one in initial sample (data1))
g <- PS_matched[exposure==0, mean(conversion)]*100   #conversion rate for unexposed users

#Calculate exposure and conversion as percentage
exposure_rate_on_matched <- mean(PS_matched$exposure)*100                     

ATT_PS <- f-g
ATT.lift_PS <- ATT_PS/(f-ATT_PS)*100


#Linear model 
#Causal effect can also be defined by coefficient OLS of simple/full model.
#Simple model (exposure & conversion)

casual_effect_linear_matched<- lm(conversion ~ exposure, PS_matched)
summary(casual_effect_linear_matched)

#full model 
casual_effect_linear_full_matched <- lm(full_formula, PS_matched)
summary(casual_effect_linear_full_matched)

#Testing joint sigfinicance of f0-f11 ( F test )
anova(casual_effect_linear_full_matched, casual_effect_linear_matched)
#> the f0, f11 are joinly significant


###Fit the model on pre-matched data. 
  #sub_sample$fitted_conversion <- (predict(casual_effect_linear_full, sub_sample))

###Fitting linear model for binary requires assuming binary probability bellongs to [0,1]
  #sub_sample$fitted_conversion <- ifelse(sub_sample$fitted_conversion <0,0,sub_sample$fitted_conversion)


###Distribution of fitted values 
  #qplot() + geom_density(aes(sub_sample$fitted_conversion), fill = 'blue')




