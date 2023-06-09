setwd("/Users/Mel/Desktop/All Phyla Species")
library(leaps)
library(tidyverse)
library(caret)
library(MASS)

taxalevel = tail(strsplit(getwd(),split=" ")[[1]],1)

#Load data
topology = read.csv("network_topology.csv", row.names = 1)
topology <- subset(topology, Biomass.Quantile != "All Quantiles")

# remove 'inf' values
topology = topology[!grepl("Inf", topology$Pos.Neg),]

# convert Biomass.Quantile to ordered factor
topology$Biomass.Quantile <- ordered(topology$Biomass.Quantile, 
                                   levels = c("physeq_subpop_glom_1", "physeq_subpop_glom_2", "physeq_subpop_glom_3", "physeq_subpop_glom_4"))
############

# record network topology data
vars = as.vector(colnames(topology)[-c(1:3)])
#topsum = data.frame()
for(var in vars){
  min = summary(topology[[var]])[1]
  max = summary(topology[[var]])[6]
  mean = summary(topology[[var]])[3]
  sd = sd(topology[[var]])
  
  temp <- c(taxalevel, var, min, max, mean, sd)
  topsum = rbind(topsum, temp)
}
colnames(topsum) <- c("Taxonomic Level", "Network Topology Variable", "Minimum","Maximum","Mean","SD")
write.csv(topsum, "/Users/Mel/Desktop/topsum.csv")

############

#Remove population, taxa.level, multicollinear variables, and size-dependent variables
subtop = subset(topology, select=-c(Positive.Edges,Negative.Edges,Pos.Total,Neg.Total,Population,Taxa.Level, Avg.Degree, Total.Edges, Total.Nodes))

# ##removing total nodes counts because of perfect separation across quantiles at these taxonomic levels
# if(taxalevel %in% c("Genus","Species")){
#   subtop = subset(subtop, select=-c(Total.Nodes))
# }



##################### leaps Method of determining best fit: https://towardsdatascience.com/selecting-the-best-predictors-for-linear-regression-in-r-f385bf3d93e9

#Run the regsubsets() function on all variables.
Best_Subset <-
  regsubsets(Biomass.Quantile~.,
             data =subtop,
             nbest = 1,      # 1 best model for each number of predictors
             nvmax = NULL,    # NULL for no limit on number of variables
             force.in = NULL, force.out = NULL,
             method = "exhaustive")
summary_best_subset <- summary(Best_Subset)
as.data.frame(summary_best_subset$outmat)
which.max(summary_best_subset$adjr2) #see what the package recommends in terms of the number of predictors to use for our dataset
predictors = summary_best_subset$which[which.max(summary_best_subset$adjr2),] #What are the best predictors? best predictors are indicated by ‘TRUE’.
predictors = names(predictors)[predictors]
predictors = predictors[-1]
predictors_leaps = c(predictors, "Biomass.Quantile")

predictors_leaps
subtop_leaps = subtop[predictors_leaps]

##################### AIC Method of determining best fit: http://www.sthda.com/english/articles/36-classification-methods-essentials/150-stepwise-logistic-regression-essentials-in-r/

# Split the data into training and test set
#set.seed(123)
training.samples <- subtop$Biomass.Quantile %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- subtop[training.samples, ]
test.data <- subtop[-training.samples, ]

# # Fit the model
# model <- MASS::polr(formula = Biomass.Quantile ~., data = subtop) %>%
#   MASS::stepAIC(trace = FALSE)
# 
# # Summarize the final selected model
# summary(model)
# # Make predictions
# probabilities <- model %>% predict(test.data, type = "probs")
# predicted.classes <- colnames(probabilities)[apply(probabilities,1,which.max)]
# # Model accuracy
# observed.classes <- test.data$Biomass.Quantile
# mean(predicted.classes == observed.classes)

## Full logistic regression model. (Incorporating all predictors)
full.model <- MASS::polr(Biomass.Quantile ~., data = train.data)
#coef(full.model)

## Perform stepwise variable selection. (Select the most contributive variables)
step.model <- full.model %>% MASS::stepAIC(trace = FALSE)
#coef(step.model)
predictors_aic = as.character(names(coef(step.model)))
predictors_aic = c(predictors_aic, "Biomass.Quantile")

predictors_aic
subtop_aic = subtop[predictors_aic]

## leaps model
leaps.model <- MASS::polr(Biomass.Quantile ~., data = subtop_leaps)

## Compare the full, stepwise, and leaps models
# Here, we’ll compare the performance of the full and the stepwise logistic models. The best model is defined as the model that has the lowest classification error rate in predicting the class of new test data:

# Prediction accuracy of the full logistic regression model:
# Make predictions
probabilities <- full.model %>% predict(test.data, type = "probs")
predicted.classes <- colnames(probabilities)[apply(probabilities,1,which.max)]
# Model accuracy (Full)
observed.classes <- test.data$Biomass.Quantile
mean(predicted.classes == observed.classes)

# Prediction accuracy of the stepwise logistic regression model:
# Make predictions
probabilities <- predict(step.model, test.data, type = "probs")
predicted.classes <- colnames(probabilities)[apply(probabilities,1,which.max)]
# Model accuracy (Stepwise)
observed.classes <- test.data$Biomass.Quantile
mean(predicted.classes == observed.classes)

# Prediction accuracy of the leaps logistic regression model:
# Make predictions
probabilities <- leaps.model %>% predict(test.data, type = "probs")
predicted.classes <- colnames(probabilities)[apply(probabilities,1,which.max)]
# Model accuracy (Leaps)
observed.classes <- test.data$Biomass.Quantile
mean(predicted.classes == observed.classes)


#########################


# get summary
summary(step.model)
# get coefficients (it's in matrix form)
coefficients <- summary(step.model)$coefficients
# calculate p-values
p_value <- (1 - pnorm(abs(coefficients[ ,"t value"]), 0, 1))*2
# bind back to coefficients
coefficients <- cbind(coefficients, p_value)

# calculate odds ratios
odds_ratio <- exp(coefficients[ ,"Value"])
# combine with coefficient and p_value
(coefficients <- cbind(
  coefficients[ ,c("Value", "p_value")],
  odds_ratio
))

#Calculate percent lower/higher odds 
(percent_odds <- ifelse(p_value < 0.05, paste0(as.character(round((coefficients[,"odds_ratio"]-1)*100, 2)), "%"), "N/A (p > 0.05)"))

# ##########################
# # record network topology data
# vars = as.vector(predictors_aic[-length(predictors_aic)])
# #supptable = data.frame()
# for(var in vars){
#   min = summary(subtop[[var]])[1]
#   max = summary(subtop[[var]])[6]
#   mean = summary(subtop[[var]])[3]
#   sd = sd(subtop[[var]])
#   
#   temp <- c(taxalevel, var, min, max, mean, sd)
#   supptable = rbind(supptable, temp)
# }
# colnames(supptable) <- c("Taxonomic Level", "Network Topology Variable", "Minimum","Maximum","Mean","SD")
# write.csv(supptable, "/Users/Mel/Desktop/supptable.csv")
# #########################

