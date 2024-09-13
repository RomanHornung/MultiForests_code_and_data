# This function performs one repetition of the five times repeated 
# 5-fold stratified cross-validation on a specific data set using
# one of the five compared methods.
#
# It takes the whole number 'iter', which corresponds to the iter-th line 
# of 'scenariogrid', which contains the necessary information
# on the iter-th setting.

# Input parameters:

# iter  - iter-th line of 'scenariogrid', which corresponds to the iter-th line 
#         of 'scenariogrid', which contains the necessary information
#         on the iter-th setting.
# mywd  - working directory. This needs to be the path to 
#         the directory "MultiForests_code_and_data/benchmark_study".

# Output:

# A list of length 5, where each element contains the values of the four performance metrics
# for one of the 5 cross-validation iterations

evaluatesetting <- function(iter, mywd) {
  
  setwd(mywd)
  
  # Obtain information for the iter-th setting:
  
  dataset <- scenariogrid$dataset[iter]
  method <- scenariogrid$method[iter]
  cvind <- scenariogrid$cvind[iter]
  seed <- scenariogrid$seed[iter]
  settingid <- scenariogrid$settingid[iter]
  
  
  # Set seed:
  
  set.seed(seed)
  
  
  # Load data set:
  
  load(paste0("./data/datasets/", dataset))
  
  if (nrow(dataset) <= 5000)
    num_trees <- 5000
  else
    num_trees <- 1000
  
  
  # Make CV division:
  
  Kcv <- 5
  cvdiv <- makeCVdiv(data=dataset, yname="ytarget", ncv=Kcv)
  
  
  # Perform cross-validation: 
  
  res <- list()
  
  for(i in 1:Kcv) {
    
    datatrain <- dataset[cvdiv!=i,]
    datatest <- dataset[cvdiv==i,]
    
    if(method=="rf") {
      
      res[[i]] <- trainandpredict(datatrain, datatest, trainfun=trainRandomForest, 
                                  predictfun=predictRandomForest, num_trees=num_trees)
      
    }
    
    if(method=="muwf_wgini_wsquared") {
      
      res[[i]] <- trainandpredict(datatrain, datatest, trainfun=trainMuwfWginiWsquared, 
                                  predictfun=predictMuwfWginiWsquared, num_trees=num_trees)
      
    }
    
    if(method=="muwf_wogini_wsquared") {
      
      res[[i]] <- trainandpredict(datatrain, datatest, trainfun=trainMuwfWoginiWsquared, 
                                  predictfun=predictMuwfWoginiWsquared, num_trees=num_trees)
      
    }
    
    if(method=="muwf_wgini_wosquared") {
      
      res[[i]] <- trainandpredict(datatrain, datatest, trainfun=trainMuwfWginiWosquared, 
                                  predictfun=predictMuwfWginiWosquared, num_trees=num_trees)
      
    }
    
    if(method=="muwf_wogini_wosquared") {
      
      res[[i]] <- trainandpredict(datatrain, datatest, trainfun=trainMuwfWoginiWosquared, 
                                  predictfun=predictMuwfWoginiWosquared, num_trees=num_trees)
      
    }
    
  }
  
  
  # Return the results:
  
  return(res)
  
}







# This function performs one iteration of the 5-fold stratified cross-validation:

# Input parameters:

# datatrain      - training data set
# datatest       - test data set
# trainfun       - function used for training the classifier
# predictfun     - function used for predicting with the classifier
# num_trees      - the number of trees to construct for each forest

# Output:

# The values of the calculated performance metrics

trainandpredict <- function(datatrain, datatest, trainfun, predictfun, num_trees) {
  
  # Construct forest:
  model <- trainfun(dependent.variable.name="ytarget", data=datatrain, num_trees=num_trees)
  
  preds <- predictfun(model, dependent.variable.name="ytarget", newdata=datatest)
  probpreds <- preds$probpreds
  classpreds <- preds$classpreds
  
  # Calculate ACC values:
  acctest <- measureACC(truth=datatest$ytarget, response=classpreds)
  
  # Calculate Brier values:
  briertest <- measureBrier(probabilities=probpreds, truth=datatest$ytarget)
  
  aunptest <- measureAUNP(probabilities=probpreds, truth=datatest$ytarget)
  
  aunutest <- measureAUNU(probabilities=probpreds, truth=datatest$ytarget)
  
  # Return results:
  return(list(acctest=acctest, briertest=briertest, aunptest=aunptest, aunutest=aunutest))
  
}








# Function used for training the random forests:

# Input parameters:

# dependent.variable.name - name of the outcome variable
# data                    - training data set
# num_trees               - the number of trees to construct for the forest

# Output:

# An object containing the trained forest.

trainRandomForest <- function(dependent.variable.name, data, num_trees) {
  
  require("ranger")
  
  model <- ranger::ranger(dependent.variable.name = "ytarget", num.trees=num_trees, data=data, 
                          importance="none", respect.unordered.factors="order", replace = FALSE,
                          sample.fraction = 0.7, probability=TRUE, min.node.size=5, num.threads=1)
  
  return(model)
  
}


# Function used for predicting with the trained random forests:

# Input parameters:

# model                   - trained forest returned by "trainRandomForest"
# dependent.variable.name - name of the outcome variable
# newdata                 - test data for which predictions should be obtained.

# Output:

# A list with two elements:

# probpreds  - matrix with class probability predictions, where each row contains
#              the class probabilities for one observation
# classpreds - a factor vector that contains the predicted classes

predictRandomForest <- function(model, dependent.variable.name, newdata) {
  
  require("ranger")
  
  # Predict conditional class probabilities for the observations in
  # the test data set:
  probpreds <- ranger:::predict.ranger(model, data=newdata)$predictions
  
  # In rare cases, some outcome classes were present only in the test data 
  # and not in the training data. In these situations, the predict function 
  # above did not include these classes in the matrix of probability predictions. 
  # The following code adds zero probabilities for these classes:
  if (!all(levels(newdata$ytarget) %in% colnames(probpreds))) {
    probpredsnew <- matrix(nrow=nrow(probpreds), ncol=length(levels(newdata$ytarget)), data=0)
    indsfill <- which(levels(newdata$ytarget) %in% colnames(probpreds))
    reorderind <- as.numeric(factor(levels(newdata$ytarget)[levels(newdata$ytarget) %in% colnames(probpreds)], 
                                    levels=colnames(probpreds)))
    probpreds <- probpreds[,reorderind]
    probpredsnew[,indsfill] <- probpreds
    probpreds <- probpredsnew
    colnames(probpreds) <- levels(newdata$ytarget)
  }
  
  # Obtain class predictions:
  classpreds <- factor(levels(newdata[,dependent.variable.name])[apply(probpreds, 1, nnet::which.is.max)], 
                       levels=levels(newdata[,dependent.variable.name]))
  
  return(list(probpreds=probpreds, classpreds=classpreds))
  
}




# Function used for training the version of multi forests with variants "Gini"
# and "Squared":

# Input parameters:

# dependent.variable.name - name of the outcome variable
# data                    - training data set
# num_trees               - the number of trees to construct for the forest

# Output:

# An object containing the trained forest.

trainMuwfWginiWsquared <- function(dependent.variable.name, data, num_trees) {
  
  require("diversityForestA")
  
  model <- diversityForestA::multifor(dependent.variable.name = "ytarget", data=data, 
                                      importance="none", npervar = 5, num.trees=num_trees, 
                                      replace = FALSE, sample.fraction = 0.7, probability=TRUE,
                                      min.node.size=5, num.threads=1)
  
  return(model)
  
}



# Function used for predicting with the version of multi forests with variants 
# "Gini" and "Squared":

# Input parameters:

# model                   - trained forest returned by "trainMuwfWginiWsquared"
# dependent.variable.name - name of the outcome variable
# newdata                 - test data for which predictions should be obtained.

# Output:

# A list with two elements:

# probpreds  - matrix with class probability predictions, where each row contains
#              the class probabilities for one observation
# classpreds - a factor vector that contains the predicted classes

predictMuwfWginiWsquared <- function(model, dependent.variable.name, newdata) {
  
  require("diversityForestA")
  
  # Predict conditional class probabilities for the observations in
  # the test data set:
  probpreds <- diversityForestA:::predict.multifor(model, data=newdata)$predictions
  
  if (!all(levels(newdata$ytarget) %in% colnames(probpreds))) {
    probpredsnew <- matrix(nrow=nrow(probpreds), ncol=length(levels(newdata$ytarget)), data=0)
    indsfill <- which(levels(newdata$ytarget) %in% colnames(probpreds))
    reorderind <- as.numeric(factor(levels(newdata$ytarget)[levels(newdata$ytarget) %in% colnames(probpreds)], 
                                    levels=colnames(probpreds)))
    probpreds <- probpreds[,reorderind]
    probpredsnew[,indsfill] <- probpreds
    probpreds <- probpredsnew
    colnames(probpreds) <- levels(newdata$ytarget)
  }
  
  # Obtain class predictions:
  classpreds <- factor(levels(newdata[,dependent.variable.name])[apply(probpreds, 1, nnet::which.is.max)], 
                       levels=levels(newdata[,dependent.variable.name]))
  
  return(list(probpreds=probpreds, classpreds=classpreds))
  
}





# Function used for training the version of multi forests with variants "Assign Classes"
# and "Squared":

# Input parameters:

# dependent.variable.name - name of the outcome variable
# data                    - training data set
# num_trees               - the number of trees to construct for the forest

# Output:

# An object containing the trained forest.

trainMuwfWoginiWsquared <- function(dependent.variable.name, data, num_trees) {
  
  require("diversityForestB")
  
  model <- diversityForestB::multifor(dependent.variable.name = "ytarget", data=data, 
                                      importance="none", npervar = 5, num.trees=num_trees, 
                                      replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
                                      min.node.size=5, num.threads=1)
  
  return(model)
  
}


# Function used for predicting with the version of multi forests with variants 
# "Assign Classes" and "Squared":

# Input parameters:

# model                   - trained forest returned by "trainMuwfWoginiWsquared"
# dependent.variable.name - name of the outcome variable
# newdata                 - test data for which predictions should be obtained.

# Output:

# A list with two elements:

# probpreds  - matrix with class probability predictions, where each row contains
#              the class probabilities for one observation
# classpreds - a factor vector that contains the predicted classes

predictMuwfWoginiWsquared <- function(model, dependent.variable.name, newdata) {
  
  require("diversityForestB")
  
  # Predict conditional class probabilities for the observations in
  # the test data set:
  probpreds <- diversityForestB:::predict.multifor(model, data=newdata)$predictions
  
  if (!all(levels(newdata$ytarget) %in% colnames(probpreds))) {
    probpredsnew <- matrix(nrow=nrow(probpreds), ncol=length(levels(newdata$ytarget)), data=0)
    indsfill <- which(levels(newdata$ytarget) %in% colnames(probpreds))
    reorderind <- as.numeric(factor(levels(newdata$ytarget)[levels(newdata$ytarget) %in% colnames(probpreds)], 
                                    levels=colnames(probpreds)))
    probpreds <- probpreds[,reorderind]
    probpredsnew[,indsfill] <- probpreds
    probpreds <- probpredsnew
    colnames(probpreds) <- levels(newdata$ytarget)
  }
  
  # Obtain class predictions:
  classpreds <- factor(levels(newdata[,dependent.variable.name])[apply(probpreds, 1, nnet::which.is.max)], 
                       levels=levels(newdata[,dependent.variable.name]))
  
  return(list(probpreds=probpreds, classpreds=classpreds))
  
}





# Function used for training the version of multi forests with variants "Gini"
# and "Non-Squared":

# Input parameters:

# dependent.variable.name - name of the outcome variable
# data                    - training data set
# num_trees               - the number of trees to construct for the forest

# Output:

# An object containing the trained forest.

trainMuwfWginiWosquared <- function(dependent.variable.name, data, num_trees) {
  
  require("diversityForestC")
  
  model <- diversityForestC::multifor(dependent.variable.name = "ytarget", data=data, 
                                      importance="none", npervar = 5, num.trees=num_trees, 
                                      replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
                                      min.node.size=5, num.threads=1)
  
  return(model)
  
}


# Function used for predicting with the version of multi forests with variants 
# "Gini" and "Non-Squared":

# Input parameters:

# model                   - trained forest returned by "trainMuwfWginiWosquared"
# dependent.variable.name - name of the outcome variable
# newdata                 - test data for which predictions should be obtained.

# Output:

# A list with two elements:

# probpreds  - matrix with class probability predictions, where each row contains
#              the class probabilities for one observation
# classpreds - a factor vector that contains the predicted classes

predictMuwfWginiWosquared <- function(model, dependent.variable.name, newdata) {
  
  require("diversityForestC")
  
  # Predict conditional class probabilities for the observations in
  # the test data set:
  probpreds <- diversityForestC:::predict.multifor(model, data=newdata)$predictions
  
  if (!all(levels(newdata$ytarget) %in% colnames(probpreds))) {
    probpredsnew <- matrix(nrow=nrow(probpreds), ncol=length(levels(newdata$ytarget)), data=0)
    indsfill <- which(levels(newdata$ytarget) %in% colnames(probpreds))
    reorderind <- as.numeric(factor(levels(newdata$ytarget)[levels(newdata$ytarget) %in% colnames(probpreds)], 
                                    levels=colnames(probpreds)))
    probpreds <- probpreds[,reorderind]
    probpredsnew[,indsfill] <- probpreds
    probpreds <- probpredsnew
    colnames(probpreds) <- levels(newdata$ytarget)
  }
  
  # Obtain class predictions:
  classpreds <- factor(levels(newdata[,dependent.variable.name])[apply(probpreds, 1, nnet::which.is.max)], 
                       levels=levels(newdata[,dependent.variable.name]))
  
  return(list(probpreds=probpreds, classpreds=classpreds))
  
}





# Function used for training the version of multi forests with variants "Assign Classes"
# and "Non-Squared":

# Input parameters:

# dependent.variable.name - name of the outcome variable
# data                    - training data set
# num_trees               - the number of trees to construct for the forest

# Output:

# An object containing the trained forest.

trainMuwfWoginiWosquared <- function(dependent.variable.name, data, num_trees) {
  
  require("diversityForestD")
  
  model <- diversityForestD::multifor(dependent.variable.name = "ytarget", data=data, 
                                      importance="none", npervar = 5, num.trees=num_trees, 
                                      replace = FALSE, sample.fraction = 0.7, probability=TRUE,
                                      min.node.size=5, num.threads=1)
  
  return(model)
  
}


# Function used for predicting with the version of multi forests with variants 
# "Assign Classes" and "Non-Squared":

# Input parameters:

# model                   - trained forest returned by "trainMuwfWoginiWosquared"
# dependent.variable.name - name of the outcome variable
# newdata                 - test data for which predictions should be obtained.

# Output:

# A list with two elements:

# probpreds  - matrix with class probability predictions, where each row contains
#              the class probabilities for one observation
# classpreds - a factor vector that contains the predicted classes

predictMuwfWoginiWosquared <- function(model, dependent.variable.name, newdata) {
  
  require("diversityForestD")
  
  # Predict conditional class probabilities for the observations in
  # the test data set:
  probpreds <- diversityForestD:::predict.multifor(model, data=newdata)$predictions
  
  if (!all(levels(newdata$ytarget) %in% colnames(probpreds))) {
    probpredsnew <- matrix(nrow=nrow(probpreds), ncol=length(levels(newdata$ytarget)), data=0)
    indsfill <- which(levels(newdata$ytarget) %in% colnames(probpreds))
    reorderind <- as.numeric(factor(levels(newdata$ytarget)[levels(newdata$ytarget) %in% colnames(probpreds)], 
                                    levels=colnames(probpreds)))
    probpreds <- probpreds[,reorderind]
    probpredsnew[,indsfill] <- probpreds
    probpreds <- probpredsnew
    colnames(probpreds) <- levels(newdata$ytarget)
  }
  
  # Obtain class predictions:
  classpreds <- factor(levels(newdata[,dependent.variable.name])[apply(probpreds, 1, nnet::which.is.max)], 
                       levels=levels(newdata[,dependent.variable.name]))
  
  return(list(probpreds=probpreds, classpreds=classpreds))
  
}






# Function that generates the splittings for stratified cross-validation:

# Input parameters:

# data   - data set
# yname  - label of the target variable
# ncv    - number of folds to use

# Output:

# A vector of length nrow(data), where the ith element contains the index
# of the fold in the K-fold cross-validation to which the ith observation is assigned.

makeCVdiv <- function(data, yname="y", ncv=5) {
  
  cvid <- rep(NA, length = nrow(data))
  for(i in levels(data[, yname])) cvid[data[, yname] == i] <- 
      sample(rep(1:ncv, length = sum(data[, yname] == i)))
  
  cvid
  
}








# Functions for calculating the performance metrics:


# Function calculating ACC values (taken from "mlr" R package):

measureACC <- function (truth, response) 
{
  mean(response == truth)
}


# Function calculating multi-class Brier values (taken from "mlr" R package):

measureBrier <- function (probabilities, truth) 
{
  class_mat <- model.matrix(~ 0 + truth)
  mean(rowSums((class_mat - probabilities)^2))
}


# Function calculating AUNU (taken from "mlr" R package):

measureAUNU <- function (probabilities, truth)
{
  mean(onevrestauc(probabilities, truth))
}


# Function calculating AUNP (taken from "mlr" R package):

measureAUNP <- function (probabilities, truth)
{
  #### CHANGED, ORIGINAL: sum(vapply(levels(truth), function(lvl) mean(truth == lvl), FUN.VALUE = NA_real_) * onevrestauc(prob, truth))
  sum(vapply(unique(as.character(truth)), function(lvl) mean(truth == lvl), FUN.VALUE = NA_real_) * onevrestauc(probabilities, truth))
}


# Function needed by the functinos calculating AUNU and AUNP (taken from "mlr" R package):

# returns a numeric length nlevel(truth), with one-vs-rest AUC
onevrestauc = function(prob, truth) {
  ntotal = nrow(prob)
  vapply(unique(as.character(truth)), function(cls) { #### CHANGED, ORIGINAL: vapply(levels(truth), function(cls) {
    nrest = sum(truth != cls)
    if (nrest == ntotal) {
      # this class has no members. What happened?
      #  - If we were called by mauc_aunu --> this does never happen, because we return(NA) if there are empty classes
      #  - If we were called by mauc_aunp --> this value gets multiplied with 0 and should result in 0
      # therefore we don't want to return Inf here, but a final value that does not matter in the end
      return(0)
    }
    
    r = rank(c(prob[truth == cls, cls], prob[truth != cls, cls]), ties.method = "average")
    # simplify the following:
    # (sum(r[seq_len(nrest)]) - nrest * (nrest + 1) / 2) / (nrest * (ntotal - nrest))
    (mean(r[seq_len(ntotal - nrest)]) - (ntotal - nrest + 1) / 2) / nrest
  }, FUN.VALUE = NA_real_)
}
