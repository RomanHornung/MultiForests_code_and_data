####################################################################################

# NOTE: Before the code can be executed, the R working directory *MUST* 
# be set to the directory 'MultiForests_code_and_data/benchmark_study/data' (that is
# the directory in which this R script is contained):

# Remove '#' from the line below and replace 'here/is/my/path/' by the path
# to the directory 'MultiForests_code_and_data/benchmark_study/data':

# setwd("here/is/my/path/")

####################################################################################

# Consider the subset of the OpenML-CC18 datasets that feature a maximum of
# 15000 observations and a maximum of 500 covariates:

library("OpenML")

datainfo <- read.table("./OpenML-CC18_table.txt", header=TRUE)
datainfosub <- datainfo[datainfo$n <= 15000 & datainfo$p <= 500 & datainfo$cl > 2,]




# Query multi-class data sets on OpenML that have a maximum of 20 classes, a maximum
# of 15000 observations, a maximum of 500 features and no missing values:

# tasks <- listOMLTasks(limit = NULL)
# classifTasks.infos <- subset(tasks, task.type == "Supervised Classification" &    # classification
#                               number.of.classes > 2 &                             # multi-class classification
#                                number.of.classes <= 20 &                             # not too many classes
#                                number.of.instances <= 15000 &   # remove too large datasets
#                               number.of.features <= 500 &  # remove datasets with too many features
#                               number.of.instances.with.missing.values == 0)           # no missing values

# save(classifTasks.infos, file = "./classifTasks.infos.RData")




# Exlude data sets in the same way as in the paper "Random forest versus logistic regression: a large-scale
# benchmark experiment" (code taken from that paper):

# We work with the given data from OpenML
load("./classifTasks.infos.RData")
datasets.index = sort(unique(classifTasks.infos$data.id))
clas = classifTasks.infos


# =============================
# Part 2 : Select datasets using OpenML task features
# ============================= 

print("2. Remove datasets using the tasks features", quote = FALSE)

# remove the redundancies : 470 datasets
clas = clas[order(clas$data.id),]
logic = diff(clas$data.id)>0
clas = clas[logic,]
print(paste("  Number of datasets after removing the redundacies of datasets's IDs :", dim(clas)[1]), quote = FALSE)

# Friedman-, volcanoes- und trX-Datasets : 393 datasets
clas = clas[substr(clas$name,1,9) != "volcanoes" & substr(clas$name,1,4) != "fri_" & 
              substr(clas$name,1,3) != "tr1" & substr(clas$name,1,3) != "tr2" & substr(clas$name,1,3) != "tr3" & 
              substr(clas$name,1,3) != "tr4", ]
print(paste("  Number of datasets after removing the obviously simulated datasets :", dim(clas)[1]), quote = FALSE)
print("  Datasets removed : Friedman, volcanoes, TrX-Datasets", quote = FALSE)

# remove the datasets with the same name, they correspond often to datasets with only very slight changes : 380 datasets
doublon = names(sort(table(clas$name)[table(clas$name) > 1]))
doublon = clas[clas$name %in% doublon,]
doublon = doublon[order(doublon$name), ]

diff.categorical <- function(x) {
  x = as.factor(x)
  n = length(x)
  res = rep(NA,n)
  res[1] = TRUE
  
  for (i in c(2:n)) {
    res[i] = !identical(x[i-1],x[i])
  }
  res = res*1
  return(res)
}

indexDuplicate.useful = which(diff.categorical(doublon$name)==1)
indexDuplicate.notuseful = which(diff.categorical(doublon$name)==0)
task.id.notuseful = doublon$task.id[indexDuplicate.notuseful]
indexclas.notuseful = which(clas$task.id %in% task.id.notuseful)
clas = clas[-indexclas.notuseful,]
print(paste("  Number of datasets after removing the redundancies in dataset's names:", 
            dim(clas)[1]), quote = FALSE)

# Ordering according to size (n*p)
clas = clas[order(clas$number.of.features * clas$number.of.instances), ]






# Remove the data sets that are already in OpenML-CC18:

clas <- clas[!(clas$data.id %in% datainfosub$Data_id),]
clas <- clas[!(clas$name %in% datainfosub$Name),]






# Combine the meta information from the OpenML-CC18 datasets and the remaining OpenML datasets:

datainfosub2 <- data.frame(
  Data_id = clas$data.id,
  Task_id = clas$task.id,
  Name = clas$name,
  cl = clas$number.of.classes,
  p = clas$number.of.features,
  n = clas$number.of.instances,
  stringsAsFactors = FALSE  # To ensure character columns are not converted to factors
)

part1 <- datainfosub
part1$MinMaj <- NULL

part2 <- datainfosub2

datatable <- rbind(part1, part2)




# There are 49 "timing-attack" datasets. Exclude all of them except one.
# I decided for the one that appeared first in a Google search of "timing-attack-dataset openml":

datatable <- datatable[-setdiff(grep("timing-attack-dataset", datatable$Name), 
                                grep("timing-attack-dataset-15-micro-seconds-delay-2022-09-18", 
                                     datatable$Name)),]





# Download the datasets:

for(i in 1:nrow(datatable)) {
  dataobj = getOMLDataSet(data.id = datatable$Data_id[i])
  dataset <- dataobj$data
  
  if("ytarget" %in% names(dataset))
    stop("Name 'ytarget' already taken")
  
  names(dataset)[names(dataset)==dataobj$target.features] <- "ytarget"
  
  if(class(dataset$ytarget)!="factor")
    dataset$ytarget <- factor(dataset$ytarget)
  
  save(dataset, file=paste0("./datasets/", datatable$Name[i], "_dataid_", datatable$Data_id[i], ".Rda"))
  
  cat(paste0("Iteration: ", i, " of ", nrow(datatable)), "\n")
}




# Query datasets from the Penn Machine Learning Benchmarks (PMLB) collection:

library("pmlbr")

library("dplyr")

datasetsrel <- summary_stats %>%
  filter(n_instances <= 15000 & n_features <= 500 & n_classes > 2 & n_classes <= 20, task == "classification") %>%
  pull(dataset)



# Exclude the datasets the names of which already occurred among the OpenML datasets:

datasetsrel <- datasetsrel[!(datasetsrel %in% datatable$Name)]



# Look for datasets that have similar names as OpenML datasets:

library("stringdist")

distances <- stringdistmatrix(datasetsrel, datatable$Name, method = "lv")
threshold <- 2
similar_indices <- which(distances <= threshold, arr.ind = TRUE)

similar_strings <- apply(similar_indices, 1, function(idx) {
  list(string_in_A = datasetsrel[idx[1]], string_in_B = datatable$Name[idx[2]])
})
similar_strings


# --> Exclude further datasets:

datasetsrel <- setdiff(datasetsrel, c("balance_scale", "mfeat_factors", "mfeat_fourier", 
                                      "mfeat_karhunen", "mfeat_morphological", "mfeat_zernike",
                                      "mfeat_pixel", "hayes_roth", "cleveland_nominal",
                                      "wine_quality_red", "car_evaluation", "wine_quality_white",
                                      "page_blocks"))


# Look for further datasets with less similar names:

distances <- stringdistmatrix(datasetsrel, datatable$Name, method = "lv")
threshold <- 3
similar_indices <- which(distances <= threshold, arr.ind = TRUE)

similar_strings <- apply(similar_indices, 1, function(idx) {
  list(string_in_A = datasetsrel[idx[1]], string_in_B = datatable$Name[idx[2]])
})
similar_strings


# Download the remaining datasets:

for(i in seq(along=datasetsrel)) {
  
  dataset <- fetch_data(datasetsrel[i])
  
  if("ytarget" %in% names(dataset))
    stop("Name 'ytarget' already taken")
  
  names(dataset)[names(dataset)=="target"] <- "ytarget"
  
  if(class(dataset$ytarget)!="factor")
    dataset$ytarget <- factor(dataset$ytarget)
  
  save(dataset, file=paste0("./datasets/", datasetsrel[i], "_dataid_NA.Rda"))
  
  cat(paste0("Iteration: ", i, " of ", length(datasetsrel)), "\n")
  
}



# Study remaining datasets that still have similar names as OpenML datasets:

# Each time compare the data set from PMLB with that from OpenML:

load("./datasets/cars_dataid_NA.Rda")
data1 <- dataset

load("./datasets/car_dataid_40975.Rda")
data2 <- dataset

dim(data1)
dim(data2)

# --> These are different datasets.



load("./datasets/cars_dataid_NA.Rda")
data1 <- dataset

load("./datasets/cars1_dataid_40700.Rda")
data2 <- dataset

dim(data1)
dim(data2)

head(data1)
head(data2)

# --> Basically the same dataset, remove the first one ("cars_dataid_NA.Rda").

file.remove("./datasets/cars_dataid_NA.Rda")
datasetsrel <- setdiff(datasetsrel, "cars")





load("./datasets/solar_flare_1_dataid_NA.Rda")
data1 <- dataset

load("./datasets/solar-flare_dataid_40686.Rda")
data2 <- dataset

dim(data1)
dim(data2)

head(data1)
head(data2)

# --> The same dataset, remove the first one ("solar_flare_1_dataid_NA.Rda").

file.remove("./datasets/solar_flare_1_dataid_NA.Rda")
datasetsrel <- setdiff(datasetsrel, "solar_flare_1")





load("./datasets/solar_flare_2_dataid_NA.Rda")
data1 <- dataset

load("./datasets/solar-flare_dataid_40686.Rda")
data2 <- dataset

dim(data1)
dim(data2)

sum(apply(data2, 1, function(x) paste(x, collapse="_")) %in% apply(data1, 1, function(x) paste(x, collapse="_")))

# --> These datasets are reasonably different.








# Remove datasets with categorical predictors with more than 50 categories:

alldatafiles <- list.files("./datasets")

for(i in seq(along=alldatafiles)) {
  
  load(paste("./datasets/", alldatafiles[i], sep=""))
  
  toomanylevels <- sapply(dataset, function(x) !is.numeric(x) & length(unique(x)) > 50)
  
  if(any(toomanylevels)) {
    
    datasub <- dataset[,toomanylevels,drop=FALSE]
    
    for(j in 1:ncol(datasub)) {
      eval(parse(text=paste0("dataset$", names(datasub)[j], " <- NULL")))
    }
    
    save(dataset, file=paste0("./datasets/", alldatafiles[i]))
    
  }
  
}



# Remove datasets with only one feature:

alldatafiles <- list.files("./datasets")

for(i in seq(along=alldatafiles)) {
  load(paste("./datasets/", alldatafiles[i], sep=""))
  
  if(ncol(dataset) <= 2)
    file.remove(paste("./datasets/", alldatafiles[i], sep=""))
  
}



# Remove datasets with missing values:

alldatafiles <- list.files("./datasets")

for(i in seq(along=alldatafiles)) {
  load(paste("./datasets/", alldatafiles[i], sep=""))
  
  if(any(is.na(dataset)))
    file.remove(paste("./datasets/", alldatafiles[i], sep=""))
  
}






# Apply interaction forests and random forests to each of the datasets
# to see whether any errors result:

library("diversityForest")
library("ranger")

alldatafiles <- list.files("./datasets")

error_iterations <- c()  # Vector to store iterations where errors occur

for(i in seq(along=alldatafiles)) {
  
  tryCatch({
    
    load(paste("./datasets/", alldatafiles[i], sep=""))
    
    
    ntrain <- floor(nrow(dataset)*0.8)
    trainind <- sample(1:nrow(dataset), size=ntrain)
    
    datatrain <- dataset[trainind,]
    datatest <- dataset[-trainind,]
    
    
    modeltemp <- interactionfor(dependent.variable.name = "ytarget", data = datatrain, 
                                num.trees = 20, importance = "none")
    
    preds <- predict(modeltemp, data=datatest)
    
    
    modeltemp <- ranger(dependent.variable.name = "ytarget", data = datatrain, 
                        num.trees = 200, importance = "none")
    
    preds <- predict(modeltemp, data=datatest)
    
    cat(paste0("Iteration: ", i), "\n")
    
  }, error = function(e) {
    # Record the iteration number if an error occurs
    error_iterations <<- c(error_iterations, i)
  })
  
}


# Print or inspect the error_iterations vector after the loop
print(error_iterations)

# --> No errors. --> All datasets are fine for use.






# Check the datasets for simulated datasets:


# The datasets
# "autoUniv-au4-2500_dataid_1548.Rda", "autoUniv-au6-1000_dataid_1555.Rda", 
# "autoUniv-au6-400_dataid_1551.Rda", "autoUniv-au6-750_dataid_1549.Rda",
# "autoUniv-au7-1100_dataid_1552.Rda", "autoUniv-au7-500_dataid_1554.Rda", 
# and "autoUniv-au7-700_dataid_1553.Rda" are simulated, but simulated realistically:

# AutoUniv is an advanced data generator for classifications tasks. 
# The aim is to reflect the nuances and heterogeneity of real data. 
# Data can be generated in .csv, ARFF or C4.5 formats.



# The dataset "balance-scale_dataid_11.Rda" is also simulated, but realistically:

# Citation from the paper "Modeling Cognitive Development on Balance
# Scale Phenomena":

# "The simulations reported here suggest that cascade-correlation networks can capture the
# main features of cognitive development on the balance scale."



# Unclear dataset "calendarDOW_dataid_40663.Rda":

# Look at the dataset:

load("./datasets/calendarDOW_dataid_40663.Rda")
boxplot(dataset[,-ncol(dataset)])

# --> Does not look simulated.



# Unclear dataset "diggle_table_a2_dataid_694.Rda":

# Look at the dataset:

load("./datasets/diggle_table_a2_dataid_694.Rda")
boxplot(dataset[,colnames(dataset)!="ytarget"])
boxplot(dataset[,colnames(dataset)!="ytarget"][,-(1:4)])

# --> Does not look simulated.




# There are four iris datasets:

load("./datasets/Iris_dataid_44344.Rda")
dim(dataset)
head(dataset)

load("./datasets/iris_dataid_61.Rda")
dim(dataset)
head(dataset)

load("./datasets/iris_reproduced_dataid_44154.Rda")
dim(dataset)
head(dataset)

load("./datasets/JuanFeldmanIris_dataid_42186.Rda")
dim(dataset)
head(dataset)

load("./datasets/MyIris_dataid_1413.Rda")
dim(dataset)
head(dataset)

# --> These are the same.
# --> Remove "iris_dataid_61.Rda", "iris_reproduced_dataid_44154.Rda",
# "JuanFeldmanIris_dataid_42186.Rda", and "MyIris_dataid_1413.Rda".

file.remove("./datasets/iris_dataid_61.Rda")
file.remove("./datasets/iris_reproduced_dataid_44154.Rda")
file.remove("./datasets/JuanFeldmanIris_dataid_42186.Rda")
file.remove("./datasets/MyIris_dataid_1413.Rda")




# The dataset "JapaneseVowels_dataid_375.Rda" contains time-series from
# nine individuals.

# --> Delete.

file.remove("./datasets/JapaneseVowels_dataid_375.Rda")




# The datasets "led7_dataid_40678.Rda" and "led24_dataid_40677.Rda" 
# (LED display data) are artificial, but they have been used frequently 
# for comparing classification algorithms.




# The "meta_..." datasets are all (almost) the same:


par(mfrow=c(2,2))

load("./datasets/meta_all.arff_dataid_275.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"], ylim=c(-5,5))

load("./datasets/meta_batchincremental.arff_dataid_276.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"], ylim=c(-5,5))

load("./datasets/meta_ensembles.arff_dataid_277.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"], ylim=c(-5,5))

load("./datasets/meta_instanceincremental.arff_dataid_278.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"], ylim=c(-5,5))

par(mfrow=c(1,1))

# --> Use "meta_ensembles.arff_dataid_277.Rda" because here the target
# variable does not feature very small classes.

# Delete the others:

file.remove("./datasets/meta_all.arff_dataid_275.Rda")
file.remove("./datasets/meta_batchincremental.arff_dataid_276.Rda")
file.remove("./datasets/meta_instanceincremental.arff_dataid_278.Rda")





# The datasets "MIP-2016-classification_dataid_41703.Rda" and
# "MIP-2016-PAR10-classification_dataid_41939.Rda" are (almost?)
# the same:

par(mfrow=c(1,2))

load("./datasets/MIP-2016-classification_dataid_41703.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"], ylim=c(-5,5))

load("./datasets/MIP-2016-PAR10-classification_dataid_41939.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"], ylim=c(-5,5))

par(mfrow=c(1,1))

# --> Delete "MIP-2016-PAR10-classification_dataid_41939.Rda":

file.remove("./datasets/MIP-2016-PAR10-classification_dataid_41939.Rda")






# Look at the "robot-failures-..." datasets:

par(mfrow=c(3,2))

load("./datasets/robot-failures-lp1_dataid_1516.Rda ")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"])

load("./datasets/robot-failures-lp2_dataid_1517.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"])

load("./datasets/robot-failures-lp3_dataid_1518.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"])

load("./datasets/robot-failures-lp4_dataid_1519.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"])

load("./datasets/robot-failures-lp5_dataid_1520.Rda")
dim(dataset)
table(dataset$ytarget)
boxplot(dataset[,colnames(dataset)!="ytarget"])

par(mfrow=c(1,1))

# --> The second and third, as well as the fourth and
# the fifth dataset have the same covariates; however, the
# target variable is different. --> Keep these datasets.




# As described in Section 3 of the paper "Time-Series Similarity Queries 
# Employing a Feature-Based Approach", the dataset "synthetic_control_dataid_377.Rda"
# is simulated, where it is not clear whether the simulation is realistic.

# -> Delete the dataset:

file.remove("./datasets/synthetic_control_dataid_377.Rda")



# The datasets "tae_dataid_48.Rda" and "teachingAssistant_dataid_1115.Rda"
# are the same:

load("./datasets/tae_dataid_48.Rda")
dim(dataset)
head(dataset)

load("./datasets/teachingAssistant_dataid_1115.Rda")
dim(dataset)
head(dataset)

# --> Delete "tae_dataid_48.Rda":

file.remove("./datasets/tae_dataid_48.Rda")




# The datasets "allhyper_dataid_NA.Rda" and "thyroid-allhyper_dataid_40475.Rda"
# are not the same, but seem to be versions of the same dataset:

load("./datasets/allhyper_dataid_NA.Rda")
dim(dataset)
head(dataset)

load("./datasets/thyroid-allhyper_dataid_40475.Rda")
dim(dataset)
head(dataset)

# --> Delete the smaller dataset "thyroid-allhyper_dataid_40475.Rda":

file.remove("./datasets/thyroid-allhyper_dataid_40475.Rda")





# The datasets "allhypo_dataid_NA.Rda", "thyroid-allhypo_dataid_40476.Rda",
# "allbp_dataid_40707.Rda", "thyroid-allbp_dataid_40474.Rda",
# "allrep_dataid_40708.Rda", "thyroid-allrep_dataid_40477.Rda",
# "thyroid-dis_dataid_40478.Rda", and "thyroid-ann_dataid_40497.Rda" seem to be 
# the same as the "allhyper" datasets above:

load("./datasets/allhypo_dataid_NA.Rda")
dim(dataset)
head(dataset)

load("./datasets/thyroid-allhypo_dataid_40476.Rda")
dim(dataset)
head(dataset)

load("./datasets/allbp_dataid_40707.Rda")
dim(dataset)
head(dataset)

load("./datasets/thyroid-allbp_dataid_40474.Rda")
dim(dataset)
head(dataset)

load("./datasets/allrep_dataid_40708.Rda")
dim(dataset)
head(dataset)

load("./datasets/thyroid-allrep_dataid_40477.Rda")
dim(dataset)
head(dataset)

load("./datasets/thyroid-dis_dataid_40478.Rda")
dim(dataset)
head(dataset)

load("./datasets/thyroid-ann_dataid_40497.Rda")
dim(dataset)
head(dataset)

# --> Remove them:

file.remove("./datasets/allhypo_dataid_NA.Rda")
file.remove("./datasets/thyroid-allhypo_dataid_40476.Rda")
file.remove("./datasets/allbp_dataid_40707.Rda")
file.remove("./datasets/thyroid-allbp_dataid_40474.Rda")
file.remove("./datasets/allrep_dataid_40708.Rda")
file.remove("./datasets/thyroid-allrep_dataid_40477.Rda")
file.remove("./datasets/thyroid-dis_dataid_40478.Rda")
file.remove("./datasets/thyroid-ann_dataid_40497.Rda")

# However, the dataset "thyroid-new_dataid_40682.Rda" is
# a different dataset:s

load("./datasets/thyroid-new_dataid_40682.Rda")
dim(dataset)
head(dataset)





# Unclear dataset "Touch2_dataid_42544.Rda":

# Look at the dataset:

load("./datasets/Touch2_dataid_42544.Rda")
dim(dataset)
head(dataset)
boxplot(dataset[,-ncol(dataset)])

# --> Does not look simulated.



# The dataset "waveform-5000_dataid_60.Rda" is simulated and its not clear
# whether it is simulated in a realistic way.

# --> Remove it:

file.remove("./datasets/waveform-5000_dataid_60.Rda")




# Make an overview table of the processed datasets:

alldatafiles <- list.files("./datasets")

labels <- filenames <- ""
ns <- ps <- n_cls <- 0

for(i in seq(along=alldatafiles)) {
  load(paste("./datasets/", alldatafiles[i], sep=""))
  
  filenames[i] <- alldatafiles[i]
  labels[i] <- strsplit(alldatafiles[i], split="_dataid")[[1]][1]
  
  ns[i] <- nrow(dataset)
  ps[i] <- ncol(dataset) - 1
  n_cls[i] <- length(unique(dataset$ytarget))
  
  if(i %% 10 == 0)
    cat(paste("Iteration:", i), "\n")
}

datainfo_preliminary <- data.frame(filename=filenames, label=labels, n=ns, p=ps, n_cl=n_cls)
datainfo_preliminary <- datainfo_preliminary[order(datainfo_preliminary$label),]








# Many data sets have factor variables that are not cast as factors in R.
# Manually inspect all datasets and convert factor variables that are not 
# cast as factors into factors (two further data sets are removed in the process):


datasetid <- 1

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

table(dataset$V1)




datasetid <- 2

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

tempinds <- which(sapply(dataset, function(x) length(unique(x))) <= 5 & names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.factor(dataset[, tempinds[i]])
}

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))




datasetid <- 3

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

table(dataset$BookID)
dataset$BookID <- as.factor(dataset$BookID)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 4

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 5

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 6

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)
sapply(dataset, is.ordered)

table(dataset$Age)

agefactor <- dataset$Age
dataset$Age <- factor(agefactor, levels = levels(agefactor), ordered = TRUE)

table(dataset$Time_of_survey)
levels(dataset$Time_of_survey)

factortemp <- dataset$Time_of_survey
dataset$Time_of_survey <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))




datasetid <- 7

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

dataset$Years_of_schooling <- as.numeric(as.character(dataset$Years_of_schooling))

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))




datasetid <- 8

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

tempinds <- which(sapply(dataset, function(x) length(unique(x))) <= 2 & names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.factor(dataset[, tempinds[i]])
}

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))




datasetid <- 9

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 10

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

tempinds <- which(sapply(dataset, function(x) length(unique(x))) <= 2 & names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.factor(dataset[, tempinds[i]])
}

dataset$body.style <- as.factor(dataset$body.style)
dataset$engine.type <- as.factor(dataset$engine.type)
dataset$make <- as.factor(dataset$make)
dataset$fuel.system <- as.factor(dataset$fuel.system)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))




datasetid <- 11

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 12

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 13

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 14

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 15

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 16

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 17

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 18

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 19

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 20

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 21

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

table(dataset$buying)
levels(dataset$buying)

factortemp <- dataset$buying
dataset$buying <- factor(factortemp, levels = rev(levels(factortemp)), ordered = TRUE)

factortemp <- dataset$maint
dataset$maint <- factor(factortemp, levels = rev(levels(factortemp)), ordered = TRUE)

factortemp <- dataset$lug_boot
dataset$lug_boot <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$safety
dataset$safety <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 22

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 23

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 24

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

factortemp <- dataset$cylinders
dataset$cylinders <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 25

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}

dataset$sex <- as.factor(dataset$sex)
dataset$cp <- as.factor(dataset$cp)
dataset$fbs <- as.factor(dataset$fbs)
dataset$restecg <- as.factor(dataset$restecg)
dataset$exang <- as.factor(dataset$exang)

dataset$slope <- as.factor(dataset$slope)
table(dataset$slope)

factortemp <- dataset$slope
dataset$slope <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

dataset$thal <- as.factor(dataset$thal)

sapply(dataset, class)

head(dataset)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 26

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 27

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

dataset$SEEDED <- as.factor(dataset$SEEDED)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 28

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 29

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# --> No information found on this dataset.
# --> Nothing changed.






datasetid <- 30

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 31

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}

dataset$Wife_religion <- as.factor(dataset$Wife_religion)
dataset$Wife_working <- as.factor(dataset$Wife_working)
dataset$Husband_occupation <- as.factor(dataset$Husband_occupation)
dataset$Media_exposure <- as.factor(dataset$Media_exposure)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 32

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# No information found on this dataset.
# Most variables have many numbers of different values.
# --> Most of the variables seem numeric.
# --> Leave the dataset as it is.






datasetid <- 33

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# Basically the same variables as in the previous dataset.
# --> Leave the dataset as it is.





datasetid <- 34

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# --> According to information found online the attributes in
# this dataset are numeric.






datasetid <- 35

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# No information found on this dataset.
# --> Leave the dataset as it is.

# Exceptions: Remove "ID" and "Project", as this variables are
# clearly not helpful:

dataset$ID <- NULL
dataset$Project <- NULL

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 36

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> The variables in this dataset have many different
# unique values, which is why the do not seem to be categorical.
# --> Leave the dataset as it is.





datasetid <- 37

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 38

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 39

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 40

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

head(dataset$P1stFixation)
head(dataset$P2stFixation)
head(dataset$nextWordRegress)

# --> The variables with fewer categories have "n", "Cnt", "count",
# or similar words in their names, which suggests that these
# are numbers and can thus be treated as numeric variables.





datasetid <- 41

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 42

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

table(dataset$columns)





datasetid <- 43

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 44

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 45

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 46

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

paste0(rep(c(0,1,2,6,7,8,9), each=3), rep(c("f", "m", "c"), times=3))





datasetid <- 47

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 48

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

table(dataset$year_zone)

dataset$year_zone <- factor(as.character(dataset$year_zone), levels = paste0(rep(c(0,1,2,6,7,8,9), each=3), 
                                                                             rep(c("f", "m", "c"), times=3)), 
                            ordered = TRUE)

factortemp <- dataset$year
head(factortemp)
dataset$year <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$damage_rankRJT
head(factortemp)
dataset$damage_rankRJT <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$damage_rankALL
head(factortemp)
dataset$damage_rankALL <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 49

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

dataset$hobby <- as.factor(dataset$hobby)
dataset$marital_status <- as.factor(dataset$marital_status)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 50

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# Some variables seem categorical, but leave it as it is because
# no information was found on this version of the dataset online.




datasetid <- 51

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

dataset$V2 <- as.factor(dataset$V2)
dataset$V3 <- as.factor(dataset$V3)
dataset$V6 <- as.factor(dataset$V6)
dataset$V7 <- as.factor(dataset$V7)
dataset$V9 <- as.factor(dataset$V9)

dataset$V11 <- as.factor(dataset$V11)

factortemp <- dataset$V11
levels(factortemp)
dataset$V11 <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

dataset$V12 <- as.factor(dataset$V12)
dataset$V13 <- as.factor(dataset$V13)

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 52

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# --> This dataset has a very similar name and almost the same variables as the
# previous dataset (seen on OpenML) and for the previous dataset all variables with
# few unique values were categorical covariates.
# --> Convert the variables with few unique variables into factors:

dataset$V2 <- as.factor(dataset$V2)
dataset$V3 <- as.factor(dataset$V3)
dataset$V5 <- as.factor(dataset$V5)
dataset$V6 <- as.factor(dataset$V6)
dataset$V8 <- as.factor(dataset$V8)
dataset$V10 <- as.factor(dataset$V10)
dataset$V11 <- as.factor(dataset$V11)
dataset$V12 <- as.factor(dataset$V12)

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 53

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))




datasetid <- 54

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 55

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> Given that some variables are factor and others are numeric, the variables
# likely already have the right scales of measurements.




datasetid <- 56

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> Given that some variables are factor and others are numeric, the variables
# likely already have the right scales of measurements.




datasetid <- 57

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> Given that some variables are factor and others are numeric, the variables
# likely already have the right scales of measurements.





datasetid <- 58

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> The variables are all binary and according to OpenML that are all Boolean.

# --> Convert them to categorical covariates:

tempinds <- which(sapply(dataset, function(x) length(unique(x))) <= 2 & names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.factor(dataset[, tempinds[i]])
}

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))




datasetid <- 59

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 60

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 61

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> Most variables are strings and thus factors. There are are a few variables,
# which are numbers and thus of "numeric" type.
# --> It is likely that the variables which are numbers should be treated as metric
# variables - if not, they likely would have been coded as strings.




datasetid <- 62

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> This is the same dataset as the previous one.
# --> Remove it:

file.remove("./datasets/lymphography_dataid_NA.Rda")




datasetid <- 63

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 64

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 65

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 66

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 67

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> Some variables have few unique values and thus could be categorical
# variables. However, according to OpenML: "The meaning of the features is mostly unknown."
# Thus, we leave the dataset as it is.




datasetid <- 68

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# --> According to OpenML: "The mfeatures represent 240 (15 x 16) pixel averages in 2 x 3 windows."
# --> These are to be treated as metric features.




datasetid <- 69

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 70

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# There are a number of variables in this dataset that have few unique values.
# However, there names suggest that they should be treated as metric variables.
# On OpenML it says: "One of a set of 6 datasets describing features of 
# handwritten numerals (0 - 9) extracted from a collection of Dutch utility maps. "





datasetid <- 71

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))





datasetid <- 72

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 73

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

lapply(dataset, table)

factortemp <- dataset$parents
levels(factortemp)
dataset$parents <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$has_nurs
levels(factortemp)
dataset$has_nurs <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$children
levels(factortemp)
dataset$children <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$housing
levels(factortemp)
dataset$housing <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$finance
levels(factortemp)
dataset$finance <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$social
levels(factortemp)
dataset$social <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 74

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))

# According to OpenML all variables are integers.
# --> Should be treated as numeric.





datasetid <- 75

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 76

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 77

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 78

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

dataset$island <- as.factor(dataset$island)
dataset$sex <- as.factor(dataset$sex)

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 79

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# According to OpenML, the variables which are numeric are either metric (like Age)
# or ranks.
# --> Leave the dataset as it is.





datasetid <- 80

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# Delete Label variable because this is a categorical variable
# with 27 categories (ther are only 27 observations in total):
dataset$Label <- NULL

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 81

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))






datasetid <- 82

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# --> Some of the variables are factors, but according to OpenML these
# variables all represent numbers.
# --> Change them all to numeric.

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(as.character(dataset[, tempinds[i]]))
}

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 83

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 84

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 85

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 86

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 87

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 88

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 89

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 90

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(as.character(dataset[, tempinds[i]]))
}

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

dataset$sex <- as.factor(dataset$sex)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 91

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 92

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 93

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

# --> Upon closer inspection it was seen that this dataset is the same 
# as the previous dataset (just with a different row order).
# --> Remove this dataset:

file.remove("./datasets/segmentation_dataid_NA.Rda")





datasetid <- 94

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 95

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# --> The variable are all factors. --> Transform them into factors:

tempinds <- which(sapply(dataset, function(x) length(unique(x))) <= 2 & names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.factor(dataset[, tempinds[i]])
}

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 96

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

dataset$Person

# --> This variable has too many factor levels. --> Remove:

dataset$Person <- NULL

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 97

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

factortemp <- dataset$Evolution
levels(factortemp)
dataset$Evolution <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$Previous_24_hour_flare_activity_code
levels(factortemp)
dataset$Previous_24_hour_flare_activity_code <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

dataset$C.class_flares_production_by_this_region <- as.numeric(as.character(dataset$C.class_flares_production_by_this_region))
dataset$M.class_flares_production_by_this_region <- as.numeric(as.character(dataset$M.class_flares_production_by_this_region))
dataset$X.class_flares_production_by_this_region <- as.numeric(as.character(dataset$X.class_flares_production_by_this_region))

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 98

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}


dataset$largest_spot_size <- as.factor(dataset$largest_spot_size)
dataset$spot_distribution <- as.factor(dataset$spot_distribution)
dataset$Activity <- as.factor(dataset$Activity)

factortemp <- dataset$Evolution
levels(factortemp)
dataset$Evolution <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

factortemp <- dataset$Previous_24_hour_flare_activity_code
levels(factortemp)
dataset$Previous_24_hour_flare_activity_code <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

dataset$Historically.complex <- as.factor(dataset$Historically.complex)
dataset$Did_region_become_historically_complex <- as.factor(dataset$Did_region_become_historically_complex)
dataset$Area <- as.factor(dataset$Area)
dataset$Area <- as.factor(dataset$Area)

table(dataset$Area_of_the_largest_spot)
dataset$Area_of_the_largest_spot <- NULL

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 99

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)
table(sapply(dataset, class))

tempinds <- which(names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.numeric(dataset[, tempinds[i]])
}

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# The variables in this dataset all have few unique values and 
# according to the UCI Machine Learning Repository
# (https://archive.ics.uci.edu/dataset/90/soybean+large)
# all are categorical.
# --> Transform them into categorical covariates:

tempinds <- which(sapply(dataset, function(x) length(unique(x))) <= 10 & names(dataset) != "ytarget")
for (i in seq(along = tempinds)) {
  dataset[, tempinds[i]] <- as.factor(dataset[, tempinds[i]])
}

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 100

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 101

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# The following variables are categorical according to OpenML:

dataset$V12 <- as.factor(dataset$V12)
dataset$V13 <- as.factor(dataset$V13)

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 102

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

dataset$ID <- NULL

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))





datasetid <- 103

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 104

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)




datasetid <- 105

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
table(sapply(dataset, class))

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# --> There are many variables with constant values in this dataset.
# The rest seem to be metric variables.
# --> Leave them as they are.






datasetid <- 106

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

sapply(dataset, function(x) length(unique(x)))
sort(sapply(dataset, function(x) length(unique(x))))

# No information found on this dataset online. 
# However, most variables have many unique values, which is
# why they are likely categorical. Online one variable "f1" has
# only two unique values, which is why this variable could be
# categorical (but not certain).
# --> Leave this dataset as it is.





datasetid <- 107

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 108

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)

dataset$ID <- NULL

# There are a couple of variables with missing values, as indicated
# by the category "NULL":

sapply(dataset, function(x) which(x=="NULL"))

# Remove observations with missing values:

dataset <- dataset[!(dataset$DataOut=="NULL" | dataset$Tools=="NULL"),]


dataset$Effort <- as.numeric(as.character(dataset$Effort))

factortemp <- dataset$IntComplx
levels(factortemp)
dataset$Evolution <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

dataset$DataFile <- as.numeric(as.character(dataset$DataFile))
dataset$DataEn <- as.numeric(as.character(dataset$DataEn))
dataset$DataOut <- as.numeric(as.character(dataset$DataOut))

table(dataset$ToolExpr)

# --> This variable has too many categories and these categories
# also don't seem to make a lot of sense.
# --> Delete the variable:

dataset$ToolExpr <- NULL

dataset$AppExpr <- as.factor(dataset$AppExpr)

factortemp <- dataset$AppExpr
levels(factortemp)
dataset$AppExpr <- factor(factortemp, levels = levels(factortemp), ordered = TRUE)

table(dataset$TeamSize)

# --> This variable has too many categories and these categories
# also don't seem to make a lot of sense as a predictor.
# --> Delete the variable:

dataset$TeamSize <- NULL

sapply(dataset, class)

save(dataset, file=paste0("./datasets/", datainfo_preliminary$filename[datasetid]))






datasetid <- 109

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)
table(sapply(dataset, class))





datasetid <- 110

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 111

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 112

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 113

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 114

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 115

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 116

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 117

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 118

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 119

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 120

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 121

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)





datasetid <- 122

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)






datasetid <- 123

paste0("./datasets/", datainfo_preliminary$filename[datasetid])
load(paste0("./datasets/", datainfo_preliminary$filename[datasetid]))

dim(dataset)
head(dataset)
sapply(dataset, class)











# From all datasets remove variables which are constant across all observations:


alldatafiles <- list.files("./datasets")

datainfo_preliminary <- datainfo_preliminary[datainfo_preliminary$filename %in% alldatafiles,]

for(i in seq(along=alldatafiles)) {
  
  load(paste("./datasets/", alldatafiles[i], sep=""))
  
  constbool <- sapply(dataset, function(x) length(unique(x)) == 1)
  
  if (any(constbool)) {
    dataset <- dataset[,!constbool]
    save(dataset, file=paste("./datasets/", alldatafiles[i], sep=""))
  }
  
}










# Factor variables with many categories are not ideal for random forests. 
# Examine all datasets that have factors with more than 15 categories to see 
# if these factors are ordered factors, so that we can convert them to ordered 
# factors, which are treated differently in the random forest variants than 
# unordered factors:

inds <- c()

for(i in seq(along=alldatafiles)) {
  
  load(paste("./datasets/", alldatafiles[i], sep=""))
  
  manycats <- sapply(dataset, function(x) is.factor(x) & length(unique(x)) > 15)
  
  if (any(manycats)) {
    inds <- c(inds, i)
  }
  
}


# "./datasets/analcatdata_challenger_dataid_462.Rda":

paste("./datasets/", alldatafiles[inds[1]], sep="")
load(paste("./datasets/", alldatafiles[inds[1]], sep=""))

dataset$Date <- NULL

save(dataset, file=paste("./datasets/", alldatafiles[inds[1]], sep=""))



# "./datasets/auto_dataid_NA.Rda":

paste("./datasets/", alldatafiles[inds[2]], sep="")
load(paste("./datasets/", alldatafiles[inds[2]], sep=""))

manycats <- sapply(dataset, function(x) is.factor(x) & length(unique(x)) > 15)
sum(manycats)

names(dataset)[manycats]

table(dataset$make)

# Online search revealed that make is unfortunately not an ordered factor.



# "./datasets/grub-damage_dataid_338.Rda":

paste("./datasets/", alldatafiles[inds[3]], sep="")
load(paste("./datasets/", alldatafiles[inds[3]], sep=""))

dim(dataset)

manycats <- sapply(dataset, function(x) is.factor(x) & length(unique(x)) > 15)
sum(manycats)

names(dataset)[manycats]

table(dataset$year_zone)
head(dataset$year_zone)




# "./datasets/soybean_dataid_NA.Rda":

paste("./datasets/", alldatafiles[inds[4]], sep="")
load(paste("./datasets/", alldatafiles[inds[4]], sep=""))

dim(dataset)

manycats <- sapply(dataset, function(x) is.factor(x) & length(unique(x)) > 15)
sum(manycats)

names(dataset)[manycats]




# "./datasets/teachingAssistant_dataid_1115.Rda":

paste("./datasets/", alldatafiles[inds[5]], sep="")
load(paste("./datasets/", alldatafiles[inds[5]], sep=""))

dim(dataset)

manycats <- sapply(dataset, function(x) is.factor(x) & length(unique(x)) > 15)
sum(manycats)

names(dataset)[manycats]

table(dataset$courseInstructor)
table(dataset$course)

# --> These are unfortunately no ordered factors.




# "./datasets/usp05_dataid_1047.Rda":

paste("./datasets/", alldatafiles[inds[6]], sep="")
load(paste("./datasets/", alldatafiles[inds[6]], sep=""))

dim(dataset)

manycats <- sapply(dataset, function(x) is.factor(x) & length(unique(x)) > 15)
sum(manycats)

names(dataset)[manycats]

table(dataset$Lang)
table(dataset$Tools)

# --> These are unfortunately no ordered factors.





# "./datasets/usp05_dataid_1047.Rda":

paste("./datasets/", alldatafiles[inds[7]], sep="")
load(paste("./datasets/", alldatafiles[inds[7]], sep=""))

dim(dataset)

manycats <- sapply(dataset, function(x) is.factor(x) & length(unique(x)) > 15)
sum(manycats)

names(dataset)[manycats]

table(dataset$country)

# --> These is unfortunately no ordered factor.











# Make an overview table of the processed datasets (Tables S9 to S11):

alldatafiles <- list.files("./datasets")

labels <- filenames <- ""
ns <- ps <- n_cls <- prop_cats <- 0

for(i in seq(along=alldatafiles)) {
  load(paste("./datasets/", alldatafiles[i], sep=""))
  
  filenames[i] <- alldatafiles[i]
  labels[i] <- strsplit(alldatafiles[i], split="_dataid")[[1]][1]
  
  ns[i] <- nrow(dataset)
  ps[i] <- ncol(dataset) - 1
  n_cls[i] <- length(unique(dataset$ytarget))
  prop_cats[i] <- mean(sapply(dataset[,names(dataset)!="ytarget"], function(x) length(unique(x)) <= 10))
  
  if(i %% 10 == 0)
    cat(paste("Iteration:", i), "\n")
}

datainfo <- data.frame(filename=filenames, label=labels, n=ns, p=ps, n_cl=n_cls, prop_cat=prop_cats)
datainfo <- datainfo[order(datainfo$label),]

save(datainfo, file="./datainfo.Rda")




load("./datainfo.Rda")

datatab <- data.frame(data.id=as.numeric(sapply(gsub(".Rda", "", datainfo$filename), function(x) strsplit(x, split="_dataid_")[[1]][2])),
                      label=datainfo$label,
                      n=datainfo$n,
                      p=datainfo$p, 
                      n_cl=datainfo$n_cl,
                      prop_cat=format(round(datainfo$prop_cat, 3), nsmall = 3))
rownames(datatab) <- NULL

# Replace NA values in 'data.id' with "$-$"
datatab$data.id <- ifelse(is.na(datatab$data.id), "$-$", datatab$data.id)

library("xtable")

head(datatab)

# Convert the data.frame to a LaTeX table with xtable
latex_table <- xtable(datatab, 
                      caption = "Overview of the datasets used in the real data analysis. The following information is provided: \\lq\\lq data.id'': OpenML ID of the datasets, is \\lq\\lq $--$'' for datasets from PMLB, \\lq\\lq label'': dataset label, \\lq\\lq n'': sample size, \\lq\\lq p'': number of covariates, \\lq\\lq n\\_cl'': number of outcome classes, \\lq\\lq prop\\_cat'': percentage of categorical covariates", 
                      label = "tab:data_overview",
                      digits=c(0, 0, NA, 0, 0, 0, 3)
)

# Save the LaTeX table to a file:
print(latex_table, type = "latex", file = "../../tables/TabS9S10S11.tex", include.rownames = FALSE)
