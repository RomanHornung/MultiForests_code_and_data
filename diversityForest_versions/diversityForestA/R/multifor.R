# -------------------------------------------------------------------------------
#   This file is part of 'diversityForestA'.
#
# 'diversityForestA' is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# 'diversityForestA' is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with 'diversityForestA'. If not, see <http://www.gnu.org/licenses/>.
#
#  NOTE: 'diversityForestA' is a fork of the popular R package 'ranger', written by Marvin N. Wright.
#  Most R and C++ code is identical with that of 'ranger'. The package 'diversityForestA'
#  was written by taking the original 'ranger' code and making any
#  changes necessary to implement diversity forests.
#
# -------------------------------------------------------------------------------

##' Implements multi forests as described in Hornung & Hapfelmeier (2024).
##' 
##' @title Construct a multi forest prediction rule and calculates multi-class variable importance scores as described in Hornung & Hapfelmeier (2024).
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit. Interaction terms supported only for numerical variables.
##' @param data Training data of class \code{data.frame}, or \code{matrix}, \code{dgCMatrix} (Matrix).
##' @param num.trees Number of trees. Default is 500.
##' @param importance Variable importance mode, one of the following: "both" (the default), "multiclass", "binary", "none". 
##' @param imp.min.nclass Minimum class size for the classes to be considered in importance mode 'best_min'. Default is 0.
##' @param write.forest Save \code{multifor.forest} object, required for prediction. Set to \code{FALSE} to reduce memory usage if no prediction intended.
##' @param probability Grow a probability forest as in Malley et al. (2012). Using this option (default is \code{TRUE}), class probability predictions are obtained.
##' @param min.node.size Minimal node size. Default 5 for probability and 1 for classification.
##' @param max.depth Maximal tree depth. A value of NULL or 0 (the default) corresponds to unlimited depth, 1 to tree stumps (1 split per tree).
##' @param replace Sample with replacement. 
##' @param sample.fraction Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement. This can be a vector of class-specific values. 
##' @param case.weights Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.
##' @param scale.permutation.importance Scale permutation importance by standard error as in (Breiman 2001). Only applicable if permutation variable importance mode selected. TO DO
##' @param keep.inbag Save how often observations are in-bag in each tree. 
##' @param inbag Manually set observations per tree. List of size num.trees, containing inbag counts for each observation. Can be used for stratified sampling.
##' @param holdout Hold-out mode. Hold-out all samples with case weight 0 and use these for variable importance and prediction error.
##' @param oob.error Compute OOB prediction error. Set to \code{FALSE} to save computation time, e.g. for large survival forests.
##' @param num.threads Number of threads. Default is number of CPUs available.
##' @param verbose Show computation status and estimated runtime.
##' @param seed Random seed. Default is \code{NULL}, which generates the seed from \code{R}. Set to \code{0} to ignore the \code{R} seed. 
##' @param dependent.variable.name Name of outcome variable, needed if no formula given.
##' @param mtry Number of candidate variables to sample for each split. Default is the (rounded down) square root of the number variables.
##' @param npervar Number of splits to sample per candidate variable. Default is 3.
##' @return Object of class \code{multifor} with elements
##'   \item{\code{forest}}{Saved forest (If write.forest set to TRUE). Note that the variable IDs in the \code{split.varIDs} object do not necessarily represent the column number in R.}
##'   \item{\code{predictions}}{Predicted classes/values, based on out-of-bag samples (classification and regression only).}
##'   \item{\code{variable.importance}}{Variable importance for each independent variable.}
##'   \item{\code{prediction.error}}{Overall out-of-bag prediction error. For classification this is the fraction of missclassified samples, for probability estimation the Brier score, for regression the mean squared error and for survival one minus Harrell's C-index.}
##'   \item{\code{r.squared}}{R squared. Also called explained variance or coefficient of determination (regression only). Computed on out-of-bag data.}
##'   \item{\code{confusion.matrix}}{Contingency table for classes and predictions based on out-of-bag samples (classification only).}
##'   \item{\code{unique.death.times}}{Unique death times (survival only).}
##'   \item{\code{chf}}{Estimated cumulative hazard function for each sample (survival only).}
##'   \item{\code{survival}}{Estimated survival function for each sample (survival only).}
##'   \item{\code{call}}{Function call.}
##'   \item{\code{num.trees}}{Number of trees.}
##'   \item{\code{num.independent.variables}}{Number of independent variables.}
##'   \item{\code{min.node.size}}{Value of minimal node size used.}
##'   \item{\code{treetype}}{Type of forest/tree. classification, regression or survival.}
##'   \item{\code{importance.mode}}{Importance mode used.}
##'   \item{\code{num.samples}}{Number of samples.}
##'   \item{\code{splitrule}}{Splitting rule.}
##'   \item{\code{replace}}{Sample with replacement.}
##'   \item{\code{nsplits}}{Value of \code{nsplits} used.}
##'   \item{\code{proptry}}{Value of \code{proptry} used.}
##' @examples
##' \dontrun{
##' 
##' ## Load package:
##' library("diversityForestA")
##' 
##' ## Set seed to obtain reproducible results:
##' set.seed(1234)
##'
##' ## Diversity forest with default settings (NOT recommended)
##' # Classification:
##' divfor(Species ~ ., data = iris, num.trees = 20)
##' # Regression:
##' iris2 <- iris; iris2$Species <- NULL; iris2$Y <- rnorm(nrow(iris2))
##' divfor(Y ~ ., data = iris2, num.trees = 20)
##' # Survival:
##' library("survival")
##' divfor(Surv(time, status) ~ ., data = veteran, num.trees = 20, respect.unordered.factors = "order")
##' # NOTE: num.trees = 20 is specified too small for practical 
##' # purposes - the prediction performance of the resulting 
##' # forest will be suboptimal!!
##' # In practice, num.trees = 500 (default value) or a 
##' # larger number should be used.
##' 
##' ## Diversity forest with specified values for nsplits and proptry (NOT recommended)
##' divfor(Species ~ ., data = iris, nsplits = 10, proptry = 0.4, num.trees = 20)
##' # NOTE again: num.trees = 20 is specified too small for practical purposes.
##' 
##' ## Applying diversity forest after optimizing the values of nsplits and proptry (recommended)
##' tuneres <- tunedivfor(formula = Species ~ ., data = iris, num.trees.pre = 20)
##' # NOTE: num.trees.pre = 20 is specified too small for practical 
##' # purposes - the out-of-bag error estimates of the forests 
##' # constructed during optimization will be much too variable!!
##' # In practice, num.trees.pre = 500 (default value) or a 
##' # larger number should be used.
##' divfor(Species ~ ., data = iris, nsplits = tuneres$nsplitsopt, 
##'   proptry = tuneres$proptryopt, num.trees = 20)
##' # NOTE again: num.trees = 20 is specified too small for practical purposes.
##' 
##' ## Prediction
##' train.idx <- sample(nrow(iris), 2/3 * nrow(iris))
##' iris.train <- iris[train.idx, ]
##' iris.test <- iris[-train.idx, ]
##' tuneres <- tunedivfor(formula = Species ~ ., data = iris.train, num.trees.pre = 20)
##' # NOTE again: num.trees.pre = 20 is specified too small for practical purposes.
##' rg.iris <- divfor(Species ~ ., data = iris.train, nsplits = tuneres$nsplitsopt, 
##'   proptry = tuneres$proptryopt, num.trees = 20)
##' # NOTE again: num.trees = 20 is specified too small for practical purposes.
##' pred.iris <- predict(rg.iris, data = iris.test)
##' table(iris.test$Species, pred.iris$predictions)
##' 
##' ## Variable importance
##' rg.iris <- divfor(Species ~ ., data = iris, importance = "permutation", num.trees = 20)
##' # NOTE again: num.trees = 20 is specified too small for practical purposes.
##' rg.iris$variable.importance
##' }
##'
##' @author Roman Hornung, Marvin N. Wright
##' @references
##' \itemize{
##'   \item Hornung, R. (2022). Diversity forests: Using split sampling to enable innovative complex split procedures in random forests. SN Computer Science 3(2):1, <\doi{10.1007/s42979-021-00920-1}>.
##'   \item Wright, M. N., Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. Journal of Statistical Software 77:1-17, <\doi{10.18637/jss.v077.i01}>.
##'   \item Breiman, L. (2001). Random forests. Machine Learning 45:5-32, <\doi{10.1023/A:1010933404324}>.
##'   \item Malley, J. D., Kruppa, J., Dasgupta, A., Malley, K. G., & Ziegler, A. (2012). Probability machines: consistent probability estimation using nonparametric learning machines. Methods of Information in Medicine 51:74-81, <\doi{10.3414/ME00-01-0052}>.
##'   \item Meinshausen (2006). Quantile Regression Forests. Journal of Machine Learning Research 7:983-999.
##'   }
##' @seealso \code{\link{predict.multifor}}
##' @encoding UTF-8
##' @useDynLib diversityForestA, .registration = TRUE
##' @importFrom Rcpp evalCpp
##' @import stats 
##' @import utils
##' @importFrom Matrix Matrix
##' @export
multifor <- function(formula = NULL, data = NULL, num.trees = 500,
                   importance = "both", imp.min.nclass = 0, write.forest = TRUE, probability = TRUE,
                   min.node.size = NULL, max.depth = NULL, replace = TRUE, 
                   sample.fraction = ifelse(replace, 1, 0.632), 
                   case.weights = NULL,
                   scale.permutation.importance = FALSE,
                   keep.inbag = FALSE, inbag = NULL, holdout = FALSE,
                   oob.error = TRUE,
                   num.threads = NULL,
                   verbose = TRUE, seed = NULL, 
                   dependent.variable.name = NULL, 
                   mtry = NULL, npervar = 3) {

  ## For multi forests we always order the categories of categorical variables:
  respect.unordered.factors <- "order"
  save.memory <- FALSE
  
  ## GenABEL GWA data
  if ("gwaa.data" %in% class(data)) {
    stop("Error: Ordering of SNPs currently not implemented for multi forests.")
  }
  
    snp.data <- as.matrix(0)

  ## Sparse matrix data
  if (inherits(data, "Matrix")) {
    if (!("dgCMatrix" %in% class(data))) {
      stop("Error: Currently only sparse data of class 'dgCMatrix' supported.")
    }
  
    if (!is.null(formula)) {
      stop("Error: Sparse matrices only supported with alternative interface. Use dependent.variable.name instead of formula.")
    }
  }
    
    if(is.null(mtry)) {
      mtry <- floor(sqrt(ncol(data) - 1))
    }
    
  ## Formula interface. Use whole data frame is no formula provided and depvarname given
  if (is.null(formula)) {
    if (is.null(dependent.variable.name)) {
      stop("Error: Please give formula or outcome variable name.")
    }
      response <- data[, dependent.variable.name, drop = TRUE]
    data.selected <- data
  } else {
    formula <- formula(formula)
    if (!inherits(formula, "formula")) {
      stop("Error: Invalid formula.")
    }
    data.selected <- parse.formula(formula, data, env = parent.frame())
    response <- data.selected[, 1]
  }
  
  ## Check missing values
  if (any(is.na(data.selected))) {
    offending_columns <- colnames(data.selected)[colSums(is.na(data.selected)) > 0]
    stop("Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  
  ## Outcome must be factor for multi forests:
  if (!is.factor(response)) {
    stop("Error: Outcome variable needs to be a factor.")
  }
  
  ## The outcome must have a maximum of 20 classes:
  maxlevels <- 20
  if (nlevels(droplevels(response)) > maxlevels)
    stop(paste0("Error: The outcome must have at most ", maxlevels, " classes (categories)."))
  
  ## Check response levels
    if (nlevels(response) != nlevels(droplevels(response))) {
      dropped_levels <- setdiff(levels(response), levels(droplevels(response)))
      warning("Dropped unused factor level(s) in outcome variable: ",
              paste0(dropped_levels, collapse = ", "), ".", call. = FALSE)
    }

  ## Treetype
    if (probability) {
      treetype <- 9
    } else {
      treetype <- 1
    }

  ## Dependent and status variable name. For non-survival dummy status variable name.
  if (!is.null(formula)) {
      dependent.variable.name <- names(data.selected)[1]
    independent.variable.names <- names(data.selected)[-1]
  } else {
    independent.variable.names <- colnames(data.selected)[colnames(data.selected) != dependent.variable.name]
  }
  
  ## Recode characters as factors and recode factors if 'order' mode
  if (!is.matrix(data.selected) && !inherits(data.selected, "Matrix")) {
    character.idx <- sapply(data.selected, is.character)

      ## Recode characters and unordered factors
      names.selected <- names(data.selected)
      ordered.idx <- sapply(data.selected, is.ordered)
      factor.idx <- sapply(data.selected, is.factor)
      independent.idx <- names.selected != dependent.variable.name
      recode.idx <- independent.idx & (character.idx | (factor.idx & !ordered.idx))
      
      ## Numeric response
        num.response <- as.numeric(response)

      ## Recode each column
      data.selected[recode.idx] <- lapply(data.selected[recode.idx], function(x) {
        if (!is.factor(x)) {
          x <- as.factor(x)
        } 
        
        if (nlevels(response) > 2) {
          levels.ordered <- pca.order(y = response, x = x)
        } else {
          ## Order factor levels by mean response
          means <- sapply(levels(x), function(y) {
            mean(num.response[x == y])
          })
          levels.ordered <- as.character(levels(x)[order(means)])
        }
        
        ## Return reordered factor
        factor(x, levels = levels.ordered, ordered = TRUE, exclude = NULL)
      })
      
      ## Save levels
      covariate.levels <- lapply(data.selected[independent.idx], levels)
  }
  
  ## Input data and variable names, create final data matrix
  if (is.matrix(data.selected) || inherits(data.selected, "Matrix")) {
    data.final <- data.selected
  } else {
    data.final <- data.matrix(data.selected)
  }
  variable.names <- colnames(data.final)
  
    all.independent.variable.names <- independent.variable.names

  ## Error if no covariates
  if (length(all.independent.variable.names) < 1) {
    stop("Error: No covariates found.")
  }
  
  ## Number of trees
  if (!is.numeric(num.trees) || num.trees < 1) {
    stop("Error: Invalid value for num.trees.")
  }
 
  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  
  ## Keep inbag
  if (!is.logical(keep.inbag)) {
    stop("Error: Invalid value for keep.inbag")
  }
  
  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) || num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }
  
  ## Minumum node size
  if (is.null(min.node.size)) {
    min.node.size <- 0
  } else if (!is.numeric(min.node.size) || min.node.size < 0) {
    stop("Error: Invalid value for min.node.size")
  }
  
  ## Tree depth
  if (is.null(max.depth)) {
    max.depth <- 0
  } else if (!is.numeric(max.depth) || max.depth < 0) {
    stop("Error: Invalid value for max.depth. Please give a positive integer.")
  }
  
  ## Sample fraction
  if (!is.numeric(sample.fraction)) {
    stop("Error: Invalid value for sample.fraction. Please give a value in (0,1] or a vector of values in [0,1].")
  }
  if (length(sample.fraction) > 1) {
    if (any(sample.fraction < 0) || any(sample.fraction > 1)) {
      stop("Error: Invalid value for sample.fraction. Please give a value in (0,1] or a vector of values in [0,1].")
    }
    if (sum(sample.fraction) <= 0) {
      stop("Error: Invalid value for sample.fraction. Sum of values must be >0.")
    }
    if (length(sample.fraction) != nlevels(response)) {
      stop("Error: Invalid value for sample.fraction. Expecting ", nlevels(response), " values, provided ", length(sample.fraction), ".")
    }
    if (!replace & any(sample.fraction * length(response) > table(response))) {
      idx <- which(sample.fraction * length(response) > table(response))[1]
      stop("Error: Not enough samples in class ", names(idx), 
           "; available: ", table(response)[idx], 
           ", requested: ", (sample.fraction * length(response))[idx], ".")
    }
    if (!is.null(case.weights)) {
      stop("Error: Combination of case.weights and class-wise sampling not supported.")
    }
  } else {
    if (sample.fraction <= 0 || sample.fraction > 1) {
      stop("Error: Invalid value for sample.fraction. Please give a value in (0,1] or a vector of values in [0,1].")
    }
  }
  
  ## Importance mode To Do:
  if (is.null(importance) || importance == "none") {
    importance.mode <- 0
  } else if (importance == "both") {
    importance.mode <- 6
  } else if (importance == "multiclass") {
    importance.mode <- 7
  } else if (importance == "binary") {
    importance.mode <- 8
  } else {
    stop("Error: Importance mode not supported for multi forests.")
  }
  
  ## Case weights: NULL for no weights asdf
  if (is.null(case.weights)) {
    case.weights <- c(0,0)
    use.case.weights <- FALSE
    if (holdout) {
      stop("Error: Case weights required to use holdout mode.")
    }
  } else {
    use.case.weights <- TRUE
    
    ## Sample from non-zero weights in holdout mode
    if (holdout) {
      sample.fraction <- sample.fraction * mean(case.weights > 0)
    }
    
    if (!replace && sum(case.weights > 0) < sample.fraction * nrow(data.final)) {
      stop("Error: Fewer non-zero case weights than observations to sample.")
    }
  }
  
  ## Manual inbag selection
  if (is.null(inbag)) {
    inbag <- list(c(0,0))
    use.inbag <- FALSE
  } else if (is.list(inbag)) {
    use.inbag <- TRUE
    if (use.case.weights) {
      stop("Error: Combination of case.weights and inbag not supported.")
    }
    if (length(sample.fraction) > 1) {
      stop("Error: Combination of class-wise sampling and inbag not supported.")
    }
    if (length(inbag) != num.trees) {
      stop("Error: Size of inbag list not equal to number of trees.")
    }
  } else {
    stop("Error: Invalid inbag, expects list of vectors of size num.trees.")
  }
  
  ## Splitting rule
      splitrule <- "gini"
  
  ## Prediction mode always false. Use predict.multifor() method.
  prediction.mode <- FALSE
  predict.all <- FALSE
  prediction.type <- 1
  
  ## No loaded forest object
  loaded.forest <- list()
  
  ## Use sparse matrix
  if ("dgCMatrix" %in% class(data.final)) {
    sparse.data <- data.final
    data.final <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.data <- Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
  }
  
    order.snps <- TRUE
  
  ## Clean up
  rm("data.selected")

  metricind <- which(apply(data.final, 2, function(x) length(unique(x)) >= length(unique(data.final[,dependent.variable.name]))))
  metricind <- setdiff(metricind, which(colnames(data.final)==dependent.variable.name))
  metricind <- metricind - 1

  # indsdatasets <- which(sapply(1:nrow(results_wide_sum), function(x) results_wide_sum[x,]$p)>10)
  
  ## Call C++:
  result <- divforCpp(treetype, dependent.variable.name, data.final, variable.names, mtry=0,
                      num.trees, verbose, seed, num.threads, write.forest, importance.mode,
                      min.node.size, split_select_weights=list(c(0,0)), use_split_select_weights=FALSE,
                      always_split_variable_names=c("0", "0"), use_always_split_variable_names=FALSE,
                      status_variable_name="", prediction.mode, loaded.forest, snp.data,
                      replace, probability, unordered_variable_names=c("0", "0"), use_unordered_variable_names=FALSE, 
                      save.memory, splitrule_r=1, case.weights, use.case.weights, class_weights=rep(1, nlevels(response)), 
                      predict.all, keep.inbag, sample.fraction, alpha=0.5, minprop=0.1, holdout, prediction.type, 
                      num_random_splits=npervar, sparse.data, use.sparse.data, order.snps, oob.error, max.depth, 
                      inbag, use.inbag, nsplits=mtry, npairs=0, proptry=0, divfortype=3, promispairs=list(0,0), 
                      eim_mode=0, metricind=metricind)
  
  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }
  
  ## Prepare results
  if (importance.mode != 0) {
    if (importance.mode == 6 || importance.mode == 7) {
    muw.varimps.all <- rep(NA, length(all.independent.variable.names))
    muw.varimps.all[match(metricind+1, which(colnames(data.final)!=dependent.variable.name))] <- result$var.imp.multiclass
    names(muw.varimps.all) <- all.independent.variable.names
    result$var.imp.multiclass <- muw.varimps.all
    }
    if (importance.mode == 6 || importance.mode == 8)
    names(result$var.imp.binary) <- all.independent.variable.names
  }

  ## Set predictions
  if (treetype == 1 && is.factor(response) && oob.error) {
    result$predictions <- integer.to.factor(result$predictions,
                                            levels(response))
    true.values <- integer.to.factor(unlist(data.final[, dependent.variable.name]),
                                     levels(response))
    result$confusion.matrix <- table(true.values, result$predictions, 
                                     dnn = c("true", "predicted"), useNA = "ifany")
  } else if (treetype == 9 && !is.matrix(data) && oob.error) {
    if (is.list(result$predictions)) {
      result$predictions <- do.call(rbind, result$predictions)
    } 
    if (is.vector(result$predictions)) {
      result$predictions <- matrix(result$predictions, nrow = 1)
    }
    
    ## Set colnames and sort by levels
    colnames(result$predictions) <- unique(response)
    if (is.factor(response)) {
      result$predictions <- result$predictions[, levels(droplevels(response)), drop = FALSE]
    }
  }
  
  ## Splitrule
  result$splitrule <- splitrule
  
  ## Set treetype
  if (treetype == 1) {
    result$treetype <- "Classification"
  } else if (treetype == 9) {
    result$treetype <- "Probability estimation"
  }

  result$call <- sys.call()
  result$importance.mode <- importance
  result$num.samples <- nrow(data.final)
  result$replace <- replace
  
  ## Write forest object
  if (write.forest) {
      result$forest$levels <- levels(response)

    result$forest$independent.variable.names <- independent.variable.names
    result$forest$treetype <- result$treetype
    class(result$forest) <- "multifor.forest"
    
    ## In 'ordered' mode, save covariate levels
    if (respect.unordered.factors == "order" && !is.matrix(data)) {
      result$forest$covariate.levels <- covariate.levels
    }
  }
  
  class(result) <- "multifor"

  result$mtry <- NULL
  
  return(result)
}
