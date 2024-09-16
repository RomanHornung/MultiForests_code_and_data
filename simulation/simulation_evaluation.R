####################################################################################

# NOTE: Before the code can be executed, the R working directory *MUST* 
# be set to the directory 'MultiForests_code_and_data/simulation' (that is
# the directory in which this R script is contained):

# Remove '#' from the line below and replace 'here/is/my/path/' by the path
# to the directory 'MultiForests_code_and_data/simulation':

# setwd("here/is/my/path/")

####################################################################################

# Load and pre-process the results:
#####################################

load("./intermediate_results/scenariogrid_simulation.Rda")
load("./intermediate_results/results_simulation.Rda")

reorderind <- order(scenariogrid$n, scenariogrid$K, scenariogrid$itind)
scenariogrid <- scenariogrid[reorderind,]
Results <- Results[reorderind]

scenariogrid$seed <- scenariogrid$settingid <- NULL


results <- scenariogrid[rep(1:nrow(scenariogrid), each=length(Results[[1]])),]
results$method_all <- rep(names(Results[[1]]), times=nrow(scenariogrid))

Results <- unlist(Results, recursive = FALSE)

resultsall <- results[rep(1:nrow(results), times=sapply(Results, length)),]
resultsall$rank <- unlist(lapply(Results, function(x) rank(-x)))
resultsall$vim <- unlist(Results)
rownames(resultsall) <- 1:nrow(resultsall)



resultsall$method_all <- factor(resultsall$method_all, levels=c("perm", "gini_corr", "bin_wgini_wsquared", "bin_wgini_wosquared", "bin_wogini_wsquared",
                                                                "bin_wogini_wosquared", "muw_wgini_wsquared", "muw_wgini_wosquared", "muw_wogini_wsquared", 
                                                                "muw_wogini_wosquared", "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared", 
                                                                "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared"))






# Function used for saving a table produced in R as a LaTeX table:

# Input parameters:

# combined_stats - object containing the table (has to be of a very specific
#                  format, see below)
# filename       - name of the file to which the table should be saved, can include
#                  path to that file
# table_title    - title to include in the produced LaTeX table

# Output:

# NULL

prepare_and_save_latex <- function(combined_stats, filename, table_title = NULL) {
  combined_stats$n <- as.integer(combined_stats$n) # Convert n to integer to remove decimal points
  
  # Check if a table title is provided and set the caption and caption.placement accordingly
  if (!is.null(table_title)) {
    latex_table <- xtable(combined_stats, digits=0, auto=FALSE, include.rownames=FALSE, caption = table_title)
    caption_placement <- "top"
  } else {
    latex_table <- xtable(combined_stats, digits=0, auto=FALSE, include.rownames=FALSE)
    caption_placement <- NULL  # No title, so no caption placement needed
  }
  
  # Customize the LaTeX output to add horizontal lines after each sample size
  print(latex_table, 
        type = "latex", 
        file = filename,
        include.rownames=FALSE,  # Ensure row names are not included
        add.to.row = list(pos = list(which(diff(combined_stats$n) != 0)), 
                          command = "\\hline "),
        caption.placement = caption_placement) # Conditionally place the caption based on if a title was provided
}












# Tables S1 to S8: Mean AUC values with 95% confidence intervals
################################################################


library("dplyr")
library("xtable")
library("stringr")
library("forcats")


resultstemp <- resultsall[resultsall$K==4,]
resultstemp$K <- NULL


# Provided AUC calculation function
auroc <- function(score, bool) {
  n1 <- sum(!bool)
  n2 <- sum(bool)
  U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

calculate_mean_ci_l <- function(auc_values) {
  mean_auc <- mean(auc_values)
  sd_auc <- sd(auc_values)
  n <- length(auc_values)
  se_auc <- sd_auc / sqrt(n)
  error_margin <- qnorm(0.975) * se_auc
  lower_ci <- mean_auc - error_margin
  return(lower_ci)
}

calculate_mean_ci_u <- function(auc_values) {
  mean_auc <- mean(auc_values)
  sd_auc <- sd(auc_values)
  n <- length(auc_values)
  se_auc <- sd_auc / sqrt(n)
  error_margin <- qnorm(0.975) * se_auc
  upper_ci <- mean_auc + error_margin
  upper_ci
}


# Compute AUC for specified groups and aggregate results
resultsK4_vsnoise <- resultstemp %>%
  filter(!str_detect(method_all, "muw_m_bin")) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(1:50, 51:53)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_2 = auroc(vim[c(1:50, 54:56)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_3 = auroc(vim[ c(1:50, 57:59)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(1:50, 60:62)], c(rep(FALSE, 50), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4) %>%
  rename_with(
    ~ c("two_gr", "cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_drop(method_all)) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK4_vsnoise, "../tables/TabS1.tex", table_title = "Influential vs noise variables, C = 4")




# Compute AUC for specified groups and aggregate results
resultsK4_vstwo_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 54:56)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wgini)" = "muw_m_bin_wgini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wgini)" = "muw_m_bin_wgini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wogini)" = "muw_m_bin_wogini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wogini)" = "muw_m_bin_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK4_vstwo_gr, "../tables/TabS4.tex", table_title = "cl\\_sp vs two\\_gr, C = 4")










resultstemp <- resultsall[resultsall$K==6,]
resultstemp$K <- NULL



# Compute AUC for specified groups and aggregate results
resultsK6_vsnoise <- resultstemp %>%
  filter(!str_detect(method_all, "muw_m_bin")) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(1:50, 51:53)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_2 = auroc(vim[c(1:50, 54:56)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_3 = auroc(vim[ c(1:50, 57:59)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(1:50, 60:62)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_5 = auroc(vim[c(1:50, 63:65)], c(rep(FALSE, 50), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    mean_auc_5 = mean(auc_5),
    ci_auc_5_l = calculate_mean_ci_l(auc_5),
    ci_auc_5_u = calculate_mean_ci_u(auc_5),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u),
    auc_5 = sprintf("%.2f [%.2f, %.2f]", mean_auc_5, ci_auc_5_l, ci_auc_5_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4, auc_5) %>%
  rename_with(
    ~ c("two_gr", "thr_gr", "cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_drop(method_all)) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK6_vsnoise, "../tables/TabS2.tex", table_title = "Influential vs noise variables, C = 6")





# Compute AUC for specified groups and aggregate results
resultsK6_vstwo_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wgini)" = "muw_m_bin_wgini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wgini)" = "muw_m_bin_wgini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wogini)" = "muw_m_bin_wogini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wogini)" = "muw_m_bin_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK6_vstwo_gr, "../tables/TabS5.tex", table_title = "cl\\_sp vs two\\_gr, C = 6")




# Compute AUC for specified groups and aggregate results
resultsK6_vsthr_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wgini)" = "muw_m_bin_wgini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wgini)" = "muw_m_bin_wgini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wogini)" = "muw_m_bin_wogini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wogini)" = "muw_m_bin_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK6_vsthr_gr, "../tables/TabS6.tex", table_title = "cl\\_sp vs thr\\_gr, C = 6")







resultstemp <- resultsall[resultsall$K==10,]
resultstemp$K <- NULL



# Compute AUC for specified groups and aggregate results
resultsK10_vsnoise <- resultstemp %>%
  filter(!str_detect(method_all, "muw_m_bin")) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(1:50, 51:53)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_2 = auroc(vim[c(1:50, 54:56)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_3 = auroc(vim[ c(1:50, 57:59)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(1:50, 60:62)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_5 = auroc(vim[c(1:50, 63:65)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_6 = auroc(vim[c(1:50, 66:68)], c(rep(FALSE, 50), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    mean_auc_5 = mean(auc_5),
    ci_auc_5_l = calculate_mean_ci_l(auc_5),
    ci_auc_5_u = calculate_mean_ci_u(auc_5),
    mean_auc_6 = mean(auc_6),
    ci_auc_6_l = calculate_mean_ci_l(auc_6),
    ci_auc_6_u = calculate_mean_ci_u(auc_6),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u),
    auc_5 = sprintf("%.2f [%.2f, %.2f]", mean_auc_5, ci_auc_5_l, ci_auc_5_u),
    auc_6 = sprintf("%.2f [%.2f, %.2f]", mean_auc_6, ci_auc_6_l, ci_auc_6_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4, auc_5, auc_6) %>%
  rename_with(
    ~ c("two_gr", "thr_gr", "cl_as_1", "cl_as_2", "cl_as_3", "cl_as_4"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_drop(method_all)) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK10_vsnoise, "../tables/TabS3.tex", table_title =
                         "Influential vs noise variables, C = 10")







# Compute AUC for specified groups and aggregate results
resultsK10_vstwo_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(51:53, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3", "cl_as_4"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wgini)" = "muw_m_bin_wgini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wgini)" = "muw_m_bin_wgini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wogini)" = "muw_m_bin_wogini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wogini)" = "muw_m_bin_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK10_vstwo_gr, "../tables/TabS7.tex", table_title = "cl\\_sp vs two\\_gr, C = 10")





# Compute AUC for specified groups and aggregate results
resultsK10_vsthr_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(54:56, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3", "cl_as_4"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "perm" = "perm",
                                 "gini_corr" = "gini_corr",
                                 "discr. VIM (wsquared_wgini)" = "bin_wgini_wsquared",
                                 "discr. VIM (wosquared_wgini)" = "bin_wgini_wosquared",
                                 "discr. VIM (wsquared_wogini)" = "bin_wogini_wsquared",
                                 "discr. VIM (wosquared_wogini)" = "bin_wogini_wosquared",
                                 "multi-class VIM (wsquared_wgini)" = "muw_wgini_wsquared",
                                 "multi-class VIM (wosquared_wgini)" = "muw_wgini_wosquared",
                                 "multi-class VIM (wsquared_wogini)" = "muw_wogini_wsquared",
                                 "multi-class VIM (wosquared_wogini)" = "muw_wogini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wgini)" = "muw_m_bin_wgini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wgini)" = "muw_m_bin_wgini_wosquared",
                                 "multi-class VIM - discr. VIM (wsquared_wogini)" = "muw_m_bin_wogini_wsquared",
                                 "multi-class VIM - discr. VIM (wosquared_wogini)" = "muw_m_bin_wogini_wosquared"
  ))

prepare_and_save_latex(resultsK10_vsthr_gr, "../tables/TabS8.tex", table_title = "cl\\_sp vs thr\\_gr, C = 10")














# Figure 2: Mean AUC values per considered sample size and method for C = 4.
#############################################################################


resultsall$n <- factor(paste0("n = ", resultsall$n), levels=c("n = 100", "n = 500", "n = 1000", "n = 2000"))

library("tidyr")

resultstemp <- resultsall[resultsall$K==4,]
resultstemp$K <- NULL


# Compute AUC for specified groups and aggregate results
resultsK4_vstwo_gr <- resultstemp %>%
  filter(!(method_all %in% c("bin_wgini_wsquared", "bin_wgini_wosquared", 
                             "bin_wogini_wsquared", "bin_wogini_wosquared"))) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 54:56)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK4_vstwo_gr %>%
  pivot_longer(cols = starts_with("cl_as_"), names_to = "cl_as", values_to = "value") %>%
  mutate(type = case_when(
    method_all %in% c("perm", "gini_corr")  ~ "conventional",
    str_detect(method_all, "muw_m_") ~ "multi-class diff.",
    str_detect(method_all, "muw_") ~ "multi-class",
    TRUE ~ NA_character_
  ))

results_long <- results_long %>%
  mutate(method = case_when(
    method_all %in% c("muw_wgini_wsquared", "muw_m_bin_wgini_wsquared") ~ "wsquared_wgini",
    method_all %in% c("muw_wgini_wosquared", "muw_m_bin_wgini_wosquared") ~ "wosquared_wgini",
    method_all %in% c("muw_wogini_wsquared", "muw_m_bin_wogini_wsquared") ~ "wsquared_wogini",
    method_all %in% c("muw_wogini_wosquared", "muw_m_bin_wogini_wosquared") ~ "wosquared_wogini",
    method_all == "perm" ~ "perm",
    method_all == "gini_corr" ~ "gini_corr"
  ),
  method = factor(method, levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini", "perm", "gini_corr")))

library("ggplot2")

p <- ggplot(results_long, aes(x = cl_as, y = value, group = method_all, color = method, linetype = type)) +
  geom_line() +
  geom_point(aes(shape = method), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=15),
        strip.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) + 
  labs(color="Method", shape="Method", linetype="Type", y="AUC") +
  scale_linetype_manual(values = c("conventional" = "solid", "multi-class diff." = 
                                     "dashed", "multi-class" = "dotdash")) +
  scale_shape_manual(values=c("wsquared_wgini" = 1, "wosquared_wgini" = 2, "wsquared_wogini" = 3, 
                              "wosquared_wogini" = 4, "perm" = 5, "gini_corr" = 6)) +
  scale_x_discrete(labels = c("cl_as_1" = expression(X[cl_as[1]]), "cl_as_2" = 
                                expression(X[cl_as[2]]), "cl_as_3" = expression(X[cl_as[3]]), 
                              "cl_as_4" = expression(X[cl_as[4]]))) +
  theme(legend.position = "right") +
  facet_wrap(~ n) +  # Add facet_wrap for separate plots by "n"
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))

p

# Figure 2:

ggsave("../figures/Fig2.eps", width=10, height=6)








# Figure 3: Mean AUC values per considered sample size and method for C = 6.
#############################################################################

# Calculate for K=6 and save
resultstemp <- resultsall[resultsall$K==6,]
resultstemp$K <- NULL

# Compute AUC for specified groups and aggregate results
resultsK6_vstwo_gr <- resultstemp %>%
  filter(!(method_all %in% c("bin_wgini_wsquared", "bin_wgini_wosquared", 
                             "bin_wogini_wsquared", "bin_wogini_wosquared"))) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK6_vstwo_gr %>%
  pivot_longer(cols = starts_with("cl_as_"), names_to = "cl_as", values_to = "value") %>%
  mutate(type = case_when(
    method_all %in% c("perm", "gini_corr")  ~ "conventional",
    str_detect(method_all, "muw_m_") ~ "multi-class diff.",
    str_detect(method_all, "muw_") ~ "multi-class",
    TRUE ~ NA_character_
  ))

results_long <- results_long %>%
  mutate(method = case_when(
    method_all %in% c("muw_wgini_wsquared", "muw_m_bin_wgini_wsquared") ~ "wsquared_wgini",
    method_all %in% c("muw_wgini_wosquared", "muw_m_bin_wgini_wosquared") ~ "wosquared_wgini",
    method_all %in% c("muw_wogini_wsquared", "muw_m_bin_wogini_wsquared") ~ "wsquared_wogini",
    method_all %in% c("muw_wogini_wosquared", "muw_m_bin_wogini_wosquared") ~ "wosquared_wogini",
    method_all == "perm" ~ "perm",
    method_all == "gini_corr" ~ "gini_corr"
  ),
  method = factor(method, levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", 
                                   "wosquared_wogini", "perm", "gini_corr")))

p1 <- ggplot(results_long, aes(x = cl_as, y = value, group = method_all, color = method, linetype = type)) +
  geom_line() +
  geom_point(aes(shape = method), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=17),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)
  ) + 
  ylab("AUC") +
  scale_linetype_manual(values = c("conventional" = "solid", "multi-class diff." = "dashed", 
                                   "multi-class" = "dotdash")) +
  scale_shape_manual(values=c("wsquared_wgini" = 1, "wosquared_wgini" = 2, "wsquared_wogini" = 3, 
                              "wosquared_wogini" = 4, "perm" = 5, "gini_corr" = 6)) +
  scale_x_discrete(labels = c("cl_as_1" = expression(X[cl_as[1]]), "cl_as_2" = expression(X[cl_as[2]]), 
                              "cl_as_3" = expression(X[cl_as[3]]), "cl_as_4" = expression(X[cl_as[4]]))) +
  theme(legend.position = "none") +
  ggtitle(expression("Comparison with" ~ X[two_gr])) + 
  facet_wrap(~ n)  # Add facet_wrap for separate plots by "n"




# Compute AUC for specified groups and aggregate results
resultsK6_vsthr_gr <- resultstemp %>%
  filter(!(method_all %in% c("bin_wgini_wsquared", "bin_wgini_wosquared", 
                             "bin_wogini_wsquared", "bin_wogini_wosquared"))) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK6_vsthr_gr %>%
  pivot_longer(cols = starts_with("cl_as_"), names_to = "cl_as", values_to = "value") %>%
  mutate(type = case_when(
    method_all %in% c("perm", "gini_corr")  ~ "conventional",
    str_detect(method_all, "muw_m_") ~ "multi-class diff.",
    str_detect(method_all, "muw_") ~ "multi-class",
    TRUE ~ NA_character_
  ))

results_long <- results_long %>%
  mutate(method = case_when(
    method_all %in% c("muw_wgini_wsquared", "muw_m_bin_wgini_wsquared") ~ "wsquared_wgini",
    method_all %in% c("muw_wgini_wosquared", "muw_m_bin_wgini_wosquared") ~ "wosquared_wgini",
    method_all %in% c("muw_wogini_wsquared", "muw_m_bin_wogini_wsquared") ~ "wsquared_wogini",
    method_all %in% c("muw_wogini_wosquared", "muw_m_bin_wogini_wosquared") ~ "wosquared_wogini",
    method_all == "perm" ~ "perm",
    method_all == "gini_corr" ~ "gini_corr"
  ),
  method = factor(method, levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", 
                                   "wosquared_wogini", "perm", "gini_corr")))

p2 <- ggplot(results_long, aes(x = cl_as, y = value, group = method_all, color = method, linetype = type)) +
  geom_point(aes(shape = method), size=2.5) +
  geom_line() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=17),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16)
  ) + 
  labs(color="Method", shape="Method", linetype="Type", y="AUC") +
  scale_shape_manual(values=c("wsquared_wgini" = 1, "wosquared_wgini" = 2, "wsquared_wogini" = 3, 
                              "wosquared_wogini" = 4, "perm" = 5, "gini_corr" = 6)) +
  scale_linetype_manual(values = c("conventional" = "solid", "multi-class diff." = "dashed", 
                                   "multi-class" = "dotdash")) +
  scale_x_discrete(labels = c("cl_as_1" = expression(X[cl_as[1]]), "cl_as_2" = 
                                expression(X[cl_as[2]]), "cl_as_3" = expression(X[cl_as[3]]), 
                              "cl_as_4" = expression(X[cl_as[4]]))) +
  theme(legend.position = "right") +
  ggtitle(expression("Comparison with" ~ X[thr_gr])) + 
  facet_wrap(~ n) + # Add facet_wrap for separate plots by "n"
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))


library("gridExtra")
library("grid") 

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, nullGrob(), p2, ncol = 3, widths = c(2, 0.1, 2.83))  # `nullGrob()` is an empty plot

# Figure 3:

ggsave("../figures/Fig3.eps", grid_plot, width=14, height=5.5)









# Figure 4: Mean AUC values per considered sample size and method for C = 10.
#############################################################################

resultstemp <- resultsall[resultsall$K==10,]
resultstemp$K <- NULL


# Compute AUC for specified groups and aggregate results
resultsK10_vstwo_gr <- resultstemp %>%
  filter(!(method_all %in% c("bin_wgini_wsquared", "bin_wgini_wosquared", 
                             "bin_wogini_wsquared", "bin_wogini_wosquared"))) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(51:53, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    mean_auc_4 = mean(auc_4),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3, mean_auc_4) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3", "cl_as_4"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK10_vstwo_gr %>%
  pivot_longer(cols = starts_with("cl_as_"), names_to = "cl_as", values_to = "value") %>%
  mutate(type = case_when(
    method_all %in% c("perm", "gini_corr")  ~ "conventional",
    str_detect(method_all, "muw_m_") ~ "multi-class diff.",
    str_detect(method_all, "muw_") ~ "multi-class",
    TRUE ~ NA_character_
  ))

results_long <- results_long %>%
  mutate(method = case_when(
    method_all %in% c("muw_wgini_wsquared", "muw_m_bin_wgini_wsquared") ~ "wsquared_wgini",
    method_all %in% c("muw_wgini_wosquared", "muw_m_bin_wgini_wosquared") ~ "wosquared_wgini",
    method_all %in% c("muw_wogini_wsquared", "muw_m_bin_wogini_wsquared") ~ "wsquared_wogini",
    method_all %in% c("muw_wogini_wosquared", "muw_m_bin_wogini_wosquared") ~ "wosquared_wogini",
    method_all == "perm" ~ "perm",
    method_all == "gini_corr" ~ "gini_corr"
  ),
  method = factor(method, levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", 
                                   "wosquared_wogini", "perm", "gini_corr")))

p1 <- ggplot(results_long, aes(x = cl_as, y = value, group = method_all, color = method, linetype = type)) +
  geom_line() +
  geom_point(aes(shape = method), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=16),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)
  ) +
  ylab("AUC") +
  scale_linetype_manual(values = c("conventional" = "solid", "multi-class diff." = 
                                     "dashed", "multi-class" = "dotdash")) +
  scale_shape_manual(values=c("wsquared_wgini" = 1, "wosquared_wgini" = 2, "wsquared_wogini" = 3, 
                              "wosquared_wogini" = 4, "perm" = 5, "gini_corr" = 6)) +
  scale_x_discrete(labels = c("cl_as_1" = expression(X[cl_as[1]]), "cl_as_2" = 
                                expression(X[cl_as[2]]), "cl_as_3" = expression(X[cl_as[3]]), 
                              "cl_as_4" = expression(X[cl_as[4]]))) +
  theme(legend.position = "none") +
  ggtitle(expression("Comparison with" ~ X[two_gr])) + 
  facet_wrap(~ n)  # Add facet_wrap for separate plots by "n"




# Compute AUC for specified groups and aggregate results
resultsK10_vsthr_gr <- resultstemp %>%
  filter(!(method_all %in% c("bin_wgini_wsquared", "bin_wgini_wosquared", "bin_wogini_wsquared", "bin_wogini_wosquared"))) %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(54:56, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    mean_auc_4 = mean(auc_4),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3, mean_auc_4) %>%
  rename_with(
    ~ c("cl_as_1", "cl_as_2", "cl_as_3", "cl_as_4"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK10_vsthr_gr %>%
  pivot_longer(cols = starts_with("cl_as_"), names_to = "cl_as", values_to = "value") %>%
  mutate(type = case_when(
    method_all %in% c("perm", "gini_corr")  ~ "conventional",
    str_detect(method_all, "muw_m_") ~ "multi-class diff.",
    str_detect(method_all, "muw_") ~ "multi-class",
    TRUE ~ NA_character_
  ))

results_long <- results_long %>%
  mutate(method = case_when(
    method_all %in% c("muw_wgini_wsquared", "muw_m_bin_wgini_wsquared") ~ "wsquared_wgini",
    method_all %in% c("muw_wgini_wosquared", "muw_m_bin_wgini_wosquared") ~ "wosquared_wgini",
    method_all %in% c("muw_wogini_wsquared", "muw_m_bin_wogini_wsquared") ~ "wsquared_wogini",
    method_all %in% c("muw_wogini_wosquared", "muw_m_bin_wogini_wosquared") ~ "wosquared_wogini",
    method_all == "perm" ~ "perm",
    method_all == "gini_corr" ~ "gini_corr"
  ),
  method = factor(method, levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", 
                                   "wosquared_wogini", "perm", "gini_corr")))

p2 <- ggplot(results_long, aes(x = cl_as, y = value, group = method_all, color = method, linetype = type)) +
  geom_line() +
  geom_point(aes(shape = method), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=16),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16)
  ) +
  labs(color="Method", shape="Method", linetype="Type", y="AUC") +
  scale_linetype_manual(values = c("conventional" = "solid", "multi-class diff." = "dashed", 
                                   "multi-class" = "dotdash")) +
  scale_shape_manual(values=c("wsquared_wgini" = 1, "wosquared_wgini" = 2, "wsquared_wogini" = 3, 
                              "wosquared_wogini" = 4, "perm" = 5, "gini_corr" = 6)) +
  scale_x_discrete(labels = c("cl_as_1" = expression(X[cl_as[1]]), "cl_as_2" = expression(X[cl_as[2]]), 
                              "cl_as_3" = expression(X[cl_as[3]]), "cl_as_4" = expression(X[cl_as[4]]))) +
  theme(legend.position = "right") +
  ggtitle(expression("Comparison with" ~ X[thr_gr])) + 
  facet_wrap(~ n) + # Add facet_wrap for separate plots by "n"
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1),
         linetype  = guide_legend(order = 2))

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, nullGrob(), p2, ncol = 3, widths = c(2, 0.1, 2.83))  # `nullGrob()` is an empty plot

# Figure 4:

ggsave("../figures/Fig4.eps", grid_plot, width=15, height=5.5)











# Figure S2: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 4 and n = 100.
######################################################################################################

# Determine two colors that are well-distinguishable when printing in 
# black and white:

library("RColorBrewer")

# display.brewer.pal(n = 9, name = "YlGnBu")
colors <- brewer.pal(9, "YlGnBu")[c(2, 5)]



Ktemp <- 4

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 100", 
                               !(method_all %in% c("perm", "gini_corr", 
                                                   "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared", 
                                                   "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=5), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)

resi2 <- resi2 %>%
  mutate(label_y = if_else(method != "wosquared_wogini", label_y / 2, label_y))

p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5)) +
  theme(legend.position = c(0.075, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S2:

ggsave("../figures/FigS2.pdf", width = 12, height = 6)





# Figure S3: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 4 and n = 500.
######################################################################################################

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 500", 
                               !(method_all %in% c("perm", "gini_corr", "muw_m_bin_wgini_wsquared",
                                                   "muw_m_bin_wgini_wosquared", "muw_m_bin_wogini_wsquared", 
                                                   "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=5), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5)) +
  theme(legend.position = c(0.075, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S3:

ggsave("../figures/FigS3.pdf", width = 12, height = 6)





# Figure S4: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 4 and n = 1000.
######################################################################################################

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 1000", 
                               !(method_all %in% c("perm", "gini_corr", "muw_m_bin_wgini_wsquared", 
                                                   "muw_m_bin_wgini_wosquared", "muw_m_bin_wogini_wsquared", 
                                                   "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=5), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (6 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5)) +
  theme(legend.position = c(0.075, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S4:

ggsave("../figures/FigS4.pdf", width = 12, height = 6)





# Figure S5: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 4 and n = 2000.
######################################################################################################

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 2000", !(method_all %in% c("perm", "gini_corr",
                                                                  "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared",
                                                                  "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=5), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (6 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5)) +
  theme(legend.position = c(0.075, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S5:

ggsave("../figures/FigS5.pdf", width = 12, height = 6)










# Figure S6: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 6 and n = 100.
######################################################################################################


Ktemp <- 6

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 100", 
                               !(method_all %in% c("perm", "gini_corr", 
                                                   "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared", 
                                                   "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=6), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)

resi2 <- resi2 %>%
  mutate(label_y = if_else(method != "wosquared_wogini", label_y / 2, label_y))


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5)) +
  theme(legend.position = c(0.061, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S6:

ggsave("../figures/FigS6.pdf", width = 12, height = 6)







# Figure S7: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 6 and n = 500.
######################################################################################################


resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 500", !(method_all %in% c("perm", "gini_corr", 
                                                                 "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared", 
                                                                 "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=6), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)

resi2 <- resi2 %>% mutate(label_y = if_else(method=="wsquared_wgini", label_y/2, label_y))


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5)) +
  theme(legend.position = c(0.061, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S7:

ggsave("../figures/FigS7.pdf", width = 12, height = 6)






# Figure S8: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 6 and n = 1000.
######################################################################################################


resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 1000", 
                               !(method_all %in% c("perm", "gini_corr", 
                                                   "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared", 
                                                   "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=6), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5)) +
  theme(legend.position = c(0.061, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S8:

ggsave("../figures/FigS8.pdf", width = 12, height = 6)





# Figure S9: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 6 and n = 2000.
######################################################################################################


resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 2000", 
                               !(method_all %in% c("perm", "gini_corr", 
                                                   "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared", 
                                                   "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=6), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5)) +
  theme(legend.position = c(0.061, 0.9), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S9:

ggsave("../figures/FigS9.pdf", width = 12, height = 6)










# Figure S10: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 10 and n = 100.
######################################################################################################


Ktemp <- 10

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 100", 
                               !(method_all %in% c("perm", "gini_corr", 
                                                   "muw_m_bin_wgini_wsquared", "muw_m_bin_wgini_wosquared", 
                                                   "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", 
                        "wosquared_wogini"), each=7), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19, 22), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", 
                 "X[cl_as[3]]", "X[cl_as[4]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)

resi2 <- resi2 %>%
  mutate(label_y = if_else(method != "wosquared_wogini", label_y / 2, label_y))


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5, 20.5)) +
  theme(legend.position = c(0.053, 0.9), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S10:

ggsave("../figures/FigS10.pdf", width = 12, height = 6)






# Figure S11: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 10 and n = 500.
######################################################################################################

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 500", !(method_all %in% c("perm", "gini_corr", "muw_m_bin_wgini_wsquared", 
                                                                 "muw_m_bin_wgini_wosquared", 
                                                                 "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=7), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19, 22), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]", "X[cl_as[4]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5, 20.5)) +
  theme(legend.position = c(0.053, 0.9), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S11:

ggsave("../figures/FigS11.pdf", width = 12, height = 6)





# Figure S12: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 10 and n = 1000.
######################################################################################################

resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 1000", !(method_all %in% c("perm", "gini_corr", "muw_m_bin_wgini_wsquared",
                                                                  "muw_m_bin_wgini_wosquared", 
                                                                  "muw_m_bin_wogini_wsquared", "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), each=7), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19, 22), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]", "X[cl_as[4]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5, 20.5)) +
  theme(legend.position = c(0.053, 0.9), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S12:

ggsave("../figures/FigS12.pdf", width = 12, height = 6)





# Figure S13: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets with C = 10 and n = 2000.
######################################################################################################


resultstemp <- resultsall %>%
  filter(K == Ktemp) %>%
  select(-K)

resi <- resultstemp %>% filter(n=="n = 2000", 
                               !(method_all %in% c("perm", "gini_corr", "muw_m_bin_wgini_wsquared", 
                                                   "muw_m_bin_wgini_wosquared", "muw_m_bin_wogini_wsquared", 
                                                   "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% c(1:45)))



resi2 <- data.frame(
  method = factor(rep(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", 
                        "wosquared_wogini"), each=7), 
                  levels=c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", 
                           "wosquared_wogini")),
  xpos = rep(c(2.75, 7, 10, 13, 16, 19, 22), times=4),
  labels = rep(c("X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]", "X[cl_as[4]]"), times=4)
)

# Calculate the minimum y-value for each method group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(method) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "method") %>%
  select(method, xpos, labels, label_y)


p <- ggplot(resi, aes(x = factor(seq), y = vim)) + facet_wrap(~ method, scales="free_y") +
  geom_boxplot(aes(fill = type), position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill="VIM type",
       x = "Covariates",
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  geom_vline(xintercept = c(5.5, 8.5, 11.5, 14.5, 17.5, 20.5)) +
  theme(legend.position = c(0.053, 0.9), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.box.background = element_rect(fill = NA, colour = NA)) + 
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE)

p

# Figure S13:

ggsave("../figures/FigS13.pdf", width = 12, height = 6)











# Figure S14: Multi-class and discriminatory VIM values obtained for the four versions of multi-forest
# for all simulated datasets: noise covariates only.
#######################################################################################################

resi <- resultsall %>% filter(n=="n = 500", 
                              !(method_all %in% c("perm", "gini_corr", "muw_m_bin_wgini_wsquared", 
                                                  "muw_m_bin_wgini_wosquared", "muw_m_bin_wogini_wsquared", 
                                                  "muw_m_bin_wogini_wosquared")))

resi <- resi %>% mutate(method = case_when(method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared") ~ "wsquared_wgini",
                                           method_all %in% c("muw_wgini_wosquared", "bin_wgini_wosquared") ~ "wosquared_wgini",
                                           method_all %in% c("muw_wogini_wsquared", "bin_wogini_wsquared") ~ "wsquared_wogini",
                                           method_all %in% c("muw_wogini_wosquared", "bin_wogini_wosquared") ~ "wosquared_wogini"),
                        type = case_when(str_detect(method_all, "muw_") ~ "multi-class",
                                         str_detect(method_all, "bin_") ~ "discriminatory"),
                        method = factor(method, levels = c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")))

resi <- resi %>%
  group_by(n, K, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(seq %in% 1:5)

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

p <- ggplot(resi, aes(x = factor(seq), y = vim, fill = type)) + facet_wrap(~ K + method, scales="free_y") +
  geom_boxplot(position = position_dodge(width = 0.75)) + 
  theme_bw() +
  labs(fill = "VIM type",
       x = expression(X[no]),
       y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.x = element_text(size = 14),
        axis.title=element_text(size=14))

p

# Figure S14:

ggsave("../figures/FigS14.pdf", width = 12, height = 7)











# Figure 1: VIM values obtained for wsquared wgini and the permutation VIM (perm) obtained for all
# simulated datasets with n = 500.
#################################################################################################

resi <- resultsall %>% filter(n=="n = 500", method_all %in% c("muw_wgini_wsquared", "bin_wgini_wsquared"))

resi <- resi %>% mutate(method = factor(method_all, levels = c("bin_wgini_wsquared", "muw_wgini_wsquared")))

resi <- resi %>%
  group_by(n, K, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]", "X[cl_as[4]]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p1 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(aes(fill = method), position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "wsquared_wgini", fill="VIM type", x = "Covariates", y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("discriminatory", "multi-class")) +
  theme(legend.position = c(0.145, 0.94), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16),
        legend.background = element_rect(fill = NA, colour = NA),  # Transparent legend background and no border
        legend.key = element_rect(fill = NA, colour = NA),  # Transparent keys and no border
        legend.box.background = element_rect(fill = NA, colour = NA)) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)



resi <- resultsall %>% filter(n=="n = 500", method_all %in% c("perm"))

resi <- resi %>%
  group_by(n, K, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), 
             levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_as[1]]", "X[cl_as[2]]", "X[cl_as[3]]", "X[cl_as[4]]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p2 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "perm", x = "Covariates", y = "VIM values") +
  theme(legend.position = "none", axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, p2, ncol = 2)  # `nullGrob()` is an empty plot

# Figure 1:

ggsave("../figures/Fig1.eps", grid_plot, width = 12, height = 10)















# Figure S1: Class-specific distributions of the informative covariates.
######################################################################

color_palette <- brewer.pal(4, "Set1")
# color_palette <- colorRampPalette(set1_original)(10)

xseq <- seq(-4, 7, length.out = 200)

dtemp <- 0.015


plotdata <- data.frame(c=factor(rep(1:4, each=200)),
                       x=rep(xseq, 4),
                       y=c(dnorm(xseq, mean=0, sd=1), 
                           dnorm(xseq, mean=0, sd=1) + dtemp,
                           dnorm(xseq, mean=1.5, sd=1), 
                           dnorm(xseq, mean=1.5, sd=1) + dtemp))

p1 <- ggplot(data=plotdata, aes(x=x, y=y, color=c)) + theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[two_gr]))



plotdata <- data.frame(
  c=factor(rep(1:4, each=200)),
  x=rep(xseq, 4),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=1, sd=1))
)


p2 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_1]))





plotdata <- data.frame(
  c=factor(rep(1:4, each=200)),
  x=rep(xseq, 4),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=1, sd=1),
      dnorm(xseq, mean=2, sd=1))
)


p3 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_2]))





plotdata <- data.frame(
  c=factor(rep(1:4, each=200)),
  x=rep(xseq, 4),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=2.25, sd=1))
)


p4 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_3]))

# Arrange plots in a grid
C4 <- grid.arrange(p1, p2, p3, p4, ncol=3)

# Create a text grob for the title with increased font size
title_grob <- textGrob("C = 4", gp=gpar(fontsize=20))

# Use this grob as the top parameter in a new grid.arrange call
C4 <- grid.arrange(grobs = list(C4), top = title_grob)




color_palette <- brewer.pal(6, "Set1")

xseq <- seq(-4, 7, length.out = 200)


plotdata <- data.frame(c=factor(rep(1:6, each=200)),
                       x=rep(xseq, 6),
                       y=c(dnorm(xseq, mean=0, sd=1), 
                           dnorm(xseq, mean=0, sd=1) + dtemp,
                           dnorm(xseq, mean=0, sd=1) + 2*dtemp,
                           dnorm(xseq, mean=1.5, sd=1), 
                           dnorm(xseq, mean=1.5, sd=1) + dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 2*dtemp))

p1 <- ggplot(data=plotdata, aes(x=x, y=y, color=c)) + theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[two_gr]))



plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=1, sd=1), 
      dnorm(xseq, mean=1, sd=1) + dtemp,
      dnorm(xseq, mean=2, sd=1), 
      dnorm(xseq, mean=2, sd=1) + dtemp)
)


p2 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[thr_gr]))




plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=1, sd=1))
)


p3 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_1]))



plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=1, sd=1),
      dnorm(xseq, mean=2, sd=1))
)


p4 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_2]))





plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=2.25, sd=1))
)


p5 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_3]))

# Arrange plots in a grid
C6 <- grid.arrange(p1, p2, p3, p4, p5, ncol=3)

# Create a text grob for the title with increased font size
title_grob <- textGrob("C = 6", gp=gpar(fontsize=20))

# Use this grob as the top parameter in a new grid.arrange call
C6 <- grid.arrange(grobs = list(C6), top = title_grob)




set1_original <- brewer.pal(8, "Set1")
color_palette <- colorRampPalette(set1_original)(10)

xseq <- seq(-4, 7, length.out = 200)


plotdata <- data.frame(c=factor(rep(1:10, each=200)),
                       x=rep(xseq, 10),
                       y=c(dnorm(xseq, mean=0, sd=1), 
                           dnorm(xseq, mean=0, sd=1) + dtemp,
                           dnorm(xseq, mean=0, sd=1) + 2*dtemp,
                           dnorm(xseq, mean=0, sd=1) + 3*dtemp,
                           dnorm(xseq, mean=0, sd=1) + 4*dtemp,
                           dnorm(xseq, mean=1.5, sd=1), 
                           dnorm(xseq, mean=1.5, sd=1) + dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 2*dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 3*dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 4*dtemp))

p1 <- ggplot(data=plotdata, aes(x=x, y=y, color=c)) + theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[two_gr]))



plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=1, sd=1), 
      dnorm(xseq, mean=1, sd=1) + dtemp,
      dnorm(xseq, mean=1, sd=1) + 2*dtemp,
      dnorm(xseq, mean=2, sd=1), 
      dnorm(xseq, mean=2, sd=1) + dtemp,
      dnorm(xseq, mean=2, sd=1) + 2*dtemp)
)


p2 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[thr_gr]))




plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=0, sd=1) + 5*dtemp,
      dnorm(xseq, mean=0, sd=1) + 6*dtemp,
      dnorm(xseq, mean=0, sd=1) + 7*dtemp,
      dnorm(xseq, mean=0, sd=1) + 8*dtemp,
      dnorm(xseq, mean=1, sd=1))
)


p3 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_1]))



plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=0, sd=1) + 5*dtemp,
      dnorm(xseq, mean=0, sd=1) + 6*dtemp,
      dnorm(xseq, mean=0, sd=1) + 7*dtemp,
      dnorm(xseq, mean=1, sd=1),
      dnorm(xseq, mean=2, sd=1))
)


p4 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_2]))




plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=0, sd=1) + 5*dtemp,
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=0.75, sd=1) + dtemp,
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=2.25, sd=1))
)


p5 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_3]))





plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=0.75, sd=1) + dtemp,
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=1.5, sd=1) + dtemp,
      dnorm(xseq, mean=2.25, sd=1),
      dnorm(xseq, mean=3, sd=1))
)

p6 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_as_4]))

# Arrange plots in a grid
C10 <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)

# Create a text grob for the title with increased font size
title_grob <- textGrob("C = 10", gp=gpar(fontsize=20))

# Use this grob as the top parameter in a new grid.arrange call
C10 <- grid.arrange(grobs = list(C10), top = title_grob)




# Create a spacer plot
spacer <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank())


library("cowplot")

combined_plot <- plot_grid(
  C4,
  spacer,
  C6,
  spacer,
  C10,
  ncol = 1,
  rel_heights = c(1, 0.1, 1, 0.1, 1)
)


# Figure S1:

ggsave("../figures/FigS1.pdf", width=10*0.8, height=13*0.8)
