####################################################################################

# NOTE: Before the code can be executed, the R working directory *MUST* 
# be set to the directory 'MultiForests_code_and_data/benchmark_study' (that is
# the directory in which this R script is contained):

# Remove '#' from the line below and replace 'here/is/my/path/' by the path
# to the directory 'MultiForests_code_and_data/benchmark_study':

# setwd("here/is/my/path/")

####################################################################################

# Load and pre-process the results:
#####################################

load("./intermediate_results/scenariogrid_benchmark_study.Rda")
load("./intermediate_results/results_benchmark_study.Rda")
load("./data/datainfo.Rda")

results <- scenariogrid
results$acc <- sapply(Results, function(x) mean(sapply(x, function(y) y$acctest)))
results$brier <- sapply(Results, function(x) mean(sapply(x, function(y) y$briertest)))
results$aunp <- sapply(Results, function(x) mean(sapply(x, function(y) y$aunptest)))
results$aunu <- sapply(Results, function(x) mean(sapply(x, function(y) y$aunutest)))

results$n <- sapply(results$dataset, function(x) datainfo$n[datainfo$filename==x])
results$p <- sapply(results$dataset, function(x) datainfo$p[datainfo$filename==x])
results$n_cl <- sapply(results$dataset, function(x) datainfo$n_cl[datainfo$filename==x])
results$prop_cat <- sapply(results$dataset, function(x) datainfo$prop_cat[datainfo$filename==x])

results$method[results$method=="muwf_wgini_wsquared"] <- "MuwF_gini_squ"
results$method[results$method=="muwf_wgini_wosquared"] <- "MuwF_gini_WOsqu"
results$method[results$method=="muwf_wogini_wsquared"] <- "MuwF_WOgini_squ"
results$method[results$method=="muwf_wogini_wosquared"] <- "MuwF_WOgini_WOsqu"
results$method[results$method=="rf"] <- "RF"

results$method <- factor(results$method, levels=c("MuwF_gini_squ", "MuwF_gini_WOsqu", 
                                                  "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu", "RF"))

sort(tapply(results$acc, results$method, mean))
sort(tapply(results$acc, results$method, median))

sort(tapply(results$brier, results$method, mean), decreasing=TRUE)
sort(tapply(results$brier, results$method, median), decreasing=TRUE)

sort(tapply(results$aunp, results$method, mean))
sort(tapply(results$aunp, results$method, median))

sort(tapply(results$aunu, results$method, mean))
sort(tapply(results$aunu, results$method, median))



orderind <- order(results$dataset, results$cvind, results$method)
results <- results[orderind,]
Results <- Results[orderind]

results$seed <- NULL
results$settingid <- NULL


library("dplyr")

resultssum <- results %>% group_by(dataset, method) %>% summarise(n = first(n), p = first(p),
                                                                  n_cl = first(n_cl),
                                                                  prop_cat = first(prop_cat),
                                                                  acc = mean(acc),
                                                                  brier = mean(brier),
                                                                  aunp = mean(aunp),
                                                                  aunu = mean(aunu),
                                                                  .groups = 'drop')






# Table 2: Performances of the methods summarized across the 121 datasets.
##########################################################################

library("forcats")

# Calculate medians and quantiles for each metric and method
summary_df <- resultssum %>%
  group_by(method) %>%
  mutate(method = as.factor(method)) %>% 
  mutate(method = fct_relevel(method, "RF")) %>%
  summarise(
    acc_median = median(acc, na.rm = TRUE),
    acc_Q1 = quantile(acc, probs = 0.25, na.rm = TRUE),
    acc_Q3 = quantile(acc, probs = 0.75, na.rm = TRUE),
    brier_median = median(brier, na.rm = TRUE),
    brier_Q1 = quantile(brier, probs = 0.25, na.rm = TRUE),
    brier_Q3 = quantile(brier, probs = 0.75, na.rm = TRUE),
    aunp_median = median(aunp, na.rm = TRUE),
    aunp_Q1 = quantile(aunp, probs = 0.25, na.rm = TRUE),
    aunp_Q3 = quantile(aunp, probs = 0.75, na.rm = TRUE),
    aunu_median = median(aunu, na.rm = TRUE),
    aunu_Q1 = quantile(aunu, probs = 0.25, na.rm = TRUE),
    aunu_Q3 = quantile(aunu, probs = 0.75, na.rm = TRUE)
  ) %>%
  # Rename methods
  mutate(method = case_when(
    method == "MuwF_gini_squ" ~ "wsquared_wgini",
    method == "MuwF_gini_WOsqu" ~ "wosquared_wgini",
    method == "MuwF_WOgini_squ" ~ "wsquared_wogini",
    method == "MuwF_WOgini_WOsqu" ~ "wosquared_wogini",
    TRUE ~ method
  )) %>%
  # Format each metric with median [Q1, Q3]
  mutate(
    acc = paste(format(round(acc_median, 4), nsmall = 4),
                " [", format(round(acc_Q1, 4), nsmall = 4),
                ", ", format(round(acc_Q3, 4), nsmall = 4), "]", sep = ""),
    brier = paste(format(round(brier_median, 4), nsmall = 4),
                  " [", format(round(brier_Q1, 4), nsmall = 4),
                  ", ", format(round(brier_Q3, 4), nsmall = 4), "]", sep = ""),
    aunp = paste(format(round(aunp_median, 4), nsmall = 4),
                 " [", format(round(aunp_Q1, 4), nsmall = 4),
                 ", ", format(round(aunp_Q3, 4), nsmall = 4), "]", sep = ""),
    aunu = paste(format(round(aunu_median, 4), nsmall = 4),
                 " [", format(round(aunu_Q1, 4), nsmall = 4),
                 ", ", format(round(aunu_Q3, 4), nsmall = 4), "]", sep = "")
  ) %>%
  # Select only the columns needed for the final summary
  select(method, aunu, aunp, brier, acc) %>%
  rename(
    method = method,
    AUNU = aunu,
    AUNP = aunp,
    Brier = brier,
    ACC = acc
  )

# View the summary dataframe
print(summary_df)

library("xtable")

# Convert your dataframe to a LaTeX table with xtable
latex_table <- xtable(summary_df, caption = "Performances of the methods summarised across the 121 datasets", 
                      label = "tab:benchmark_summary")

# Table 2:

# Print the LaTeX table to the console or to a file
print(latex_table, type = "latex", file = "../tables/Tab2.tex", include.rownames = FALSE)








# Table 3: Testing for significant differences between the performance of IF and the other methods:
###################################################################################################

library("rstatix")

paunus <- c(wilcox.test(resultssum$aunu[resultssum$method=="RF"], resultssum$aunu[resultssum$method=="MuwF_gini_squ"], paired=TRUE)$p.value,
            wilcox.test(resultssum$aunu[resultssum$method=="RF"], resultssum$aunu[resultssum$method=="MuwF_gini_WOsqu"], paired=TRUE)$p.value,
            wilcox.test(resultssum$aunu[resultssum$method=="RF"], resultssum$aunu[resultssum$method=="MuwF_WOgini_squ"], paired=TRUE)$p.value,
            wilcox.test(resultssum$aunu[resultssum$method=="RF"], resultssum$aunu[resultssum$method=="MuwF_WOgini_WOsqu"], paired=TRUE)$p.value)

raunus <- rep(NA, 4)
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_squ"))
raunus[1] <- wilcox_effsize(data = datatemp, formula=aunu ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_WOsqu"))
raunus[2] <- wilcox_effsize(data = datatemp, formula=aunu ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_squ"))
raunus[3] <- wilcox_effsize(data = datatemp, formula=aunu ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_WOsqu"))
raunus[4] <- wilcox_effsize(data = datatemp, formula=aunu ~ method, paired = TRUE)$effsize
raunus <- round(raunus, 2)

paunusadjust <- p.adjust(paunus, method = "holm")
names(paunusadjust) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")
names(raunus) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")


paunps <- c(wilcox.test(resultssum$aunp[resultssum$method=="RF"], resultssum$aunp[resultssum$method=="MuwF_gini_squ"], paired=TRUE)$p.value,
            wilcox.test(resultssum$aunp[resultssum$method=="RF"], resultssum$aunp[resultssum$method=="MuwF_gini_WOsqu"], paired=TRUE)$p.value,
            wilcox.test(resultssum$aunp[resultssum$method=="RF"], resultssum$aunp[resultssum$method=="MuwF_WOgini_squ"], paired=TRUE)$p.value,
            wilcox.test(resultssum$aunp[resultssum$method=="RF"], resultssum$aunp[resultssum$method=="MuwF_WOgini_WOsqu"], paired=TRUE)$p.value)

raunps <- rep(NA, 4)
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_squ"))
raunps[1] <- wilcox_effsize(data = datatemp, formula=aunp ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_WOsqu"))
raunps[2] <- wilcox_effsize(data = datatemp, formula=aunp ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_squ"))
raunps[3] <- wilcox_effsize(data = datatemp, formula=aunp ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_WOsqu"))
raunps[4] <- wilcox_effsize(data = datatemp, formula=aunp ~ method, paired = TRUE)$effsize
raunps <- round(raunps, 2)

paunpsadjust <- p.adjust(paunps, method = "holm")
names(paunpsadjust) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")
names(raunps) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")



paccs <- c(wilcox.test(resultssum$acc[resultssum$method=="RF"], resultssum$acc[resultssum$method=="MuwF_gini_squ"], paired=TRUE)$p.value,
           wilcox.test(resultssum$acc[resultssum$method=="RF"], resultssum$acc[resultssum$method=="MuwF_gini_WOsqu"], paired=TRUE)$p.value,
           wilcox.test(resultssum$acc[resultssum$method=="RF"], resultssum$acc[resultssum$method=="MuwF_WOgini_squ"], paired=TRUE)$p.value,
           wilcox.test(resultssum$acc[resultssum$method=="RF"], resultssum$acc[resultssum$method=="MuwF_WOgini_WOsqu"], paired=TRUE)$p.value)

raccs <- rep(NA, 4)
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_squ"))
raccs[1] <- wilcox_effsize(data = datatemp, formula=acc ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_WOsqu"))
raccs[2] <- wilcox_effsize(data = datatemp, formula=acc ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_squ"))
raccs[3] <- wilcox_effsize(data = datatemp, formula=acc ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_WOsqu"))
raccs[4] <- wilcox_effsize(data = datatemp, formula=acc ~ method, paired = TRUE)$effsize
raccs <- round(raccs, 2)

paccsadjust <- p.adjust(paccs, method = "holm")
names(paccsadjust) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")
names(raccs) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")




pbriers <- c(wilcox.test(resultssum$brier[resultssum$method=="RF"], resultssum$brier[resultssum$method=="MuwF_gini_squ"], paired=TRUE)$p.value,
             wilcox.test(resultssum$brier[resultssum$method=="RF"], resultssum$brier[resultssum$method=="MuwF_gini_WOsqu"], paired=TRUE)$p.value,
             wilcox.test(resultssum$brier[resultssum$method=="RF"], resultssum$brier[resultssum$method=="MuwF_WOgini_squ"], paired=TRUE)$p.value,
             wilcox.test(resultssum$brier[resultssum$method=="RF"], resultssum$brier[resultssum$method=="MuwF_WOgini_WOsqu"], paired=TRUE)$p.value)

rbriers <- rep(NA, 4)
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_squ"))
rbriers[1] <- wilcox_effsize(data = datatemp, formula=brier ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_gini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_gini_WOsqu"))
rbriers[2] <- wilcox_effsize(data = datatemp, formula=brier ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_squ"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_squ"))
rbriers[3] <- wilcox_effsize(data = datatemp, formula=brier ~ method, paired = TRUE)$effsize
datatemp <- resultssum[resultssum$method %in% c("RF", "MuwF_WOgini_WOsqu"),]
datatemp$method <- factor(as.character(datatemp$method), levels=c("RF", "MuwF_WOgini_WOsqu"))
rbriers[4] <- wilcox_effsize(data = datatemp, formula=brier ~ method, paired = TRUE)$effsize
rbriers <- round(rbriers, 2)

pbriersadjust <- p.adjust(pbriers, method = "holm")
names(pbriersadjust) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")
names(rbriers) <- c("MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu")


paunusadjust <- ifelse(paunusadjust < 0.001, "< 0.001", format(round(paunusadjust, 3), nsmall=3))
paunpsadjust <- ifelse(paunpsadjust < 0.001, "< 0.001", format(round(paunpsadjust, 3), nsmall=3))
paccsadjust <- ifelse(paccsadjust < 0.001, "< 0.001", format(round(paccsadjust, 3), nsmall=3))
pbriersadjust <- ifelse(pbriersadjust < 0.001, "< 0.001", format(round(pbriersadjust, 3), nsmall=3))

raunus <- format(raunus, nsmall=2)
raunps <- format(raunps, nsmall=2)
raccs <- format(raccs, nsmall=2)
rbriers <- format(rbriers, nsmall=2)



# p-values:
cat(paste(paste(names(paunusadjust), " vs.\\ RF: p =", paunusadjust, sep=""), collapse=", "), "\n")
cat(paste(paste(names(paunpsadjust), " vs.\\ RF: p =", paunpsadjust, sep=""), collapse=", "), "\n")
cat(paste(paste(names(paccsadjust), " vs.\\ RF: p =", paccsadjust, sep=""), collapse=", "), "\n")
cat(paste(paste(names(pbriersadjust), " vs.\\ RF: p =", pbriersadjust, sep=""), collapse=", "), "\n")

# Effect sizes:
cat(paste(paste(names(raunus), " vs.\\ RF: r =", raunus, sep=""), collapse=", "), "\n")
cat(paste(paste(names(raunps), " vs.\\ RF: r =", raunps, sep=""), collapse=", "), "\n")
cat(paste(paste(names(raccs), " vs.\\ RF: r =", raccs, sep=""), collapse=", "), "\n")
cat(paste(paste(names(rbriers), " vs.\\ RF: r =", rbriers, sep=""), collapse=", "), "\n")

restab <- cbind(c("wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini"), paste0(raunus, " (", paunusadjust, ")"),
                paste0(raunps, " (", paunpsadjust, ")"),
                paste0(raccs, " (", paccsadjust, ")"),
                paste0(rbriers, " (", pbriersadjust, ")"))
colnames(restab) <- c("method", "AUNU", "AUNP", "ACC", "Brier")


# Convert your dataframe to a LaTeX table with xtable
latex_table <- xtable(restab, 
                      caption = "Results of Wilcoxon tests between the values of the cross-validated metrics obtained with the different MuF versions and conventional RF", 
                      label = "tab:wilcox_test")


# Table 3:

# Print the LaTeX table to the console or to a file
print(latex_table, type = "latex", file = "../tables/Tab3.tex", include.rownames = FALSE)








# Figure 5: Ranks of the methods with respect to the different performance metrics.
#################################################################################

aunuranks <- t(apply(-cbind(resultssum$aunu[resultssum$method=="MuwF_gini_squ"], resultssum$aunu[resultssum$method=="MuwF_gini_WOsqu"], 
                            resultssum$aunu[resultssum$method=="MuwF_WOgini_squ"], resultssum$aunu[resultssum$method=="MuwF_WOgini_WOsqu"],
                            resultssum$aunu[resultssum$method=="RF"]), 1, rank))
aunpranks <- t(apply(-cbind(resultssum$aunp[resultssum$method=="MuwF_gini_squ"], resultssum$aunp[resultssum$method=="MuwF_gini_WOsqu"], 
                            resultssum$aunp[resultssum$method=="MuwF_WOgini_squ"], resultssum$aunp[resultssum$method=="MuwF_WOgini_WOsqu"],
                            resultssum$aunp[resultssum$method=="RF"]), 1, rank))
accranks <- t(apply(-cbind(resultssum$acc[resultssum$method=="MuwF_gini_squ"], resultssum$acc[resultssum$method=="MuwF_gini_WOsqu"], 
                           resultssum$acc[resultssum$method=="MuwF_WOgini_squ"], resultssum$acc[resultssum$method=="MuwF_WOgini_WOsqu"],
                           resultssum$acc[resultssum$method=="RF"]), 1, rank))
brierranks <- t(apply(cbind(resultssum$brier[resultssum$method=="MuwF_gini_squ"], resultssum$brier[resultssum$method=="MuwF_gini_WOsqu"], 
                            resultssum$brier[resultssum$method=="MuwF_WOgini_squ"], resultssum$brier[resultssum$method=="MuwF_WOgini_WOsqu"],
                            resultssum$brier[resultssum$method=="RF"]), 1, rank))

resultssum$aunuranks <- as.vector(t(aunuranks))
resultssum$aunpranks <- as.vector(t(aunpranks))
resultssum$accranks <- as.vector(t(accranks))
resultssum$brierranks <- as.vector(t(brierranks))

ggdata <- reshape(as.data.frame(resultssum), varying=c("aunuranks", "aunpranks", "accranks", "brierranks"), 
                  v.names="rank", 
                  timevar="metric", times=c("aunuranks", "aunpranks", "accranks", "brierranks"),
                  direction="long")

ggdata$metric <- recode(ggdata$metric,
                        aunuranks = "AUNU",
                        aunpranks = "AUNP",
                        accranks = "ACC",
                        brierranks = "Brier")

ggdata$metric <- factor(ggdata$metric, levels=c("AUNU", "AUNP", "ACC", "Brier"))

ggdata$rankold <- ggdata$rank
ggdata$rank[ggdata$rankold <= 1.5] <- "1, 1.5"
ggdata$rank[ggdata$rankold >= 2 & ggdata$rankold <= 2.5] <- "2, 2.5"
ggdata$rank[ggdata$rankold >= 3 & ggdata$rankold <= 3.5] <- "3, 3.5"
ggdata$rank[ggdata$rankold >= 4 & ggdata$rankold <= 4.5] <- "4, 4.5"
ggdata$rank[ggdata$rankold == 5] <- "5"
ggdata$rank <- factor(ggdata$rank, levels=c("1, 1.5", "2, 2.5", "3, 3.5", "4, 4.5", "5"))

ggdata$method <- factor(ggdata$method, levels=c("RF", "MuwF_gini_squ", "MuwF_gini_WOsqu", "MuwF_WOgini_squ", "MuwF_WOgini_WOsqu"))
levels(ggdata$method) <- c("RF", "wsquared_wgini", "wosquared_wgini", "wsquared_wogini", "wosquared_wogini")

library("ggplot2")

p <- ggplot(data=ggdata, aes(fill=rank, x=method)) + theme_bw() + 
  theme(strip.text.x=element_text(size=11), axis.text.x = element_text(colour = "black", size=11, angle = 45, hjust = 1),
        axis.title=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x="Method", y="Number of datasets")  +
  geom_bar(position="stack", colour="black") + facet_wrap(~metric) + scale_fill_brewer(name="Rank", palette="Greens")
p

# Figure 5:

ggsave(file="../figures/Fig5.eps", width=8, height=8)










# Figures 6 and S16: Dataset-specific cross-validated performance metric values: 
# wsquared_wgini versus RF. 
#################################################################################

library("tidyr")

comparison_data <- resultssum %>%
  filter(method %in% c("MuwF_gini_squ", "RF")) %>%
  select(dataset, n, p, n_cl, prop_cat, method, brier, acc, aunp, aunu) %>%
  pivot_wider(names_from = method, values_from = c(brier, acc, aunp, aunu))

plot_data <- comparison_data %>% mutate(many_class = n_cl >= 10)


library("scales")

# "Transformed using negative complementary square root transformation"
custom_trans <- trans_new(
  name = "proportion",
  transform = function(x) -sqrt(1-x),
  inverse = function(x) 1 - x^2
)

p1 <- ggplot(plot_data, aes(x = aunu_MuwF_gini_squ, y = aunu_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60"),
                     labels = c(expression(C < 10), expression(C >= 10))) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19),
                     labels = c(expression(C < 10), expression(C >= 10))) +
  scale_x_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  scale_y_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  labs(title = "AUNU",
       color = "Number of classes",
       shape = "Number of classes",
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position=c(0.3, 0.82), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))


p1.2 <- ggplot(plot_data, aes(x = aunp_MuwF_gini_squ, y = aunp_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  scale_y_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  labs(title = "AUNP", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))




p2 <- ggplot(plot_data, aes(x = acc_MuwF_gini_squ, y = acc_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.2, 0.5, 0.8, 0.95, 0.99, 1)) +
  scale_y_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.2, 0.5, 0.8, 0.95, 0.99, 1)) +
  labs(title = "ACC", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))




# "Transformed using Square Root Scale"
prop_trans <- trans_new(
  name = "proportion",
  transform = function(x) sqrt(x),
  inverse = function(x) x^2
)

p3 <- ggplot(plot_data, aes(x = brier_MuwF_gini_squ, y = brier_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = prop_trans, limits = c(NA, 1), breaks = c(0.01, 0.05, 0.15, 0.25, 0.5, 1)) +
  scale_y_continuous(trans = prop_trans, limits = c(NA, 1), breaks = c(0.01, 0.05, 0.15, 0.25, 0.5, 1)) +
  labs(title = "Brier", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))
p3


library("patchwork")
ps <- p1 | p2 | p3
ps

# Figure 6:

ggsave("../figures/Fig6.eps", width=15*0.8, height=5*0.8)



# Figure S16 (includes AUNP):

ps <- (p1 | p1.2) / (p2 | p3)
ps

ggsave("../figures/FigS16.eps", width=10, height=9)











# Figure S15: Dataset-specific performance metric values: wsquared_wgini 
# versus RF - distinguishing between datasets that contain only numerical 
# covariates and those that also include categorical covariates.
#########################################################################


comparison_data <- resultssum %>%
  filter(method %in% c("MuwF_gini_squ", "RF")) %>%
  select(dataset, n, p, n_cl, prop_cat, method, brier, acc, aunp, aunu) %>%
  pivot_wider(names_from = method, values_from = c(brier, acc, aunp, aunu))

plot_data <- comparison_data %>% mutate(many_class = n_cl >= 10) %>% filter(prop_cat == 0)



# "Transformed using negative complementary square root transformation"
custom_trans <- trans_new(
  name = "proportion",
  transform = function(x) -sqrt(1-x),
  inverse = function(x) 1 - x^2
)

p1 <- ggplot(plot_data, aes(x = aunu_MuwF_gini_squ, y = aunu_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60"),
                     labels = c(expression(C < 10), expression(C >= 10))) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19),
                     labels = c(expression(C < 10), expression(C >= 10))) +
  scale_x_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  scale_y_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  labs(title = "AUNU", 
       color = "Number of classes",
       shape = "Number of classes",
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position=c(0.3, 0.82), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=13), 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))


p2 <- ggplot(plot_data, aes(x = acc_MuwF_gini_squ, y = acc_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.2, 0.5, 0.8, 0.95, 0.99, 1)) +
  scale_y_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.2, 0.5, 0.8, 0.95, 0.99, 1)) +
  labs(title = "ACC", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))


# "Transformed using Square Root Scale"
prop_trans <- trans_new(
  name = "proportion",
  transform = function(x) sqrt(x),
  inverse = function(x) x^2
)

p3 <- ggplot(plot_data, aes(x = brier_MuwF_gini_squ, y = brier_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = prop_trans, limits = c(NA, 1), breaks = c(0.01, 0.05, 0.15, 0.25, 0.5, 1)) +
  scale_y_continuous(trans = prop_trans, limits = c(NA, 1), breaks = c(0.01, 0.05, 0.15, 0.25, 0.5, 1)) +
  labs(title = "Brier", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))




# First, filter for the methods of interest and then spread the data for easy comparison
comparison_data <- resultssum %>%
  filter(method %in% c("MuwF_gini_squ", "RF")) %>%
  select(dataset, n, p, n_cl, prop_cat, method, brier, acc, aunp, aunu) %>%
  pivot_wider(names_from = method, values_from = c(brier, acc, aunp, aunu))

plot_data <- comparison_data %>% mutate(many_class = n_cl >= 10) %>% filter(prop_cat > 0)



# "Transformed using negative complementary square root transformation"
custom_trans <- trans_new(
  name = "proportion",
  transform = function(x) -sqrt(1-x),
  inverse = function(x) 1 - x^2
)

p4 <- ggplot(plot_data, aes(x = aunu_MuwF_gini_squ, y = aunu_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  scale_y_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.5, 0.8, 0.95, 0.99, 1)) +
  labs(title = "AUNU", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))


p5 <- ggplot(plot_data, aes(x = acc_MuwF_gini_squ, y = acc_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.2, 0.5, 0.8, 0.95, 0.99, 1)) +
  scale_y_continuous(trans = custom_trans, limits = c(NA, 1), breaks = c(0.2, 0.5, 0.8, 0.95, 0.99, 1)) +
  labs(title = "ACC", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))


# "Transformed using Square Root Scale"
prop_trans <- trans_new(
  name = "proportion",
  transform = function(x) sqrt(x),
  inverse = function(x) x^2
)

p6 <- ggplot(plot_data, aes(x = brier_MuwF_gini_squ, y = brier_RF)) + theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(color=many_class, shape=many_class)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray60")) +
  scale_shape_manual(values = c("TRUE" = 2, "FALSE" = 19)) +
  scale_x_continuous(trans = prop_trans, limits = c(NA, 1), breaks = c(0.01, 0.05, 0.15, 0.25, 0.5, 1)) +
  scale_y_continuous(trans = prop_trans, limits = c(NA, 1), breaks = c(0.01, 0.05, 0.15, 0.25, 0.5, 1)) +
  labs(title = "Brier", 
       x = "wsquared_wgini", 
       y = "RF") + 
  theme(legend.position="none", 
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        plot.title = element_text(size = 14))



library("cowplot")

# Combine the first three plots side-by-side without labels
top_row <- plot_grid(p1, p2, p3, labels = NULL, ncol = 3)

# Combine the last three plots side-by-side without labels
bottom_row <- plot_grid(p4, p5, p6, labels = NULL, ncol = 3)


# Create titles with left alignment
top_title <- ggdraw() + draw_label("  Datasets with only numerical variables", fontface = 'bold', hjust = 0, x = 0)
bottom_title <- ggdraw() + draw_label("  Datasets with at least one categorical variable", fontface = 'bold', hjust = 0, x = 0)

# Create a spacer plot for increasing the spacing between the plot elements
spacer_top <- ggdraw() + theme_void()  # Create an empty plot
spacer_between <- ggdraw() + theme_void()  # Create an empty plot
spacer_bottom <- ggdraw() + theme_void()  # Create an empty plot


combined_plot <- plot_grid(
  top_title,
  spacer_top,
  top_row,
  spacer_between,
  bottom_title,
  spacer_bottom,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.1, 0.03, 1, 0.08, 0.1, 0.03, 1)
)

# Print the combined plot
combined_plot

# Figure S15:

ggsave("../figures/FigS15.eps", width=14*0.8, height=10*0.8)
