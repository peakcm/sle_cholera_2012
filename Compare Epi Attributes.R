#### Header ####
# Use signed-rank test to compare epidemic attributes

#### Load libraries ####
library(Hmisc)
library(dplyr)

#### Notes ####
# The cholera and ebola epidemics in chiefdoms differed with respect to cases, CAR, order of onset, presence of cases.
# Among chiefdoms that had BOTH cholera and ebola, there was no correlation between the duration of the outbreaks
# If we say a chiefdom with 0 cases had an outbreak of 0 duration, then the p-value inflates to 0.06376, so we conclude that the durations are similar, barely.

#### Load data ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")
# Note df_cumulative has lots of demographic data

#### Load Reff data ####
ebola_Reff_summary <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_Reff_summary.csv")
names(ebola_Reff_summary) <- c("CHCODE", "days_Reff_over_1_ebola", "total_cases_ebola", "total_cases_generated_ebola", "total_exported_cases_ebola", "total_internal_cases_generated_ebola", "total_imported_cases_ebola", "export_import_ratio_ebola")

cholera_Reff_summary <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/cholera_Reff_summary.csv")
names(cholera_Reff_summary) <- c("CHCODE", "days_Reff_over_1_cholera", "total_cases_cholera", "total_cases_generated_cholera", "total_exported_cases_cholera", "total_internal_cases_generated_cholera", "total_imported_cases_cholera", "export_import_ratio_cholera")

#### Function ####
fcn_test <- function(x, y, plot = TRUE){
  scatterplot = plot(x, y)
  ro = round(rcorr(x,y, type="spearman")$r[1,2], 4)

  if (plot){scatterplot}
  
  return(c(ro = ro))
}

#### Test function ####
fcn_test(x = seq(1,20), y = seq(1,20))
fcn_test(x = seq(1,20), y = seq(1,20) + 20*runif(20))
fcn_test(x = seq(1,20), y = rnorm(20, mean = 5))

#### Look at all pairs ####

df_cumulative_rcorr <- df_cumulative %>% select(CHCODE, Pop2014, Pop2012, total_ebola, cholera, total_ebola_onset, cholera_onset, total_ebola_duration, cholera_duration, Prop_Urban, Education, Prop_Improved_Water_Source, Prop_Improved_Toilet, SES_Score, Household_Size, Shared_Toilet, Wealth_Index) %>%
  mutate(total_ebola_CAR = total_ebola / Pop2014,
         cholera_CAR = cholera / Pop2012) %>% apply(2, as.numeric)

df_cumulative_rcorr <- data.frame(df_cumulative_rcorr) %>% 
  left_join(ebola_Reff_summary[,c("CHCODE", "days_Reff_over_1_ebola", "export_import_ratio_ebola")], by = "CHCODE") %>%
  left_join(cholera_Reff_summary[,c("CHCODE", "days_Reff_over_1_cholera", "export_import_ratio_cholera")], by = "CHCODE")

out <- rcorr(as.matrix(df_cumulative_rcorr), type = "spearman")
View(out$r)

# Adjust for multiple testing
num_disease_metrics <- 14
num_DHS_metrics <- 8
num_tests <- (num_disease_metrics*(num_disease_metrics-1) + num_disease_metrics*num_DHS_metrics)/2

out_sig <- out$r * (out$P < 0.05/num_tests)
out_sig[out_sig==0] <- NA
# for (i in 1:nrow(out_sig)){
#   for (j in 1:ncol(out_sig)){
#     if (i <= j){out_sig[i,j] <- NA}
#   }
# }
cat("The minimum strength of association detected", min(abs(out_sig), na.rm = TRUE))
cat("The adjusted p threshold is ", round(0.05/num_tests, 5))
View(out_sig)

# Larger population sizes were associated with more cases, earlier onset, longer duration, and more days with Reff>1
sort(out_sig[, "Pop2014"], decreasing = TRUE)

# More cases were associated with longer outbreaks, more cases of the other disease, more source-like behavior, earlier onset, more days with Reff>1
sort(out_sig[, "cholera"], decreasing = TRUE)
sort(out_sig[, "total_ebola"], decreasing = TRUE)

sort(out_sig[, "cholera"], decreasing = TRUE)

#### Save dataset ####
write.csv(df_cumulative_rcorr, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/df_cumulative_rcorr.csv")

#### CAR ####
# Cases
fcn_test(x = df_cumulative$cholera, y = df_cumulative$confirmed_ebola)

fcn_test(x = log10(1+df_cumulative$cholera), y = log10(1+df_cumulative$confirmed_ebola))

# Population
fcn_test(x = df_cumulative$Pop2012, y = df_cumulative$Pop2014)

# CAR
fcn_test(x = df_cumulative$cholera/df_cumulative$Pop2012, y = df_cumulative$confirmed_ebola/df_cumulative$Pop2014)

#### CAR, restricted to chiefdoms with both cholera and ebola ####
# Cases
fcn_test(x = df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "cholera"], y = df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "confirmed_ebola"])

fcn_test(x = log10(df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "cholera"]), y = log10(df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "confirmed_ebola"]))

# CAR
fcn_test(x = (df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "cholera"]/df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "Pop2012"]), y = (df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "confirmed_ebola"]/df_cumulative[df_cumulative$cholera > 0 & df_cumulative$confirmed_ebola > 0, "Pop2014"]))

#### Onset ####
# Onset
fcn_test(x = df_cumulative$cholera_onset, y = df_cumulative$confirmed_ebola_onset)

#### Duration ####
# Duration
ggplot(df_cumulative, aes(x = cholera_duration, y = confirmed_ebola_duration, size = Pop2014)) + 
  geom_point() +
  xlab("Duration of Cholera in Chiefdom (Days)") + 
  ylab("Duration of Ebola in Chiefdom (Days)") +
  scale_size_continuous(name = "Chiefdom\nPopulation") + 
  theme_bw()
fcn_test(x = df_cumulative$cholera_duration, y = df_cumulative$confirmed_ebola_duration)
fcn_test(x = df_cumulative[df_cumulative$cholera_duration > 0 & df_cumulative$confirmed_ebola_duration > 0, "cholera_duration"], y = df_cumulative[df_cumulative$cholera_duration > 0 & df_cumulative$confirmed_ebola_duration > 0, "confirmed_ebola_duration"])

#### Presence ####
# Chi-square test
lm(df_cumulative$confirmed_ebola > 0 ~ df_cumulative$cholera > 0)
chisq.test(df_cumulative$confirmed_ebola > 0, df_cumulative$cholera > 0)




