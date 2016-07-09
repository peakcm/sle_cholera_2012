#### Header ####
# Use signed-rank test to compare epidemic attributes

#### Load libraries ####
library(Hmisc)
library(dplyr)
library(tidyr)
library(corrplot)

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

#### Add cumulative attack rates ####
df_cumulative$total_ebola_CAR <- df_cumulative$total_ebola / df_cumulative$Pop2014
df_cumulative$cholera_CAR <- df_cumulative$cholera / df_cumulative$Pop2012

#### Add Coastal indicator ####
coastal_chiefs = c("2206", "2406", "2404", "2405", "4199", "4299", "3312", "3302", "3305", "3313", "3209", "3203", "3208", "3207", "3407", "3403", "3410")
df_cumulative <- df_cumulative %>% mutate(Coastal = ifelse(CHCODE %in% coastal_chiefs, yes = 1, no = 0))

#### Select fields ####
demographic_fields <- c( "Pop2012","Prop_Urban", "Education", "Prop_Improved_Water_Source", "Prop_Improved_Toilet", "SES_Score", "Household_Size", "Shared_Toilet", "Wealth_Index", "weighted_density_per_km", "Pop_per_sq_km", "Coastal")
disease_fields <- c("total_ebola", "cholera", "total_ebola_onset", "cholera_onset", "total_ebola_duration", "cholera_duration", "total_ebola_CAR", "cholera_CAR")

df_cumulative_rcorr <- df_cumulative %>%
  mutate(total_ebola_CAR = total_ebola / Pop2014,
         cholera_CAR = cholera / Pop2012) %>% apply(2, as.numeric)

#### Look at all pairs ####
df_cumulative_rcorr <- df_cumulative_rcorr[,c("CHCODE", demographic_fields, disease_fields)]

df_cumulative_rcorr <- data.frame(df_cumulative_rcorr) %>% 
  left_join(ebola_Reff_summary[,c("CHCODE", "days_Reff_over_1_ebola", "export_import_ratio_ebola")], by = "CHCODE") %>%
  left_join(cholera_Reff_summary[,c("CHCODE", "days_Reff_over_1_cholera", "export_import_ratio_cholera")], by = "CHCODE")

names(df_cumulative_rcorr) <- c("CHCODE", "Population", "Urban", "Education", "Improved Water", "Improved Sanitation", "SES", "Household Size","Shared Toilet", "Wealth Index", "Density (Weighted)", "Density","Coastal", "Ebola Cases", "Cholera Cases", "Ebola Onset", "Cholera Onset", "Ebola Duration", "Cholera Duration", "Ebola Attack Rate", "Cholera Attack Rate", "Ebola Days R>1", "Ebola Export/Import Ratio", "Cholera Days R>1", "Cholera Export/Import Ratio")

out <- rcorr(as.matrix(df_cumulative_rcorr), type = "spearman")
View(out$r)
plot(out$r[,"Density"], out$r[,"Density (Weighted)"]) # Density metrics behave similarly

#### Function to prepare for plotting ####
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
FlattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

out_flat <- FlattenCorrMatrix(out$r, out$P)

#### corrplot of disease-disease associations ####
selection <- c("Cholera Cases", "Cholera Attack Rate", "Cholera Duration", "Cholera Onset", "Ebola Cases", "Ebola Attack Rate", "Ebola Duration", "Ebola Onset")

keepers <- as.numeric(sapply(selection, function(x) which(row.names(out$r) %in% x)))

# For correlations
subset <- data.frame(out$r, row.names = NULL)
subset <- subset[keepers, keepers]
row.names(subset) <- selection
names(subset) <- selection
subset <- as.matrix(subset)

# For p-values
subset.p <- data.frame(out$P, row.names = NULL)
subset.p <- subset.p[keepers, keepers]
subset.p <- as.matrix(subset.p)

n_tests_disease_disease =  length(selection)^2/2
# n_tests_disease_demographic = 0
sig.level = 0.05 / sum(n_tests_disease_disease + n_tests_disease_demographic)

setwd("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/")
pdf("20160709_Disease_Disease_Spearman.pdf", 5,5)
corrplot(subset, type="lower",  tl.col="black", tl.srt=45, p.mat = subset.p, sig.level = sig.level, insig = "blank",
         tl.cex = 0.5, cl.cex = 0.4, mar = c(1,1,5,1))
dev.off()

#### corrplot of disease-demographic associations ####
selection_disease <- c("Cholera Cases", "Cholera Attack Rate", "Cholera Duration", "Cholera Onset", "Ebola Cases", "Ebola Attack Rate", "Ebola Duration", "Ebola Onset")
keepers_disease <- as.numeric(sapply(selection_disease, function(x) which(row.names(out$r) %in% x)))

selection_demographic <- c("Population","Density", "Education", "Urban", "Improved Water", "Improved Sanitation", "Household Size", "Coastal")
keepers_demographic <- as.numeric(sapply(selection_demographic, function(x) which(row.names(out$r) %in% x)))

selection <- c(selection_disease, selection_demographic)
keepers <- as.numeric(sapply(selection, function(x) which(row.names(out$r) %in% x)))

# For correlations
subset <- data.frame(out$r, row.names = NULL)
subset <- subset[keepers_disease, keepers_demographic]
row.names(subset) <- selection_disease
names(subset) <- selection_demographic
subset <- as.matrix(subset)

# For p-values
subset.p <- data.frame(out$P, row.names = NULL)
subset.p <- subset.p[keepers_disease, keepers_demographic]
subset.p <- as.matrix(subset.p)

# n_tests_disease_disease =  0
n_tests_disease_demographic = (length(selection_disease)*length(selection_demographic))
sig.level = 0.05 / sum(n_tests_disease_disease + n_tests_disease_demographic)

setwd("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/")
pdf("20160709_Disease_Demographic_Spearman.pdf", 5,5)
corrplot(subset, tl.col="black", tl.srt=45, p.mat = subset.p, sig.level = sig.level, insig = "blank",cl.pos = "n",
         tl.cex = 0.5, cl.cex = 0.4, diag = FALSE, mar = c(5,1,3,1))
dev.off()

### Convert to long format ####
# Adjust for multiple testing
num_disease_metrics <- length(disease_fields)
num_DHS_metrics <- length(demographic_fields) - 1
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

out_sig <- data.frame(out_sig)
out_sig$row.name <- row.names(out_sig)
out_sig_gather <- gather(out_sig, key = row.name)
names(out_sig_gather) <- c("Var1", "Var2", "corr")

out_sig_gather <- out_sig_gather %>% filter(is.na(corr) == 0, 
                                            Var1 != Var2,
                                            Var1 != "CHCODE",
                                            Var2 != "CHCODE")
out_sig_gather <- out_sig_gather[0 == (out_sig_gather$Var1 %in% demographic_fields & out_sig_gather$Var2 %in% demographic_fields),] # Not interested in correlations between demographic fields

# Add indicators for cross-disease and demographic-disesease combos
out_sig_gather <- out_sig_gather %>%
  mutate(cross_disease = ifelse((grepl("cholera", Var1) & grepl("ebola", Var2)) | (grepl("ebola", Var1) & grepl("cholera", Var2)), yes = 1, no = 0)) %>%
  mutate(demo_disease = ifelse(xor(xor(grepl("cholera", Var1), grepl("ebola", Var1)), xor(grepl("cholera", Var2), grepl("ebola", Var2))), yes = 1, no = 0))

# View results
View(out_sig_gather)

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




