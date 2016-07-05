#### Header ####
# Use signed-rank test to compare epidemic attributes

#### Load libraries ####
library(Hmisc)

#### Notes ####
# The cholera and ebola epidemics in chiefdoms differed with respect to cases, CAR, order of onset, presence of cases.
# Among chiefdoms that had BOTH cholera and ebola, there was no correlation between the duration of the outbreaks
# If we say a chiefdom with 0 cases had an outbreak of 0 duration, then the p-value inflates to 0.06376, so we conclude that the durations are similar, barely.

#### Load data ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")

#### Function ####
fcn_test <- function(x, y, plot = TRUE){
  scatterplot = plot(x, y)
  p = round(rcorr(x,y, type="pearson")$P[1,2], 4)

  if (plot){scatterplot}
  
  return(c(p = p))
}

#### Test function ####
fcn_test(x = seq(1,20), y = seq(1,20))
fcn_test(x = seq(1,20), y = seq(1,20) + rnorm(20, mean = 10))
fcn_test(x = seq(1,20), y = rnorm(20, mean = 5))

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




