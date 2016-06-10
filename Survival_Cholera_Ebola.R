# Survival analysis for chiefdom cholera and ebola outbreaks

#### Load libraries and functions ####
library(survival)
source("~/Documents/2014 Cholera OCV/Data - Analysis/R Codes/ggsurv.R")
source("~/Documents/2014 Cholera OCV/Data - Analysis/R Codes/sle_cholera_2012/fcn_lookup.R")

#### Load data ####
cholera_cumulative <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/cholera_cumulative_GIS.csv")
cholera_cumulative$Cholera_Onset <- cholera_cumulative$Cholera_Onset - 5

ebola_cumulative <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_cumulative.csv")
ebola_cumulative$First_Ebola_Case <- ebola_cumulative$First_Ebola_Case + 1

#### Load Reff data ####
cholera_Reff <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/cholera_Reff_summary.csv")

ebola_Reff <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_Reff_summary.csv")

#### Create a master dataset ####
chiefs <- sort(unique(c(ebola_cumulative$CHCODE, cholera_cumulative$CHCODE)))
surv_data <- data.frame(CHCODE = chiefs)

surv_data$Cholera_Cases <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_cumulative[,"CHCODE"], value_column = cholera_cumulative[,"Cholera_Cases"]))
surv_data$t_dur_cholera <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_cumulative[,"CHCODE"], value_column = cholera_cumulative[,"t_dur_cholera"]))
surv_data$t_onset_cholera <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_cumulative[,"CHCODE"], value_column = cholera_cumulative[,"Cholera_Onset"]))
surv_data$days_Reff_over_1_Cholera <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_Reff[,"id"], value_column = cholera_Reff[,"days_Reff_over_1"]))

surv_data$Ebola_Cases <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_cumulative[,"CHCODE"], value_column = ebola_cumulative[,"Ebola_Cases"]))
surv_data$t_dur_ebola <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_cumulative[,"CHCODE"], value_column = ebola_cumulative[,"t_dur_ebola"]))
surv_data$t_onset_ebola <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_cumulative[,"CHCODE"], value_column = ebola_cumulative[,"First_Ebola_Case"]))
surv_data$days_Reff_over_1_Ebola <- apply(surv_data, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_Reff[,"id"], value_column = ebola_Reff[,"days_Reff_over_1"]))

#### Create a long data master ####
surv_data_long <- data.frame(CHCODE = c(ebola_cumulative$CHCODE, cholera_cumulative$CHCODE))
surv_data_long$Disease <- c(rep("Ebola", length(ebola_cumulative$CHCODE)), rep("Cholera", length(cholera_cumulative$CHCODE)))
surv_data_long$t_dur <- c(ebola_cumulative$t_dur_ebola, cholera_cumulative$t_dur_cholera)
surv_data_long$t_onset <- c(ebola_cumulative$First_Ebola_Case, cholera_cumulative$t_dur_cholera)

surv_data_long$Cholera_Cases <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_cumulative[,"CHCODE"], value_column = cholera_cumulative[,"Cholera_Cases"]))
surv_data_long$t_dur_cholera <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_cumulative[,"CHCODE"], value_column = cholera_cumulative[,"t_dur_cholera"]))
surv_data_long$t_onset_cholera <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_cumulative[,"CHCODE"], value_column = cholera_cumulative[,"Cholera_Onset"]))
surv_data_long$days_Reff_over_1_Cholera <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  cholera_Reff[,"id"], value_column = cholera_Reff[,"days_Reff_over_1"]))

surv_data_long$Ebola_Cases <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_cumulative[,"CHCODE"], value_column = ebola_cumulative[,"Ebola_Cases"]))
surv_data_long$t_dur_ebola <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_cumulative[,"CHCODE"], value_column = ebola_cumulative[,"t_dur_ebola"]))
surv_data_long$t_onset_ebola <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_cumulative[,"CHCODE"], value_column = ebola_cumulative[,"First_Ebola_Case"]))
surv_data_long$days_Reff_over_1_Ebola <- apply(surv_data_long, 1, function(x) fcn_lookup(query_1 = x["CHCODE"], query_2 = NA, reference =  ebola_Reff[,"id"], value_column = ebola_Reff[,"days_Reff_over_1"]))

surv_data_long$event <- 1

surv_data_long[surv_data_long$t_dur == 0,"t_dur"] <- NA
surv_data_long[surv_data_long$t_onset == 0,"t_onset"] <- NA

#### All chiefdoms ####
surv_duration_c <- Surv(time = surv_data$t_dur_cholera, event = rep(1, length(surv_data$t_dur_cholera)))

surv_duration_e <- Surv(time = surv_data$t_dur_ebola, event = rep(1, length(surv_data$t_dur_ebola)))

#### Chiefdoms with Reff>1 ####
surv_duration_c <- Surv(time = surv_data[surv_data$days_Reff_over_1_Cholera > 1,]$t_dur_cholera, event = rep(1, length(surv_data[surv_data$days_Reff_over_1_Cholera > 1,]$t_dur_cholera)))

surv_duration_e <- Surv(time = surv_data[surv_data$days_Reff_over_1_Ebola > 1,]$t_dur_ebola, event = rep(1, length(surv_data[surv_data$days_Reff_over_1_Ebola > 1,]$t_dur_ebola)))

#### Chiefdoms with >=25 cases ####
surv_duration_c <- Surv(time = surv_data[surv_data$Cholera_Cases >= 25,]$t_dur_cholera, event = rep(1, length(surv_data[surv_data$Cholera_Cases >= 25,]$t_dur_cholera)))

surv_duration_e <- Surv(time = surv_data[surv_data$Ebola_Cases >= 25,]$t_dur_ebola, event = rep(1, length(surv_data[surv_data$Ebola_Cases >= 25,]$t_dur_ebola)))

#### Chiefdoms with >=50 cases ####
surv_duration_c <- Surv(time = surv_data[surv_data$Cholera_Cases >= 50,]$t_dur_cholera, event = rep(1, length(surv_data[surv_data$Cholera_Cases >= 50,]$t_dur_cholera)))

surv_duration_e <- Surv(time = surv_data[surv_data$Ebola_Cases >= 50,]$t_dur_ebola, event = rep(1, length(surv_data[surv_data$Ebola_Cases >= 50,]$t_dur_ebola)))

#### Plot ####
# Cholera
fit1_c <- survfit(surv_duration_c ~ 1)
ggsurv(fit1_c) + ylab("Proportion of Cholera Outbreaks ongoing after x days") + xlab("days") +
  xlim(0, 425) + theme_bw()

# Ebola
fit1_e <- survfit(surv_duration_e ~ 1)
ggsurv(fit1_e) + ylab("Proportion of Ebola Outbreaks ongoing after x days") + xlab("days") +
  xlim(0, 425) + theme_bw()

#### Compare survival curves for t_dur for Ebola and Cholera (all cases) ####
surv_data_long_use = surv_data_long

survdiff(Surv(time = t_dur, event = event) ~ Disease, data = surv_data_long_use)  # log rank test
p = 1-pchisq(survdiff(Surv(time = t_dur, event = event) ~ Disease, data = surv_data_long_use)$chisq, 1)

ggsurv(survfit(Surv(time = t_dur, event = event) ~ Disease, data = surv_data_long_use)) + ylab("Proportion of Outbreaks Ongoing") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since Chiefdom Onset")  + ggtitle(paste("Epidemic Duration"))  + annotate("text", x = 10*30, y = 0.8, label = paste("Log-Rank p = ", as.character(round(p, 4))), size = 2) + theme(text = element_text(size=6))

ggsave(filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/Survival_EpiDur_all.pdf", width = 3, height = 2 )

ggsurv(survfit(Surv(time = t_dur, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Cholera",])) + ylab("Proportion of Outbreaks Ongoing") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since Chiefdom Onset") + ggtitle(paste("Cholera Duration"))  + theme(text = element_text(size=6))

ggsurv(survfit(Surv(time = t_dur, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Ebola",])) + ylab("Proportion of Outbreaks Ongoing") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since Chiefdom Onset") + ggtitle(paste("Ebola Duration")) + theme(text = element_text(size=6))

#### Compare survival curves for t_dur for Ebola and Cholera (min Reff) ####
x = 0
surv_data_long_use = surv_data_long[(surv_data_long$Disease == "Ebola" & surv_data_long$Ebola_Cases > x) | (surv_data_long$Disease == "Cholera" & surv_data_long$Cholera_Cases > x),]

min_days_Reff_over_1 = 2
surv_data_long_use = surv_data_long[(surv_data_long$Disease == "Ebola" & surv_data_long$days_Reff_over_1_Ebola >= min_days_Reff_over_1) | (surv_data_long$Disease == "Cholera" & surv_data_long$days_Reff_over_1_Cholera >= min_days_Reff_over_1),]

survdiff(Surv(time = t_dur, event = event) ~ Disease, data = surv_data_long_use)  # log rank test
p = 1-pchisq(survdiff(Surv(time = t_dur, event = event) ~ Disease, data = surv_data_long_use)$chisq, 1)

ggsurv(survfit(Surv(time = t_dur, event = event) ~ Disease, data = surv_data_long_use)) + ylab("Proportion of Outbreaks Ongoing") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since Chiefdom Onset")  + ggtitle(paste("Epidemic Duration, at least", min_days_Reff_over_1, "Day(s) of Reff > 1"))  + annotate("text", x = 10*30, y = 0.8, label = paste("Log-Rank p = ", as.character(round(p, 4))), size = 2) + theme(text = element_text(size=6))

ggsave(filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/Survival_EpiDur_ReffMin.pdf", width = 3, height = 2 )

ggsurv(survfit(Surv(time = t_dur, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Cholera",])) + ylab("Proportion of Outbreaks Ongoing") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since Chiefdom Onset") + ggtitle(paste("Cholera Onset, at least", min_days_Reff_over_1, "Day(s) of Reff > 1")) + theme(text = element_text(size=6))

ggsurv(survfit(Surv(time = t_dur, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Ebola",])) + ylab("Proportion of Outbreaks Ongoing") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since Chiefdom Onset") + ggtitle(paste("Ebola Onset, at least", min_days_Reff_over_1, "Day(s) of Reff > 1")) + theme(text = element_text(size=6))

#### Compare survival curves for t_onset for Ebola and Cholera (all cases) ####
summary(ebola_cumulative$First_Ebola_Case)
summary(cholera_cumulative$Cholera_Onset)
hist(cholera_cumulative$Cholera_Onset, breaks = 40)

surv_data_long_use = surv_data_long

t_onset = 0 # ~160 If you want to restrict the cholera outbreak to only the main epidemic
surv_data_long_use = surv_data_long_use[surv_data_long_use$Disease == "Ebola" | (surv_data_long_use$Disease == "Cholera" & surv_data_long_use$t_onset > t_onset),]
surv_data_long_use[surv_data_long_use$Disease == "Cholera" & is.na(surv_data_long_use$t_onset) == 0,"t_onset"] = surv_data_long_use[surv_data_long_use$Disease == "Cholera" & is.na(surv_data_long_use$t_onset) == 0,"t_onset"] - t_onset

survdiff(Surv(time = t_onset, event = event) ~ Disease, data = surv_data_long_use)  # log rank test
p = 1-pchisq(survdiff(Surv(time = t_onset, event = event) ~ Disease, data = surv_data_long_use)$chisq, 1)

ggsurv(survfit(Surv(time = t_onset, event = event) ~ Disease, data = surv_data_long_use)) + ylab("Proportion of Outbreaks Yet to Commence") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since first national case")  + ggtitle(paste("Epidemic Onset")) + annotate("text", x = 10*30, y = 0.8, label = paste("Log-Rank p = ", as.character(round(p, 4))), size = 2) + theme(text = element_text(size=6))

ggsave(filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/Survival_EpiOnset_all.pdf", width = 3, height = 2 )

ggsurv(survfit(Surv(time = t_onset, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Cholera",])) + ylab("Proportion of Outbreaks Yet to Commence") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months") + ggtitle(paste("Cholera Onset")) + theme(text = element_text(size=6))

ggsurv(survfit(Surv(time = t_onset, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Ebola",])) + ylab("Proportion of Outbreaks Yet to Commence") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months since first national case") + ggtitle(paste("Ebola Onset")) + theme(text = element_text(size=6))

#### Compare survival curves for t_onset for Ebola and Cholera (cholera outbreaks in Fall) ####
summary(ebola_cumulative$First_Ebola_Case)
summary(cholera_cumulative$Cholera_Onset)
hist(cholera_cumulative$Cholera_Onset, breaks = 40)

x = 0
surv_data_long_use = surv_data_long[(surv_data_long$Disease == "Ebola" & surv_data_long$Ebola_Cases > x) | (surv_data_long$Disease == "Cholera" & surv_data_long$Cholera_Cases > x),]

t_onset = 152 # ~152 If you want to restrict the cholera outbreak to only the main epidemicbeginning in June
surv_data_long_use = surv_data_long_use[surv_data_long_use$Disease == "Ebola" | (surv_data_long_use$Disease == "Cholera" & surv_data_long_use$t_onset > t_onset),]
surv_data_long_use[surv_data_long_use$Disease == "Cholera" & is.na(surv_data_long_use$t_onset) == 0,"t_onset"] = surv_data_long_use[surv_data_long_use$Disease == "Cholera" & is.na(surv_data_long_use$t_onset) == 0,"t_onset"] - t_onset

survdiff(Surv(time = t_onset, event = event) ~ Disease, data = surv_data_long_use)  # log rank test
p = 1-pchisq(survdiff(Surv(time = t_onset, event = event) ~ Disease, data = surv_data_long_use)$chisq, 1)

ggsurv(survfit(Surv(time = t_onset, event = event) ~ Disease, data = surv_data_long_use)) + ylab("Proportion of Outbreaks Yet to Commence") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months") + ggtitle(paste("Epidemic Onset, cholera cases beginning in June")) + annotate("text", x = 10*30, y = 0.8, label = paste("Log-Rank p = ", as.character(round(p, 4))), size = 2) + theme(text = element_text(size=6))

ggsave(filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/Survival_EpiOnset_FallCholera.pdf", width = 3, height = 2 )

ggsurv(survfit(Surv(time = t_onset, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Cholera",])) + ylab("Proportion of Outbreaks Yet to Commence") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months") + ggtitle(paste("Cholera Onset, cholera cases beginning in June")) + theme(text = element_text(size=6))

ggsurv(survfit(Surv(time = t_onset, event = event) ~ 1, data = surv_data_long_use[surv_data_long_use$Disease == "Ebola",])) + ylab("Proportion of Outbreaks Yet to Commence") +  theme_bw() + scale_x_continuous(breaks = seq(0, 420, 30), limits = c(0, 425), labels = seq(0, 14, 1), name = "Months") + ggtitle(paste("Cholera Onset, cholera cases beginning in June")) + theme(text = element_text(size=6))


