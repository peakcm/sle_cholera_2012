# White Adaptation to Wallinga Teunis method

#### Load Libraries ####
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggvis)
library(forecast)
library(RColorBrewer)
library(simcf)     # For panel functions and simulators and lag functions

#### Load Functions ####
source("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial.R")

#### Helper functions ####
fcn_lookup <- function(query_1, query_2 = NA, reference, value_column, transformation = "none"){
  out <- 0
  if (is.na(query_2)==0){
    row_1 <- which(reference[,1] %in% query_1) 
    row_2 <- which(reference[,2] %in% query_2)
    if (length(row_1) > 0 & length(row_2) > 0){
      out <- value_column[intersect(row_1, row_2)]
    }
  } else {
    out <- value_column[which(reference == query_1)]
  }
  if (length(out)==1){
    if (transformation == "as.character"){
      return(as.character(out))
    } else {return(out)}
  } else {
    return(0)
  }
}

#### Load workspace ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial_Cholera.RData")

#### Load cholera data ####
setwd('/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/Admin 3/Admin 3 Daily')

cumulative <- read.csv('Cumulative_Data_Chiefdom.csv')
cumulative <- cumulative[order(-cumulative$Cumulative_Cases),]

cholera_df <- read.csv('Daily_Cases.csv')
cholera_df$Day <- as.Date(as.character(as.POSIXct(cholera_df$Day, format = "%m/%d/%y")))
cholera_df <- cholera_df[cholera_df$Day < "2013-02-21",]

cholera_df$AR10000 <- NA
for (chief in unique(cholera_df$CHCODE)){
  cholera_df[cholera_df$CHCODE == chief, "AR10000"] <- cholera_df[cholera_df$CHCODE == chief, "Case"] / cumulative[cumulative$CHCODE == chief, "Population"] * 10000
}

#### Load df_daily ####
df_daily <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_cholera_daily.csv")

df_daily$date <- as.Date(df_daily$date, format = "%Y-%m-%d")

df_daily <- df_daily[df_daily$date < as.Date("2013-03-01"),]

#### Spatial networks ####

# Option 1: Gravity model from Wesolowski PLOS Currents Outbreaks
gravity_model <- function(pop_i, pop_j, dist_ij, k = 3.93, alpha=0.47, beta=0.46, gamma=-1.78){
  k*(pop_i^alpha + pop_j^beta) * dist_ij^gamma
}

setwd('/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Spatial Analysis/Neighborhood Matrices/Admin 3/')
neighborhood_ <- read.csv("Distance_km_nothresh_Network.csv")
names(neighborhood_) <- sapply(names(neighborhood_), function(x) substring(x, 2))
neighborhood_$CHCODE <- as.numeric(names(neighborhood_))
chiefs <- neighborhood_$CHCODE

weights <- melt(neighborhood_, id.vars = "CHCODE")
names(weights) <- c("CHCODE_source", "CHCODE_destination", "Distance_km")
weights$CHCODE_destination <- as.numeric(as.character(weights$CHCODE_destination))
weights <- weights[order(weights$CHCODE_source, weights$CHCODE_destination),]

weights$Pop_source <- sapply(weights$CHCODE_source, function(x) cumulative[cumulative$CHCODE == x, "Population"])
weights$Pop_destination <- sapply(weights$CHCODE_destination, function(x) cumulative[cumulative$CHCODE == x, "Population"])

weights$gravity_Wesolowski <- apply(weights, 1, function(x) gravity_model(x["CHCODE_source"], x["CHCODE_destination"], x["Distance_km"]))

row.names(weights) <- 1:nrow(weights)

# Option 2: Inverse Distance
weights$inverse_distance <- sapply(weights$Distance_km, function(x) min(1,max(0, 1-x/max(weights$Distance_km))))
weights$inverse_distance_squared <- weights$inverse_distance^2

# Option 3: CDR
setwd('/Users/peakcm/Documents/SLE_Mobility/sle_ericson/')
neighborhood_CDR <- read.csv("data_3day_melt_chiefdom_na.rm.csv")
neighborhood_CDR <- melt(neighborhood_CDR, id.vars = c("admin3_current", "admin3_previous", "date"))
neighborhood_CDR <- neighborhood_CDR[neighborhood_CDR$variable %in% c("count"),]
neighborhood_CDR$value <- as.numeric(neighborhood_CDR$value)
neighborhood_CDR <- dcast(neighborhood_CDR, admin3_current + admin3_previous  ~ variable, fun.aggregate = sum)

CDR_chiefs <- unique(weights[is.na(weights$CDR)==0,"CHCODE_source"])

    # Consider just excluding chiefdoms without CDR data. Only 1693 cases there. 7.5% of total cases
    summary(cumulative[cumulative$CHCODE %in% CDR_chiefs,c("Population","Cumulative_Cases")])
    sum(cumulative[cumulative$CHCODE %in% CDR_chiefs,c("Cumulative_Cases")])
    
    summary(cumulative[(cumulative$CHCODE %in% CDR_chiefs)==0,c("Population","Cumulative_Cases")])
    sum(cumulative[(cumulative$CHCODE %in% CDR_chiefs)==0,c("Cumulative_Cases")])
    
    # Consider making a regression model to predict travel related to cheifdoms with missing CDR data

weights$CDR <- NA
for (i in 1:nrow(weights)){
  if (weights[i, c("CHCODE_source")] %in% neighborhood_CDR[,c("admin3_previous")]){
    if (weights[i, c("CHCODE_destination")] %in% neighborhood_CDR[neighborhood_CDR$admin3_previous == weights[i, c("CHCODE_source")],c("admin3_current")]){
      weights[i,"CDR"] <- neighborhood_CDR[neighborhood_CDR$admin3_previous == weights[i, c("CHCODE_source")] & neighborhood_CDR$admin3_current == weights[i, c("CHCODE_destination")],"count"]
    }
  }
}

weights$CDR_source_pop_weighted <- weights$CDR / weights$Pop_source
    
# Option 3: No connections between chiefdoms
weights$identity <- 0
weights[weights$CHCODE_source == weights$CHCODE_destination, "identity"] <- 1

#### Generation Interval ####
# Azman generation time from EID 2016
summary(rgamma(100000, rate=0.1, shape = 0.5))
tau <- (hist(rgamma(100000, rate=0.1, shape = 0.5), breaks = 0:1000)$counts)/100000
plot(tau)
plot(cumsum(tau), xlim = c(0,30), type = "l")

#### Play dataset ####
days <- 1:6
id.var_play <- rep(c("A", "B"), each = length(days))
date.var_play <- rep(days, length(unique(id.var)))
case.var_play <- c(0,20,40,80,20,0,
                     0,0,0,1,0,0)
tau_play <- c(0.1, 0.3, 0.2, 0.2, 0.1, 0.1)

weight_play <- data.frame(id_source = c("A","A","B","B"), id_destination = c("A", "B", "A", "B"), weight = c(1,0.5,0.5,1))

weight_play <- data.frame(id_source = c("A","A","B","B"), id_destination = c("A", "B", "A", "B"), weight = c(1,1,1,1))

play_output <- WT_Spatial(id.var = id.var_play, date.var = date.var_play, case.var = case.var_play, weights = weight_play, generation_interval = tau_play)

#### Select Data to use ####
# Chiefdoms with any cases
chiefs <- cumulative[order(-cumulative$Cumulative_Cases),"CHCODE"]

# Chiefdoms represented in the raw CDR dat
chiefs <- unique(weights[is.na(weights$CDR)==0,"CHCODE_source"])

# Chiefdoms with high case loads
chiefs <- cumulative[order(-cumulative$Cumulative_Cases),"CHCODE"][1:10]

#### Use WT_Spatial to calculate Reff Function ####
cholera_df_use <- WT_Spatial(id.var = cholera_df$CHCODE, date.var = cholera_df$Day, case.var = cholera_df$Case, weights = weights[,c("CHCODE_source", "CHCODE_destination", "inverse_distance_squared")], generation_interval = tau)

#### Summarize Reff calculations ####
cholera_df_use_summary <- WT_Spatial_Summary(cholera_df_use)

View(cholera_df_use_summary)
hist(cholera_df_use_summary$days_Reff_over_1, breaks = 20)
hist(cholera_df_use_summary$export_import_ratio, breaks = 20)

names(cholera_df_use)[1:3] <- c("CHCODE", "Day", "Case")

#### Moving Averages ####
# R-effective 7-day moving average
cholera_df_use = cholera_df_use %>% group_by(CHCODE) %>%
  mutate(Reff_7dayMA = ma(Reff, order = 7))

#### Add rainfall data ####
cholera_df_use$rain_avg <- NA
cholera_df_use$rain_avg_lag7 <- NA
cholera_df_use$area <- NA

for (chief in unique(cholera_df_use$CHCODE)){
  days <- data.frame(cholera_df_use)[cholera_df_use$CHCODE == chief,"Day"]
  rain <- df_daily[df_daily$CHCODE == chief & df_daily$date %in% days,"rain_avg"]
  area <- df_daily[df_daily$CHCODE == chief & df_daily$date %in% days,"area"]
  if (length(rain)==length(days)){
    cholera_df_use[cholera_df_use$CHCODE == chief, "rain_avg"] <- rain
    cholera_df_use[cholera_df_use$CHCODE == chief, "rain_avg_lag7"] <- c(rain[8:length(days)], rep(NA,7))
    cholera_df_use[cholera_df_use$CHCODE == chief, "area"] <- area
    cat(".")
  } else {cat("Error\n")}
}

summary(cholera_df_use$rain_avg)
summary(cholera_df_use$rain_avg_lag7)

cholera_df_use = cholera_df_use %>% group_by(CHCODE) %>%
  mutate(rain_7dayMA = ma(rain_avg, order = 7)) %>%
  mutate(rain_7dayMA_lag7 = ma(rain_7dayMA, order = 7))

ggplot(cholera_df_use, aes(x = Day, group = CHCODE, color = CHCODE)) +
  geom_point(aes(y = rain_avg)) +
  geom_line(aes(y = Reff_7dayMA))

ggplot(cholera_df_use, aes(x = rain_7dayMA_lag7, y = Reff_7dayMA)) +
  geom_point() 

#### Aggregate by Region or Nation ####
# Add region data
cholera_df_use$region <- NA
cholera_df_use[cholera_df_use$CHCODE < 2000,"region"] <- "East"
cholera_df_use[cholera_df_use$CHCODE > 2000 & cholera_df_use$CHCODE < 3000,"region"] <- "North"
cholera_df_use[cholera_df_use$CHCODE > 3000 & cholera_df_use$CHCODE < 4000,"region"] <- "South"
cholera_df_use[cholera_df_use$CHCODE > 4000,"region"] <- "West"

cholera_df_use$region <- factor(cholera_df_use$region, levels = c("West", "South", "North", "East"))

# Region average
cholera_df_use_region <- data.frame(cbind(Day = rep(as.Date(as.character(unique(cholera_df_use$Day))),4)))
cholera_df_use_region$region <- rep(c("East", "North", "South", "West"), each = length(unique(cholera_df_use$Day)))
cholera_df_use_region$Case <- NA
cholera_df_use_region$Reff_weighted <- NA
cholera_df_use_region$rain_weighted <- NA

for (i in 1:nrow(cholera_df_use_region)){
  day <- as.Date(cholera_df_use_region[i,"Day"])
  if (day == min(cholera_df_use_region$Day)){
    region <- cholera_df_use_region[i,"region"][1]
  }
  Reff <- unlist(cholera_df_use$Reff[cholera_df_use$Day %in% day & cholera_df_use$region %in% region])
  Case <- unlist(cholera_df_use$Case[cholera_df_use$Day %in% day & cholera_df_use$region %in% region])
  rain_avg <-unlist(cholera_df_use$rain_avg[cholera_df_use$Day %in% day & cholera_df_use$region %in% region]) 
  area <- unlist(cholera_df_use$area[cholera_df_use$Day %in% day & cholera_df_use$region %in% region]) 
  
  cholera_df_use_region[i, "Case"] <- sum(Case)
  if (sum(Case) > 0){
    cholera_df_use_region[i, "Reff_weighted"] <- weighted.mean(x = Reff , w = Case)
  } else {cholera_df_use_region[i, "Reff_weighted"] <- 0}
  cholera_df_use_region[i,"rain_weighted"] <- weighted.mean(x = rain_avg , w = area)
}

cholera_df_use_region$Reff_weighted_7dayMA <- ma(cholera_df_use_region$Reff_weighted, order = 7)
cholera_df_use_region$rain_weighted_7dayMA <- ma(cholera_df_use_region$rain_weighted, order = 7)

cholera_df_use_region$Day <- as.Date(cholera_df_use_region$Day)

cholera_df_use_region$region <- factor(cholera_df_use_region$region, levels = c("West", "South", "North", "East"))

# Country average
cholera_df_use_country <- data.frame(cbind(Day = as.Date(as.character(unique(cholera_df_use$Day)))))
cholera_df_use_country$Case <- NA
cholera_df_use_country$Reff_weighted <- NA
cholera_df_use_country$rain_weighted <- NA

for (i in 1:nrow(cholera_df_use_country)){
  
  day <- as.Date(cholera_df_use_country[i,"Day"])
  
  Reff <- unlist(cholera_df_use$Reff[cholera_df_use$Day %in% day])
  Case <- unlist(cholera_df_use$Case[cholera_df_use$Day %in% day])
  rain <- unlist(cholera_df_use$rain_avg[cholera_df_use$Day %in% day])
  area <- unlist(cholera_df_use$area[cholera_df_use$Day %in% day])
  
  cholera_df_use_country[i, "Case"] <- sum(Case)
  if (sum(Case) > 0){
    cholera_df_use_country[i, "Reff_weighted"] <- weighted.mean(x = Reff , w = Case)
    cholera_df_use_country[i, "rain_weighted"] <- weighted.mean(x = rain , w = area)
  } else {cholera_df_use_country[i, "Reff_weighted"] <- 0}
}

cholera_df_use_country$Reff_weighted_7dayMA <- ma(cholera_df_use_country$Reff_weighted, order = 7)
cholera_df_use_country$rain_weighted_7dayMA <- ma(cholera_df_use_country$rain_weighted, order = 7)

cholera_df_use_country$Day <- as.Date(cholera_df_use_country$Day)

#### Days with Rt above 1 ####
View(cholera_df_use[cholera_df_use$Reff > 1,])
length(unlist(unique(cholera_df_use[cholera_df_use$Reff > 1,"CHCODE"])))
# 106 chiefdoms had Reff > 1 (with identity matrix)
# 37 chiefdoms had Reff > 1 (with inverse_distance_squared matrix)

ggplot(cholera_df_use[cholera_df_use$Reff > 1,], aes(x=Day)) + geom_bar() + ggtitle("Number of Chiefdoms with Reff > 1") 
# At the peak, 16 chiefdoms had Reff>1 (with identity matrix)
# At the peak, 15 chiefdoms had Reff>1 (with inverse_distance_squared matrix), At the left tail <5, right tail, <3 chiefdoms had Reff>1

ggplot(cholera_df_use[cholera_df_use$Case > 0,], aes(x=Day)) + geom_bar() + ggtitle("Number of Chiefdoms with Cases") 
# At the peak, between 30 and 44 chiefdoms had cases simultanously

table <- data.frame(table(cholera_df_use[cholera_df_use$Reff > 1,]$CHCODE))
ggplot(table, aes(x = factor(Var1), y = Freq)) + geom_bar(stat = "identity") + xlab("CHCODE") + ggtitle("Number of days with Reff > 1 by Chiefdom")

#### Store Data ####
# cholera_df_use_inverse_distance_squared <- cholera_df_use
# cholera_df_use_identity <- cholera_df_use

#### Select Data ####
# cholera_df_use <- cholera_df_use_inverse_distance_squared
# cholera_df_use <- cholera_df_use_identity

#### Chiefdom Specific Plot ####
ggplot() +
  theme_bw() +
  geom_bar(data = cholera_df_use, aes(x=Day, y = AR10000, group = factor(CHCODE), fill = factor(CHCODE)), stat = "identity", position = "dodge", alpha = 0.5, color = 0) +
  geom_line(data = cholera_df_use[cholera_df_use$Reff_7dayMA > 0 ,], aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), color = factor(CHCODE))) +
  geom_hline(yintercept = 1, col = "black", lty = "dashed") +
  xlab("Day") + 
  ylim(c(0, 3)) +
  # ylab("Reproductive Number R") + 
  ggtitle("Daily Reproductive Number") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

select_chiefs <- c("4299","4199","3410", "2406", "2511", "1212")
ggplot() +
#   geom_bar(data = cholera_df_use[cholera_df_use$CHCODE %in% select_chiefs & cholera_df_use$Reff_7dayMA > 0 ,], aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), fill = factor(region)), stat = "identity", position = "dodge", width = 1, alpha = 1) +
  geom_line(data = cholera_df_use[cholera_df_use$Reff_7dayMA > 0.8 ,], aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), color = factor(region))) +
  # facet_grid(region~.)+
  geom_hline(yintercept = 1, col = "black", lty = "dashed") +
  xlab("Day") + 
  ylim(c(0.8, 2)) +
  # ylab("Reproductive Number R") + 
  ggtitle("Daily Reproductive Number") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#### Regional Rt curves ####
regional_Rt <- ggplot() +
  theme_bw() +

  geom_area(data = cholera_df_use_region, aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = .8, color = "black", size = 0.2) +
  geom_hline(yintercept = 1, col = "black", lty = "dashed", size = 0.3) +
  geom_line(data = cholera_df_use_region, aes(x = Day, y=2*rain_weighted_7dayMA/max(cholera_df_use_region$rain_weighted_7dayMA, na.rm = TRUE), group = factor(region)), size = .2, color = "blue", alpha = 0.3) +
  
  xlab("Day") + 
  facet_grid(region~.) +
  scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b", name = "Month") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(R[t])) +
  scale_fill_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  theme(text = element_text(size=8)) +
  guides(fill = "none")
regional_Rt

ggsave(plot = regional_Rt, filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/Regional_Rt.pdf", width = 3, height = 1.5 )

ggplot() +
  theme_bw() +
  geom_area(data = cholera_df_use_region[cholera_df_use_region$region == "North",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  geom_area(data = cholera_df_use_region[cholera_df_use_region$region == "East",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  geom_area(data = cholera_df_use_region[cholera_df_use_region$region == "West",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  geom_area(data = cholera_df_use_region[cholera_df_use_region$region == "South",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  
  geom_line(data = cholera_df_use_region[cholera_df_use_region$region == "North",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)), alpha = 0.8, size = 1.25) +
  geom_line(data = cholera_df_use_region[cholera_df_use_region$region == "East",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)), alpha = 0.8, size = 1.25) +
  geom_line(data = cholera_df_use_region[cholera_df_use_region$region == "West",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)),  alpha = 0.8, size = 1.25) +
  geom_line(data = cholera_df_use_region[cholera_df_use_region$region == "South",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)), alpha = 0.8, size = 1.25) +
  
  # geom_line(data = cholera_df_use, aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), color = factor(region)), size = .3, color = "white", alpha = 0.4) +
  geom_hline(yintercept = 1, col = "black", lty = "dashed") +
  xlab("Day") + 
  ylim(c(0.0, 2)) +
  # ylab("Reproductive Number R") + 
  # ggtitle("Daily Reproductive Number") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(R[t], "7 Day Moving Average")) +
  scale_fill_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  scale_color_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  guides(fill = "none")

#### National Rt curve ####
ggplot() +
  theme_bw() +
  geom_bar(data = cholera_df_use_country, aes(x = Day, y = Case/sum(cumulative$Population)*100000), stat = "identity", alpha = 0.5) +
  geom_line(data = cholera_df_use_country[cholera_df_use_country$Reff_weighted_7dayMA > 0,], aes(x = Day, y = Reff_weighted_7dayMA)) +
  ylim(0, 7)+
  geom_hline(yintercept = 1, col = "black", lty = "dashed") +
  ylab("Effective Reproductive Numer (7-day Moving Average)\nCases per 100,000 Population")

cholera_df_use_region_weekly <- cholera_df_use_region %>% 
  group_by(region, Week = cut(Day, "week")) %>% mutate(Case_weekly = sum(Case))

national_Rt_epicurve <- ggplot() +
  theme_bw() +
  geom_area(data = cholera_df_use_region_weekly, aes(x = Day, y = Case_weekly/1000, fill = region), stat = "identity", alpha = 0.8) +
  geom_line(data = cholera_df_use_country[cholera_df_use_country$Reff_weighted_7dayMA > 0,], aes(x = Day, y = Reff_weighted_7dayMA), size = 0.3) +
  # geom_line(data = cholera_df_use_country, aes(x = Day, y = 2*rain_weighted_7dayMA/max(cholera_df_use_country$rain_weighted_7dayMA, na.rm=TRUE)), size = 0.3, color = "blue") +
  scale_x_date(date_breaks = "2 month", date_labels = "%b", name = "Month") +
  geom_hline(yintercept = 1, col = "black", lty = "dashed", size = 0.2) +
  scale_fill_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  theme(text = element_text(size=8)) +
  guides(fill = "none") +
  ylab("Effective Reproductive Numer\nCases (1000s)")
national_Rt_epicurve

ggsave(plot = national_Rt_epicurve, filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/national_Rt_epicurve.pdf", width = 3, height = 3 )

#### GGvis by chiefdom ####
cholera_df_use %>%
  ggvis(x = ~Day, y = ~Reff_7dayMA, stroke = ~CHCODE) %>%
  filter(Reff >0) %>%
  filter(CHCODE %in% eval(input_radiobuttons(unique(cholera_df_use$CHCODE)))) %>%
  scale_datetime("x", domain = c(min(cholera_df_use_country$Day), max(cholera_df_use_country$Day)), nice = "month", label = "Day") %>%
  scale_datetime("y", domain = c(min(cholera_df_use$Reff_7dayMA, na.rm = TRUE), max(cholera_df_use$Reff_7dayMA, na.rm = TRUE)), label = "Reff") %>%
  layer_paths()

# #### Correlation between rainfall and Rt ####
# cholera_df_use$Rain_Chief <- apply(cholera_df_use, 1, function(x) fcn_lookup(x["CHCODE"], x["Day"], rainfall_TRMM_chiefdom[,c("CHCODE", "Date")], value_column = rainfall_TRMM_chiefdom$Daily_Rain))
# 
# plot(cholera_df_use$Reff, cholera_df_use$Rain_Chief)

#### Save workspace ####
save.image("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial_Cholera.RData")

#### Export cumulative data for arcmap ####
cholera_df_use_summary$id <- as.numeric(cholera_df_use_summary$id)
write.csv(cholera_df_use_summary, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/cholera_Reff_summary.csv", row.names = FALSE)

#####################
#Calculate variance #
#####################

# 
# Vij <- data.frame(quotientA4) / subset[,"Case"] #Now each row adds to 1
# Var <- rep(0, nrow(Vij))
# for (j in 1:nrow(Vij)){ #column
#   for (i in 1:length(Vij)){ #row
#     Var[j] = Var[j] + ( Vij[i,j]/subset[j,"Case"] * (1-Vij[j,i]/subset[j,"Case"]) * subset[i,"Case"] )
#   }
# }
# Var = Var/subset[,"Case"]
# #plot(data[,1], Var, type="l", xlab="Day", ylab="Variance")
# CI <- 1.96*sqrt(Var)
# 
# subset$AR10000 <- subset$Case / cumulative[cumulative$CHCODE == chief, "Population"] * 10000
# subset$Reff <- Reff
# subset$Reff_min <- subset$Reff - CI
# subset$Reff_max <- subset$Reff + CI
# subset[subset$Reff_min < 0, "Reff_min"] <- 0

