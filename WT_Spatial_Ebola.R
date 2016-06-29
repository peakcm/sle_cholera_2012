# WT_Spatial_Ebola

#### Load Libraries ####
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggvis)
library(forecast)
library(RColorBrewer)

#### Load workspace ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial_Ebola.RData")

#### Load Functions ####
source("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial.R")

#### Load ebola data ####
setwd('/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/')
cumulative <- read.csv('ebola_cumulative.csv')
cumulative <- cumulative[order(-cumulative$Ebola_Cases),]

df_daily <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_cholera_daily.csv")
df_daily$Day <- as.Date(as.character(as.POSIXct(df_daily$date, format = "%Y-%m-%d")))

df_daily$AR10000 <- NA
# for (chief in unique(ebola_df$CHCODE)){
#   ebola_df[ebola_df$CHCODE == chief, "AR10000"] <- ebola_df[ebola_df$CHCODE == chief, "Case"] / cumulative[cumulative$CHCODE == chief, "Population"] * 10000
# }

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
# Althaus 2015 Lancet ID
summary(rgamma(100000, rate=0.1697, shape = 2.5931))
tau <- (hist(rgamma(100000, rate=0.1697, shape = 2.5931), breaks = 0:1000)$counts)/100000
plot(tau)
plot(cumsum(tau), xlim = c(0,30), type = "l")

#### Select Data to use ####
# Chiefdoms with any cases
chiefs <- cumulative[order(-cumulative$Cumulative_Cases),"CHCODE"]

# Chiefdoms represented in the raw CDR dat
chiefs <- unique(weights[is.na(weights$CDR)==0,"CHCODE_source"])

# Chiefdoms with high case loads
chiefs <- cumulative[order(-cumulative$Cumulative_Cases),"CHCODE"][1:10]

#### Use WT_Spatial to calculate Reff Function ####
start <- min(df_daily[df_daily$confirmed_ebola > 0,"Day"])
end <- max(df_daily[df_daily$confirmed_ebola > 0,"Day"])
df_daily_ebola <- df_daily[df_daily$Day > start & df_daily$Day < end,]

ebola_df_use <- WT_Spatial(id.var = df_daily_ebola$CHCODE, date.var = df_daily_ebola$Day, case.var = df_daily_ebola$confirmed_ebola, weights = weights[,c("CHCODE_source", "CHCODE_destination", "inverse_distance_squared")], generation_interval = tau)

#### Summarize Reff calculations ####
ebola_df_use_summary <- WT_Spatial_Summary(ebola_df_use)

View(ebola_df_use_summary)

#### Moving Averages ####
names(ebola_df_use)[1:3] <- c("CHCODE", "Day", "Case")

# R-effective 7-day moving average
ebola_df_use = ebola_df_use %>% group_by(CHCODE) %>%
  mutate(Reff_7dayMA = ma(Reff, order = 7))

#### Add rainfall data ####
ebola_df_use$rain_avg <- NA
ebola_df_use$rain_avg_lag7 <- NA
ebola_df_use$area <- NA

for (chief in unique(ebola_df_use$CHCODE)){
  days <- data.frame(ebola_df_use)[ebola_df_use$CHCODE == chief,"Day"]
  rain <- df_daily[df_daily$CHCODE == chief & df_daily$Day %in% days,"rain_avg"]
  area <- df_daily[df_daily$CHCODE == chief & df_daily$Day %in% days,"area"]
  if (length(rain)==length(days)){
    ebola_df_use[ebola_df_use$CHCODE == chief, "rain_avg"] <- rain
    ebola_df_use[ebola_df_use$CHCODE == chief, "rain_avg_lag7"] <- c(rain[8:length(days)], rep(NA,7))
    ebola_df_use[ebola_df_use$CHCODE == chief, "area"] <- area
    cat(".")
  } else {cat("Error\n")}
}

summary(ebola_df_use$rain_avg)
summary(ebola_df_use$rain_avg_lag7)

ebola_df_use = ebola_df_use %>% group_by(CHCODE) %>%
  mutate(rain_7dayMA = ma(rain_avg, order = 7)) %>%
  mutate(rain_7dayMA_lag7 = ma(rain_7dayMA, order = 7))

ggplot(ebola_df_use, aes(x = Day, group = CHCODE, color = CHCODE)) +
  geom_point(aes(y = rain_avg)) +
  geom_line(aes(y = Reff_7dayMA))

ggplot(ebola_df_use, aes(x = rain_7dayMA_lag7, y = Reff_7dayMA)) +
  geom_point() 

#### Aggregate by Region or Nation ####
# Add region data
ebola_df_use$region <- NA
ebola_df_use[ebola_df_use$CHCODE < 2000,"region"] <- "East"
ebola_df_use[ebola_df_use$CHCODE > 2000 & ebola_df_use$CHCODE < 3000,"region"] <- "North"
ebola_df_use[ebola_df_use$CHCODE > 3000 & ebola_df_use$CHCODE < 4000,"region"] <- "South"
ebola_df_use[ebola_df_use$CHCODE > 4000,"region"] <- "West"

ebola_df_use$region <- factor(ebola_df_use$region, levels = c("West", "South", "North", "East"))

# Region average
ebola_df_use_region <- data.frame(cbind(Day = rep(as.Date(as.character(unique(ebola_df_use$Day))),4)))
ebola_df_use_region$region <- rep(c("East", "North", "South", "West"), each = length(unique(ebola_df_use$Day)))
ebola_df_use_region$Case <- NA
ebola_df_use_region$Reff_weighted <- NA
ebola_df_use_region$rain_weighted <- NA

for (i in 1:nrow(ebola_df_use_region)){
  day <- as.Date(ebola_df_use_region[i,"Day"])
  if (day == min(ebola_df_use_region$Day)){
    region <- ebola_df_use_region[i,"region"][1]
  }
  Reff <- unlist(ebola_df_use$Reff[ebola_df_use$Day %in% day & ebola_df_use$region %in% region])
  Case <- unlist(ebola_df_use$Case[ebola_df_use$Day %in% day & ebola_df_use$region %in% region])
  rain_avg <-unlist(ebola_df_use$rain_avg[ebola_df_use$Day %in% day & ebola_df_use$region %in% region]) 
  area <- unlist(ebola_df_use$area[ebola_df_use$Day %in% day & ebola_df_use$region %in% region]) 
  
  ebola_df_use_region[i, "Case"] <- sum(Case)
  if (sum(Case) > 0){
    ebola_df_use_region[i, "Reff_weighted"] <- weighted.mean(x = Reff , w = Case)
  } else {ebola_df_use_region[i, "Reff_weighted"] <- 0}
  ebola_df_use_region[i,"rain_weighted"] <- weighted.mean(x = rain_avg , w = area)
  
}

ebola_df_use_region$Reff_weighted_7dayMA <- ma(ebola_df_use_region$Reff_weighted, order = 7)
ebola_df_use_region$rain_weighted_7dayMA <- ma(ebola_df_use_region$rain_weighted, order = 7)

ebola_df_use_region$Day <- as.Date(ebola_df_use_region$Day)

ebola_df_use_region$region <- factor(ebola_df_use_region$region, levels = c("West", "South", "North", "East"))

ebola_df_use_region$Case_7dayMA <- ma(ebola_df_use_region$Case, order = 7)

# Country average
ebola_df_use_country <- data.frame(cbind(Day = as.Date(as.character(unique(ebola_df_use$Day)))))
ebola_df_use_country$Case <- NA
ebola_df_use_country$Reff_weighted <- NA
ebola_df_use_country$rain_weighted <- NA
ebola_df_use_country$Case_7dayMA <- NA

for (i in 1:nrow(ebola_df_use_country)){
  
  day <- as.Date(ebola_df_use_country[i,"Day"])
  
  Reff <- unlist(ebola_df_use$Reff[ebola_df_use$Day %in% day])
  Case <- unlist(ebola_df_use$Case[ebola_df_use$Day %in% day])
  rain <- unlist(ebola_df_use$rain_avg[ebola_df_use$Day %in% day])
  area <- unlist(ebola_df_use$area[ebola_df_use$Day %in% day])
  
  ebola_df_use_country[i, "Case"] <- sum(Case)
  if (sum(Case) > 0){
    ebola_df_use_country[i, "Reff_weighted"] <- weighted.mean(x = Reff , w = Case)
    ebola_df_use_country[i, "rain_weighted"] <- weighted.mean(x = rain , w = area)
  } else {ebola_df_use_country[i, "Reff_weighted"] <- 0}
}

ebola_df_use_country$Reff_weighted_7dayMA <- ma(ebola_df_use_country$Reff_weighted, order = 7)
ebola_df_use_country$rain_weighted_7dayMA <- ma(ebola_df_use_country$rain_weighted, order = 7)
ebola_df_use_country$Case_7dayMA <- ma(ebola_df_use_country$Case, order = 7)

ebola_df_use_country$Day <- as.Date(ebola_df_use_country$Day)

#### Days with Rt above 1 ####
View(ebola_df_use[ebola_df_use$Reff > 1,])
length(unlist(unique(ebola_df_use[ebola_df_use$Reff > 1,"CHCODE"])))
# 106 chiefdoms had Reff > 1 (with identity matrix)
# 37 chiefdoms had Reff > 1 (with inverse_distance_squared matrix)

ggplot(ebola_df_use[ebola_df_use$Reff > 1,], aes(x=Day)) + geom_bar() + ggtitle("Number of Chiefdoms with Reff > 1") 
# At the peak, 16 chiefdoms had Reff>1 (with identity matrix)
# At the peak, 15 chiefdoms had Reff>1 (with inverse_distance_squared matrix), At the left tail <5, right tail, <3 chiefdoms had Reff>1

ggplot(ebola_df_use[ebola_df_use$Case > 0,], aes(x=Day)) + geom_bar() + ggtitle("Number of Chiefdoms with Cases") 
# At the peak, between 30 and 44 chiefdoms had cases simultanously

table <- data.frame(table(ebola_df_use[ebola_df_use$Reff > 1,]$CHCODE))
ggplot(table, aes(x = factor(Var1), y = Freq)) + geom_bar(stat = "identity") + xlab("CHCODE") + ggtitle("Number of days with Reff > 1 by Chiefdom")

#### Store Data ####
# ebola_df_use_inverse_distance_squared <- ebola_df_use
# ebola_df_use_identity <- ebola_df_use

#### Select Data ####
# ebola_df_use <- ebola_df_use_inverse_distance_squared
# ebola_df_use <- ebola_df_use_identity

#### Chiefdom Specific Plot ####
ggplot() +
  theme_bw() +
  geom_bar(data = ebola_df_use, aes(x=Day, y = AR10000, group = factor(CHCODE), fill = factor(CHCODE)), stat = "identity", position = "dodge", alpha = 0.5, color = 0) +
  geom_line(data = ebola_df_use[ebola_df_use$Reff_7dayMA > 0 ,], aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), color = factor(CHCODE))) +
  geom_hline(yintercept = 1, col = "black", lty = "dashed") +
  xlab("Day") + 
  ylim(c(0, 3)) +
  # ylab("Reproductive Number R") + 
  ggtitle("Daily Reproductive Number") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

select_chiefs <- c("4299","4199","3410", "2406", "2511", "1212")
ggplot() +
  #   geom_bar(data = ebola_df_use[ebola_df_use$CHCODE %in% select_chiefs & ebola_df_use$Reff_7dayMA > 0 ,], aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), fill = factor(region)), stat = "identity", position = "dodge", width = 1, alpha = 1) +
  geom_line(data = ebola_df_use[ebola_df_use$Reff_7dayMA > 0.8 ,], aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), color = factor(region))) +
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
  
  geom_area(data = ebola_df_use_region, aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = .8, color = "black", size = 0.2) +
  geom_hline(yintercept = 1, col = "black", lty = "dashed", size = 0.3) +
  geom_line(data = ebola_df_use_region, aes(x = Day, y=2*rain_weighted_7dayMA/max(ebola_df_use_region$rain_weighted_7dayMA, na.rm = TRUE), group = factor(region)), size = .2, color = "blue", alpha = 0.3) +
  
  xlab("Day") + 
  facet_grid(region~.) +
  scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 2)) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b", name = "Month") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(R[t])) +
  scale_fill_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  theme(text = element_text(size=8)) +
  guides(fill = "none")
regional_Rt

ggsave(plot = regional_Rt, filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/Regional_Rt_Ebola.pdf", width = 3, height = 1.5 )

ggplot() +
  theme_bw() +
  geom_area(data = ebola_df_use_region[ebola_df_use_region$region == "North",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  geom_area(data = ebola_df_use_region[ebola_df_use_region$region == "East",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  geom_area(data = ebola_df_use_region[ebola_df_use_region$region == "West",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  geom_area(data = ebola_df_use_region[ebola_df_use_region$region == "South",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), fill = factor(region)), stat = "identity", alpha = 0.3) +
  
  geom_line(data = ebola_df_use_region[ebola_df_use_region$region == "North",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)), alpha = 0.8, size = 1.25) +
  geom_line(data = ebola_df_use_region[ebola_df_use_region$region == "East",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)), alpha = 0.8, size = 1.25) +
  geom_line(data = ebola_df_use_region[ebola_df_use_region$region == "West",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)),  alpha = 0.8, size = 1.25) +
  geom_line(data = ebola_df_use_region[ebola_df_use_region$region == "South",], aes(x = Day, y=Reff_weighted_7dayMA, group = factor(region), color = factor(region)), alpha = 0.8, size = 1.25) +
  
  # geom_line(data = ebola_df_use, aes(x = Day, y=Reff_7dayMA, group = factor(CHCODE), color = factor(region)), size = .3, color = "white", alpha = 0.4) +
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
  geom_bar(data = ebola_df_use_country, aes(x = Day, y = Case/sum(cumulative$Population)*100000), stat = "identity", alpha = 0.5) +
  geom_line(data = ebola_df_use_country[ebola_df_use_country$Reff_weighted_7dayMA > 0,], aes(x = Day, y = Reff_weighted_7dayMA)) +
  ylim(0, 7)+
  geom_hline(yintercept = 1, col = "black", lty = "dashed") +
  ylab("Effective Reproductive Numer (7-day Moving Average)\nCases per 100,000 Population")

ebola_df_use_region_weekly <- ebola_df_use_region %>% 
  group_by(region, Week = cut(Day, "week")) %>%
  mutate(Case_weekly = sum(Case)) %>%
  mutate(Case_weekly_7dayMA = sum(Case_7dayMA))

national_Rt_epicurve <- ggplot() +
  theme_bw() +
  geom_area(data = ebola_df_use_region_weekly, aes(x = Day, y = Case_weekly/100/2, fill = region), stat = "identity", alpha = 0.8) +
  geom_line(data = ebola_df_use_country[ebola_df_use_country$Reff_weighted_7dayMA > 0,], aes(x = Day, y = Reff_weighted_7dayMA), size = 0.3) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b", name = "Month") +
  geom_hline(yintercept = 1, col = "black", lty = "dashed", size = 0.2) +
  scale_fill_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  theme(text = element_text(size=8)) +
  ylim(0, 3) +
  guides(fill = "none") +
  ylab("Effective Reproductive Number")
national_Rt_epicurve

ggsave(plot = national_Rt_epicurve, filename ="/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/national_Rt_epicurve_ebola.pdf", width = 3, height = 3 )

#### GGvis by chiefdom ####
ebola_df_use %>%
  ggvis(x = ~Day, y = ~Reff_7dayMA, stroke = ~CHCODE) %>%
  filter(Reff >0) %>%
  filter(CHCODE %in% eval(input_radiobuttons(unique(ebola_df_use$CHCODE)))) %>%
  scale_datetime("x", domain = c(min(ebola_df_use_country$Day), max(ebola_df_use_country$Day)), nice = "month", label = "Day") %>%
  scale_datetime("y", domain = c(min(ebola_df_use$Reff_7dayMA, na.rm = TRUE), max(ebola_df_use$Reff_7dayMA, na.rm = TRUE)), label = "Reff") %>%
  layer_paths()

#### Save workspace ####
save.image("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial_Ebola.RData")

#### Export cumulative data for arcmap ####
ebola_df_use_summary$id <- as.numeric(ebola_df_use_summary$id)
write.csv(ebola_df_use_summary, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_Reff_summary.csv", row.names = FALSE)

