#### Header ####
library(ggplot2)
library(reshape2)
library(epitools)
library(dplyr)
library(forecast)
require(sp)
require(rgdal)
require(maps)
library(fields)
library(geoR)
library(RColorBrewer)

#### Import cholera data ####
cholera_cumulative <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/Admin 3/Cumulative_Data_Chiefdom_New.csv")

cholera_weekly <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/Admin 3/Admin 3 Weekly/Weekly_Data.csv")

cholera_daily <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/Admin 3/Admin 3 Daily/Daily_Cases.csv")

#### Import ebola data ####
ebola_daily_suspected <- read.csv("/Users/peakcm/Dropbox/Ebola/Spatial Analysis SL PNAS/PNAS_Cases_Suspected.csv")

ebola_daily_confirmed <- read.csv("/Users/peakcm/Dropbox/Ebola/Spatial Analysis SL PNAS/PNAS_Cases_Confirmed.csv")

#### Import population data ####
population_2012 <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Raw/Population/SL_Chiefdom_Codes_Pop.csv")

population_2014 <- read.csv("/Users/peakcm/Dropbox/Ebola/Spatial Analysis SL PNAS/PNAS_Population.csv")

#### Import chiefdom shapefiles ####
setwd("/Users/peakcm/Documents/2014 Cholera OCV/Data - Raw/GIS Shapefiles/Shapefiles/Admin 3/")
Admin3 <- readOGR(".","SLE_Adm3_1m_gov_WesternAreaMerged")
Admin3.map <- fortify(Admin3)
ggplot(Admin3.map, aes(x = long, y = lat, group=group)) + 
  geom_path() +
  coord_fixed()

#### Import TRMM rainfall data ####
# Note that 4299 does NOT have a point within it

rainfall_TRMM <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Raw/Rainfall/TRMM Precipitation Data/TRMM_SL_2012_2013_Rainfall.csv")
names(rainfall_TRMM) <- c("Number", "Date", "Lat", "Long", "Rain")

# Wrangle into GIS format
rainfall_TRMM.points <- rainfall_TRMM
coordinates(rainfall_TRMM.points) <- c("Long", "Lat")
proj4string(rainfall_TRMM.points) <- proj4string(Admin3)

rainfall_TRMM.points$CHCODE <- over(rainfall_TRMM.points, Admin3)$CHCODE
rainfall_TRMM.points <- rainfall_TRMM.points[is.na(rainfall_TRMM.points$CHCODE)==0,]
View(rainfall_TRMM.points)

# Summarize
rainfall_TRMM <- as.data.frame(rainfall_TRMM.points)
rainfall_TRMM$Date <- as.Date(rainfall_TRMM$Date, format = "%d-%B-%y")
rainfall_TRMM$Region <- factor(substr(as.character(rainfall_TRMM$CHCODE),1,1), levels = c(4,3,2,1), labels = c("West", "South", "North", "East"))
rainfall_TRMM <- rainfall_TRMM[order(rainfall_TRMM$Date),]
  
rainfall_TRMM_chiefdom <- rainfall_TRMM %>% group_by(Date, CHCODE) %>%
  summarize(Daily_Rain = mean(Rain))
rainfall_TRMM_chiefdom <- rainfall_TRMM_chiefdom %>% group_by(CHCODE) %>%
  mutate(Daily_Rain_7dayMA = ma(Daily_Rain, order = 7))
rainfall_TRMM_chiefdom <- as.data.frame(rainfall_TRMM_chiefdom)

rainfall_TRMM_region <- rainfall_TRMM %>% group_by(Date, Region) %>%
  summarize(Daily_Rain = mean(Rain))
rainfall_TRMM_region <- rainfall_TRMM_region %>% group_by(Region) %>%
  mutate(Daily_Rain_7dayMA = ma(Daily_Rain, order = 7)) %>%
  mutate(Daily_Rain_14dayMA = ma(Daily_Rain, order = 14)) %>%
  mutate(Daily_Rain_28dayMA = ma(Daily_Rain, order = 28))
  
rainfall_TRMM_nation <- rainfall_TRMM %>% group_by(Date) %>%
  summarize(Daily_Rain = mean(Rain))
rainfall_TRMM_nation <- rainfall_TRMM_nation %>%
  mutate(Daily_Rain_7dayMA = ma(Daily_Rain, order = 7)) %>%
  mutate(Daily_Rain_14dayMA = ma(Daily_Rain, order = 14)) %>%
  mutate(Daily_Rain_28dayMA = ma(Daily_Rain, order = 28))

# Chiefdom
ggplot() +
  theme_bw() +
  geom_line(data = rainfall_TRMM_chiefdom[rainfall_TRMM_chiefdom$CHCODE == 4199,], aes(x = Date, y = Daily_Rain, group = CHCODE), color = "grey") +
  geom_line(data = rainfall_TRMM_chiefdom[rainfall_TRMM_chiefdom$CHCODE == 4199,], aes(x = Date, y = Daily_Rain_7dayMA, group = CHCODE)) +
  scale_x_date(limits = as.Date(c("2012-01-01", "2013-02-20")), date_breaks = "2 month", date_labels = "%b")

# Region
ggplot() +
  theme_bw() +
  geom_line(data = rainfall_TRMM_region, aes(x = Date, y = Daily_Rain_7dayMA, group = Region, color = Region, fill = Region)) +
  facet_grid(Region~.) +
  scale_fill_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  scale_color_brewer(type = "qual", palette = 8) + #qual palette 8 is good
  scale_x_date(limits = as.Date(c("2012-01-01", "2013-02-20")), date_breaks = "2 month", date_labels = "%b")

# Nation
ggplot() +
  theme_bw() +
  geom_line(data = rainfall_TRMM_nation, aes(x = Date, y = Daily_Rain_7dayMA)) +
  scale_x_date(limits = as.Date(c("2012-01-01", "2013-02-20")), date_breaks = "2 month", date_labels = "%b")

#### Helper functions ####
fcn_lookup <- function(query_1, query_2 = NA, reference, value_column, transformation = "none"){
  out <- 0
  if (is.na(query_2)==0){
    row_1 <- which(reference[,1] %in% query_1) 
    row_2 <- which(reference[,2] == query_2)
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

#### Manipulate ebola suspected data ####
ebola_daily_suspected$date <- as.Date(ebola_daily_suspected$Date.of.symptom.onset, format = "%d-%b-%y")
ebola_daily_suspected$date_num <- as.numeric(ebola_daily_suspected$date ) - min(as.numeric(ebola_daily_suspected$date ))

ebola_daily_suspected <- ebola_daily_suspected[,c("Age", "Sex", "District", "Chiefdom", "CHCODE", "date", "date_num")]

ebola_daily_suspected_melt <- melt(ebola_daily_suspected, id.vars = c("Chiefdom", "District", "CHCODE", "date", "date_num"))

ebola_daily_suspected_dcast_daily <- dcast(ebola_daily_suspected_melt, formula = CHCODE + Chiefdom + District + date + date_num~ variable, fun.aggregate = length)
ebola_daily_suspected_dcast_daily <- ebola_daily_suspected_dcast_daily[,c("CHCODE", "Chiefdom", "District", "date", "date_num", "Age")]
names(ebola_daily_suspected_dcast_daily) <- c("CHCODE", "Chiefdom", "District", "date", "date_num", "cases")

ebola_daily_suspected_dcast_cumulative <- dcast(ebola_daily_suspected_melt, formula = CHCODE + Chiefdom + District ~ variable, fun.aggregate = length)
ebola_daily_suspected_dcast_cumulative <- ebola_daily_suspected_dcast_cumulative[,c("CHCODE", "Chiefdom", "District", "Age")]
names(ebola_daily_suspected_dcast_cumulative) <- c("CHCODE", "Chiefdom", "District", "cases")

#### Manipulate ebola confirmed data ####
ebola_daily_confirmed$date <- as.Date(ebola_daily_confirmed$Date.of.symptom.onset, format = "%d-%b-%y")
ebola_daily_confirmed$date_num <- as.numeric(ebola_daily_confirmed$date ) - min(as.numeric(ebola_daily_confirmed$date ))

ebola_daily_confirmed <- ebola_daily_confirmed[,c("Age", "Sex", "District", "Chiefdom", "CHCODE", "date", "date_num")]

ebola_daily_confirmed_melt <- melt(ebola_daily_confirmed, id.vars = c("Chiefdom", "District", "CHCODE", "date", "date_num"))

ebola_daily_confirmed_dcast_daily <- dcast(ebola_daily_confirmed_melt, formula = CHCODE + Chiefdom + District + date + date_num~ variable, fun.aggregate = length)
ebola_daily_confirmed_dcast_daily <- ebola_daily_confirmed_dcast_daily[,c("CHCODE", "Chiefdom", "District", "date", "date_num", "Age")]
names(ebola_daily_confirmed_dcast_daily) <- c("CHCODE", "Chiefdom", "District", "date", "date_num", "cases")

ebola_daily_confirmed_dcast_cumulative <- dcast(ebola_daily_confirmed_melt, formula = CHCODE + Chiefdom + District ~ variable, fun.aggregate = length)
ebola_daily_confirmed_dcast_cumulative <- ebola_daily_confirmed_dcast_cumulative[,c("CHCODE", "Chiefdom", "District", "Age")]
names(ebola_daily_confirmed_dcast_cumulative) <- c("CHCODE", "Chiefdom", "District", "cases")

#### Create ebola cumulative total datasets ####
ebola_daily_total_dcast_cumulative <- data.frame("CHCODE" = unique(c(ebola_daily_suspected_dcast_cumulative$CHCODE, ebola_daily_confirmed_dcast_cumulative$CHCODE)))
ebola_daily_total_dcast_cumulative$cases <- 0
for (chief in ebola_daily_total_dcast_cumulative$CHCODE){
  suspected <- ebola_daily_suspected_dcast_cumulative[ebola_daily_suspected_dcast_cumulative$CHCODE == chief, "cases"]
  confirmed <- ebola_daily_confirmed_dcast_cumulative[ebola_daily_confirmed_dcast_cumulative$CHCODE == chief, "cases"]
  total = 0
  if (length(suspected) > 0){total = total + suspected}
  if (length(confirmed) > 0){total = total + confirmed}
  ebola_daily_total_dcast_cumulative[ebola_daily_total_dcast_cumulative$CHCODE == chief, "cases"] <- total
}

#### Manipulate ebola cumulative datasets ####
ebola_daily_confirmed_dcast_cumulative$First_Case <- NA
ebola_daily_confirmed_dcast_cumulative$Last_Case <- NA
ebola_daily_confirmed_dcast_cumulative$t_dur <- NA

for (row in 1:nrow(ebola_daily_confirmed_dcast_cumulative)){
    chief <- ebola_daily_confirmed_dcast_cumulative[row, "CHCODE"]
    days <- ebola_daily_confirmed[ebola_daily_confirmed$CHCODE == chief, "date_num"]
    ebola_daily_confirmed_dcast_cumulative[row, "First_Case"] <- days[1]
    ebola_daily_confirmed_dcast_cumulative[row, "Last_Case"] <- days[length(days)]
    ebola_daily_confirmed_dcast_cumulative[row, "t_dur"] <- ebola_daily_confirmed_dcast_cumulative[row, "Last_Case"] - ebola_daily_confirmed_dcast_cumulative[row, "First_Case"]
}

#### Manipulate cholera data ####
cholera_daily$date <- as.Date(cholera_daily$Day, format = "%m/%d/%y")

cholera_cumulative$onset_num <- as.numeric(as.Date(cholera_cumulative$First_Case, format = "%m/%d/%Y"))
cholera_cumulative$onset_num <- cholera_cumulative$onset_num - as.numeric(as.Date("1/1/2012", format = "%m/%d/%Y"))

cholera_cumulative$Last_Case <- NA
cholera_cumulative$tdur <- NA
for (row in 1:nrow(cholera_cumulative)){
  if (cholera_cumulative[row, "Any_Case"] == 1){
    chief <- cholera_cumulative[row, "CHCODE"]
    days <- cholera_daily[cholera_daily$CHCODE == chief & cholera_daily$Case > 0, "Day"]
    if (as.Date(cholera_cumulative[row, "First_Case"], format = "%m/%d/%Y") == as.Date(days[1], format = "%m/%d/%y")){
      cholera_cumulative[row, "Last_Case"] <- as.Date(days[length(days)], format = "%m/%d/%y")
      cholera_cumulative[row, "tdur"] <- cholera_cumulative[row, "Last_Case"]  - as.numeric(as.Date(cholera_cumulative[row, "First_Case"], format = "%m/%d/%Y"))
    } else {cat(".")}
  }
}

cholera_cumulative$Last_Case <- as.numeric(cholera_cumulative$Last_Case - min(cholera_cumulative[is.na(cholera_cumulative$Last_Case)==0, "Last_Case"]))

#### Combine daily datasets ####
CHCODEs <- unique(c(ebola_daily_suspected_dcast_cumulative$CHCODE, ebola_daily_confirmed_dcast_cumulative$CHCODE, cholera_daily$CHCODE))
date_start <- min(cholera_daily$date)
date_end <- max(ebola_daily_suspected$date)
dates <- seq(from = date_start, to = date_end, by = 1)
date_num <- seq(from = 0, to = length(dates)-1, by = 1)

df_daily <- data.frame("CHCODE" = rep(CHCODEs, each = length(dates)), "date" = rep(dates, times = length(CHCODEs)), "date_num" = rep(date_num, times = length(CHCODEs)))

# Add confirmed ebola data
df_daily$confirmed_ebola <- apply(df_daily[,c("CHCODE", "date")], MARGIN = 1, function(x) fcn_lookup(x[1], x[2], reference = ebola_daily_confirmed_dcast_daily[,c("CHCODE", "date")], ebola_daily_confirmed_dcast_daily$cases))

# Add suspected ebola data
df_daily$suspected_ebola <- apply(df_daily[,c("CHCODE", "date")], MARGIN = 1, function(x) fcn_lookup(x[1], x[2], reference = ebola_daily_suspected_dcast_daily[,c("CHCODE", "date")], ebola_daily_suspected_dcast_daily$cases))

# Add total ebola
df_daily$total_ebola <- df_daily$confirmed_ebola + df_daily$suspected_ebola

# Add cholera data
df_daily$cholera <- 0
cholera_daily_trim <- cholera_daily[cholera_daily$Case > 0,]
df_daily$cholera <- apply(df_daily[,c("CHCODE", "date")], MARGIN = 1, function(x) fcn_lookup(x[1], x[2], reference = cholera_daily_trim[,c("CHCODE", "date")], cholera_daily_trim$Case))

# Add days since the start of each diseasea class
days_since_start_cholera <- seq(from = 0, to = length(dates)-1, by = 1)

days_since_start_ebola_suspected <- c(rep(0, which(dates == min(ebola_daily_suspected$date)) - 1), seq(0, length(dates) - which(dates == min(ebola_daily_suspected$date))))

days_since_start_ebola_confirmed <- c(rep(0, which(dates == min(ebola_daily_confirmed$date)) - 1), seq(0, length(dates) - which(dates == min(ebola_daily_confirmed$date))))

days_since_start_ebola_total <- c(rep(0, min(which(dates == min(ebola_daily_suspected$date)), which(dates == min(ebola_daily_confirmed$date))) - 1),
                                  seq(0, length(dates) - min(which(dates == min(ebola_daily_suspected$date)), which(dates == min(ebola_daily_confirmed$date)))))

df_daily <- cbind(df_daily, "days_since_start_cholera" = rep(days_since_start_cholera, length(unique(CHCODEs))), "days_since_start_ebola_suspected" = rep(days_since_start_ebola_suspected, length(unique(CHCODEs))), "days_since_start_ebola_confirmed" = rep(days_since_start_ebola_confirmed, length(unique(CHCODEs))), "days_since_start_ebola_total" = rep(days_since_start_ebola_total, length(unique(CHCODEs))))

# Add rolling count of cases
df_daily$cholera_cumulative <- 0
df_daily$confirmed_ebola_cumulative <- 0
df_daily$suspected_ebola_cumulative <- 0
df_daily$total_ebola_cumulative <- 0

df_daily$cholera_cumulative_proportion <- 0
df_daily$confirmed_ebola_cumulative_proportion <- 0
df_daily$suspected_ebola_cumulative_proportion <- 0
df_daily$total_ebola_cumulative_proportion <- 0

for (chief in unique(df_daily$CHCODE)){
  df_daily[df_daily$CHCODE == chief, "cholera_cumulative"] <- cumsum(df_daily[df_daily$CHCODE == chief, "cholera"])
  df_daily[df_daily$CHCODE == chief, "confirmed_ebola_cumulative"] <- cumsum(df_daily[df_daily$CHCODE == chief, "confirmed_ebola"])
  df_daily[df_daily$CHCODE == chief, "suspected_ebola_cumulative"] <- cumsum(df_daily[df_daily$CHCODE == chief, "suspected_ebola"])
  df_daily[df_daily$CHCODE == chief, "total_ebola_cumulative"] <- cumsum(df_daily[df_daily$CHCODE == chief, "total_ebola"])
  
  if (max(df_daily[df_daily$CHCODE == chief, "cholera_cumulative"]) > 0){
    df_daily[df_daily$CHCODE == chief, "cholera_cumulative_proportion"] <- df_daily[df_daily$CHCODE == chief, "cholera_cumulative"]/max(df_daily[df_daily$CHCODE == chief, "cholera_cumulative"])
  }
  if (max(df_daily[df_daily$CHCODE == chief, "confirmed_ebola_cumulative"]) > 0){
    df_daily[df_daily$CHCODE == chief, "confirmed_ebola_cumulative_proportion"] <- df_daily[df_daily$CHCODE == chief, "confirmed_ebola_cumulative"]/max(df_daily[df_daily$CHCODE == chief, "confirmed_ebola_cumulative"])
  }
  if (max(df_daily[df_daily$CHCODE == chief, "suspected_ebola_cumulative"]) > 0){
    df_daily[df_daily$CHCODE == chief, "suspected_ebola_cumulative_proportion"] <- df_daily[df_daily$CHCODE == chief, "suspected_ebola_cumulative"]/max(df_daily[df_daily$CHCODE == chief, "suspected_ebola_cumulative"])
  }
  if (max(df_daily[df_daily$CHCODE == chief, "total_ebola_cumulative"]) > 0){
    df_daily[df_daily$CHCODE == chief, "total_ebola_cumulative_proportion"] <- df_daily[df_daily$CHCODE == chief, "total_ebola_cumulative"]/max(df_daily[df_daily$CHCODE == chief, "total_ebola_cumulative"])
  }
}

#### Confirm combined daily dataset is correct ####
chief <- 3101
plot(df_daily[df_daily$CHCODE == chief,]$date, df_daily[df_daily$CHCODE == chief,]$confirmed_ebola, type = "l")
points(ebola_daily_confirmed_dcast_daily[ebola_daily_confirmed_dcast_daily$CHCODE == chief,]$date, ebola_daily_confirmed_dcast_daily[ebola_daily_confirmed_dcast_daily$CHCODE == chief,]$cases)

chief <- 4299
plot(df_daily[df_daily$CHCODE == chief,]$date, df_daily[df_daily$CHCODE == chief,]$suspected_ebola, type = "l")
points(ebola_daily_suspected_dcast_daily[ebola_daily_suspected_dcast_daily$CHCODE == chief,]$date, ebola_daily_suspected_dcast_daily[ebola_daily_suspected_dcast_daily$CHCODE == chief,]$cases)

chief <- 3410
plot(df_daily[df_daily$CHCODE == chief,]$date, df_daily[df_daily$CHCODE == chief,]$cholera, type = "l")
points(cholera_daily[cholera_daily$CHCODE == chief,]$date, cholera_daily[cholera_daily$CHCODE == chief,]$Case)

#### Combine cumulative datasets ####
df_cumulative <- data.frame("CHCODE" = CHCODEs)
df_cumulative$Chiefdom <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population_2014[,c("CHCODE")], value_column = population_2014[,c("Chiefdom")], transformation = "as.character"))
df_cumulative$District <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population_2014[,c("CHCODE")], value_column = population_2014[,c("District")], transformation = "as.character"))
df_cumulative$Pop2014 <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population_2014[,c("CHCODE")], value_column = population_2014[,c("Total2014Inferred")]))
df_cumulative$Pop2012 <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population_2012[,c("CHCODE")], value_column = population_2012[,c("Pop_2012_est")]))

df_cumulative$confirmed_ebola <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = ebola_daily_confirmed_dcast_cumulative[,c("CHCODE")], value_column = ebola_daily_confirmed_dcast_cumulative[,c("cases")]))
df_cumulative$suspected_ebola <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = ebola_daily_suspected_dcast_cumulative[,c("CHCODE")], value_column = ebola_daily_suspected_dcast_cumulative[,c("cases")]))
df_cumulative$total_ebola <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = ebola_daily_total_dcast_cumulative[,c("CHCODE")], value_column = ebola_daily_total_dcast_cumulative[,c("cases")]))
df_cumulative$cholera <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = cholera_cumulative[,c("CHCODE")], value_column = cholera_cumulative[,c("Cumulative_Cases")]))

df_cumulative$confirmed_ebola_onset <- NA
df_cumulative$suspected_ebola_onset <- NA
df_cumulative$total_ebola_onset <- NA
df_cumulative$cholera_onset <- NA

df_cumulative$confirmed_ebola_duration <- NA
df_cumulative$suspected_ebola_duration <- NA
df_cumulative$total_ebola_duration <- NA
df_cumulative$cholera_duration <- NA

for (row in 1:nrow(df_cumulative)){
  chief <- df_cumulative[row, "CHCODE"]
  if (df_cumulative[row,"confirmed_ebola"] > 0){
    df_cumulative[row, "confirmed_ebola_onset"] <- min(df_daily[df_daily$CHCODE == chief & df_daily$confirmed_ebola, "days_since_start_ebola_confirmed"])
    df_cumulative[row, "confirmed_ebola_duration"] <- max(df_daily[df_daily$CHCODE == chief & df_daily$confirmed_ebola, "days_since_start_ebola_confirmed"]) - df_cumulative[row, "confirmed_ebola_onset"]
  }
  if (df_cumulative[row,"suspected_ebola"] > 0){
    df_cumulative[row, "suspected_ebola_onset"] <- min(df_daily[df_daily$CHCODE == chief & df_daily$suspected_ebola, "days_since_start_ebola_suspected"])
    df_cumulative[row, "suspected_ebola_duration"] <- max(df_daily[df_daily$CHCODE == chief & df_daily$suspected_ebola, "days_since_start_ebola_suspected"]) - df_cumulative[row, "suspected_ebola_onset"]
  }
  if (df_cumulative[row,"total_ebola"] > 0){
    df_cumulative[row, "total_ebola_onset"] <- min(df_daily[df_daily$CHCODE == chief & df_daily$total_ebola, "days_since_start_ebola_total"])
    df_cumulative[row, "total_ebola_duration"] <- max(df_daily[df_daily$CHCODE == chief & df_daily$total_ebola, "days_since_start_ebola_total"]) - df_cumulative[row, "total_ebola_onset"]
  }
  if (df_cumulative[row,"cholera"] > 0){
    df_cumulative[row, "cholera_onset"] <- min(df_daily[df_daily$CHCODE == chief & df_daily$cholera, "days_since_start_cholera"])
    df_cumulative[row, "cholera_duration"] <- max(df_daily[df_daily$CHCODE == chief & df_daily$cholera, "days_since_start_cholera"]) - df_cumulative[row, "cholera_onset"]
  }
}

df_cumulative$Region <- NA
df_cumulative[df_cumulative$CHCODE < 2000,"Region"] <- "East"
df_cumulative[df_cumulative$CHCODE > 2000 & df_cumulative$CHCODE < 3000,"Region"] <- "North"
df_cumulative[df_cumulative$CHCODE > 3000 & df_cumulative$CHCODE < 4000,"Region"] <- "South"
df_cumulative[df_cumulative$CHCODE > 4000 ,"Region"] <- "West"

#### Explore cumulative dataset ####
plot(df_cumulative$Pop2012, df_cumulative$Pop2014)

plot(df_cumulative$cholera, df_cumulative$total_ebola, log = "xy")
plot(df_cumulative$cholera, df_cumulative$suspected_ebola, log = "xy")
plot(df_cumulative$cholera, df_cumulative$confirmed_ebola, log = "xy")

plot(df_cumulative$cholera/df_cumulative$Pop2014, df_cumulative$total_ebola/df_cumulative$Pop2014, log = "xy")
plot(df_cumulative$cholera/df_cumulative$Pop2012, df_cumulative$total_ebola/df_cumulative$Pop2014, log = "xy")

hist(df_cumulative$confirmed_ebola_onset, breaks = 20)
hist(df_cumulative$cholera_onset, breaks = 20)
plot(df_cumulative$cholera_onset, df_cumulative$confirmed_ebola_onset)
ggplot(df_cumulative, aes(x = cholera_onset, y = confirmed_ebola_onset, color = log(Pop2014+1))) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")
ggplot(df_cumulative, aes(x = cholera_onset, y = confirmed_ebola_onset, color = log(Pop2014+1), shape = Region)) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")

ggplot(df_cumulative[df_cumulative$cholera_onset > 100,], aes(x = cholera_onset, y = confirmed_ebola_onset, color = log(Pop2014+1))) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")
ggplot(df_cumulative[df_cumulative$cholera_onset > 100,], aes(x = cholera_onset, y = confirmed_ebola_onset, color = log(Pop2014+1), shape = Region)) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")

hist(df_cumulative$confirmed_ebola_duration, breaks = 20)
hist(df_cumulative$cholera_duration, breaks = 20)
plot(df_cumulative$cholera_duration, df_cumulative$confirmed_ebola_duration)
ggplot(df_cumulative, aes(x = cholera_duration, y = confirmed_ebola_duration, color = log(Pop2014+1))) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")
ggplot(df_cumulative, aes(x = cholera_duration, y = confirmed_ebola_duration, color = log(Pop2014+1), shape = Region)) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")

summary(lm(formula = df_cumulative$confirmed_ebola_duration ~ df_cumulative$cholera_duration + df_cumulative$Pop2014))
summary(lm(formula = df_cumulative$confirmed_ebola_duration ~ df_cumulative$cholera_duration + df_cumulative$Pop2014 + factor(df_cumulative$Region)))

# Which chiefdoms were affected
cholera_affected <- df_cumulative[df_cumulative$cholera > 0,"CHCODE"]
ebola_affected <- df_cumulative[df_cumulative$confirmed_ebola > 0,"CHCODE"]

length(cholera_affected)
length(ebola_affected)

# Risk ratio for being affected by ebola given that you were affected by cholera
a <- intersect(cholera_affected, ebola_affected) # Affected by both cholera and ebola
b <- setdiff(cholera_affected, ebola_affected) # Affected by cholera but not ebola
c <- setdiff(ebola_affected, cholera_affected) # Affected by ebola but not cholera
d <- setdiff(df_cumulative$CHCODE, union(cholera_affected, ebola_affected)) # Not affected by cholera or ebola
table <- epitable(length(a), length(b), length(c), length(d))
epitab(table, method = "riskratio")

cat("Chiefdoms affected by cholera were not any more likely to report ebola (RR = 0.99)")

#### Compare epi curves by chiefdom ####
ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$confirmed_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = days_since_start_ebola_confirmed, y = confirmed_ebola)) +
  geom_line(aes(x = days_since_start_cholera, y = cholera), color = "red") +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$confirmed_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = confirmed_ebola), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$confirmed_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera_cumulative), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = confirmed_ebola_cumulative), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$confirmed_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = days_since_start_cholera, y = cholera_cumulative_proportion), color = "red", alpha = 0.8) +
  geom_line(aes(x = days_since_start_ebola_confirmed, y = confirmed_ebola_cumulative_proportion), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$confirmed_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera_cumulative_proportion), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = confirmed_ebola_cumulative_proportion), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$confirmed_ebola > 100,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera_cumulative_proportion), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = confirmed_ebola_cumulative_proportion), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

#### Compare epi curves by country ####
df_daily_country_confirmed_ebola <- df_daily %>% group_by(date) %>% summarize(confirmed_ebola = sum(confirmed_ebola)) 
df_daily_country_confirmed_ebola$confirmed_ebola_7dayMA <- ma(df_daily_country_confirmed_ebola$confirmed_ebola,order = 7 )

df_daily_country_suspected_ebola <- df_daily %>% group_by(date) %>% summarize(suspected_ebola = sum(suspected_ebola))
df_daily_country_suspected_ebola$suspected_ebola_7dayMA <- ma(df_daily_country_suspected_ebola$suspected_ebola,order = 7 )


df_daily_country_cholera <- df_daily %>% group_by(date) %>% summarize(cholera = sum(cholera))
df_daily_country_cholera$cholera_7dayMA <- ma(df_daily_country_cholera$cholera,order = 7 )

# Confirmed ebola started on "2014-05-19". Freetown 38 days later on "2014-06-26"
# cholera on 2012-01-01, Freetown 174 days later on "2012-06-23"

ggplot() +
  theme_bw() +
  geom_bar(data = df_daily_country_cholera, aes(x = date, y = cholera), color = "pink", stat = "identity", alpha = 0.5) +
  geom_bar(data = df_daily_country_confirmed_ebola, aes(x=date-365*2, y=confirmed_ebola), color = "skyblue", stat = "identity", alpha = 0.5) +
  geom_line(data = df_daily_country_cholera, aes(x = date, y = cholera_7dayMA), color = "red") +
  geom_line(data = df_daily_country_confirmed_ebola, aes(x=date-365*2, y=confirmed_ebola_7dayMA), color = "blue") +
  scale_x_date(limits = c(as.Date("2012-01-01"), as.Date("2014-01-01")), name = "Month", date_breaks = "1 month", date_labels = "%b") +
  ylab("Cases") +
  ggtitle("Cholera (red) and Confirmed Ebola (blue) cases by month\n{Ebola shifted by -2 years}")
  

#### Save workspace ####
save.image("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")

#### Load workspace ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")

#### Export data as csv files for ArcGIS ####
names(ebola_daily_confirmed_dcast_cumulative)
names(ebola_daily_confirmed_dcast_cumulative) <- c("CHCODE", "Chiefdom", "District", "Ebola_Cases", "Last_Ebola_Case", "First_Ebola_Case", "t_dur_ebola")
write.csv(ebola_daily_confirmed_dcast_cumulative, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_cumulative.csv")

cholera_cumulative_GIS <- cholera_cumulative[cholera_cumulative$Any_Case == 1,c("CHCODE","Cumulative_Cases", "onset_num", "Last_Case")]
cholera_cumulative_GIS$Last_Case <- as.numeric(as.Date(cholera_cumulative_GIS$Last_Case) - as.Date("2012-01-01"))
names(cholera_cumulative_GIS) <- c("CHCODE", "Cholera_Cases", "Cholera_Onset", "Cholera_Last_Case")
cholera_cumulative_GIS$t_dur_cholera <- cholera_cumulative_GIS$Cholera_Last_Case - cholera_cumulative_GIS$Cholera_Onset
write.csv(cholera_cumulative_GIS, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/cholera_cumulative_GIS.csv")

