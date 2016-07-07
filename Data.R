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
library(foreign)

#### Import raw cholera data ####
cholera_line_list <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Original Work/Epidemics/Data_Files/Cholera_Line_List.csv")

#### Import raw ebola data ####
ebola_daily_suspected <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Original Work/Epidemics/Data_Files/PNAS_Cases_Suspected.csv")

ebola_daily_confirmed <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Original Work/Epidemics/Data_Files/PNAS_Cases_Confirmed.csv")

#### Import population data ####
population <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Original Work/Epidemics/Data_Files/PNAS_Population_Adapted.csv")

#### Import chiefdom shapefiles ####
setwd("/Users/peakcm/Documents/2014 Cholera OCV/Original Work/Epidemics/Data_Files/Admin 3")
Admin3 <- readOGR(".","SLE_Adm3_1m_gov_WesternAreaMerged")
Admin3.map <- fortify(Admin3)
ggplot(Admin3.map, aes(x = long, y = lat, group=group)) + 
  geom_path() +
  coord_fixed()

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

#### Manipulate line list cholera data to daily and cumulative cholera data ####
cholera_line_list$date <- as.Date(cholera_line_list$Date_seen, format = "%m/%d/%y")
cholera_line_list$date_num <- as.numeric(cholera_line_list$date ) - min(as.numeric(cholera_line_list$date ))

cholera_line_list <- cholera_line_list[order(cholera_line_list$date, cholera_line_list$Temp_CHCODE),]

cholera_line_list <- cholera_line_list[,c("Age", "Age.group", "Sex", "Temp_CHCODE", "date", "date_num")]
names(cholera_line_list)[4]
names(cholera_line_list)[4] <- "CHCODE"

cholera_line_list_melt <- melt(cholera_line_list, id.vars = c("CHCODE", "date", "date_num"))

cholera_daily <- dcast(cholera_line_list_melt, formula = CHCODE + date + date_num~ variable, fun.aggregate = length)
cholera_daily <- cholera_daily[,c("CHCODE", "date", "date_num", "Age")]
names(cholera_daily) <- c("CHCODE", "date", "date_num", "cases")

cholera_cumulative <- dcast(cholera_line_list_melt, formula = CHCODE  ~ variable, fun.aggregate = length)
cholera_cumulative <- cholera_cumulative[,c("CHCODE", "Age")]
names(cholera_cumulative) <- c("CHCODE", "cases")

#### Manipulate line list ebola suspected data to daily and cumulative ####
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

#### Manipulate line list ebola confirmed data to daily and cumulative ####
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
# Confirmed
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

# Suspected
ebola_daily_suspected_dcast_cumulative$First_Case <- NA
ebola_daily_suspected_dcast_cumulative$Last_Case <- NA
ebola_daily_suspected_dcast_cumulative$t_dur <- NA

for (row in 1:nrow(ebola_daily_suspected_dcast_cumulative)){
  chief <- ebola_daily_suspected_dcast_cumulative[row, "CHCODE"]
  days <- ebola_daily_suspected[ebola_daily_suspected$CHCODE == chief, "date_num"]
  ebola_daily_suspected_dcast_cumulative[row, "First_Case"] <- days[1]
  ebola_daily_suspected_dcast_cumulative[row, "Last_Case"] <- days[length(days)]
  ebola_daily_suspected_dcast_cumulative[row, "t_dur"] <- ebola_daily_suspected_dcast_cumulative[row, "Last_Case"] - ebola_daily_suspected_dcast_cumulative[row, "First_Case"]
}

# Total
ebola_daily_total_dcast_cumulative$First_Case <- NA
ebola_daily_total_dcast_cumulative$Last_Case <- NA
ebola_daily_total_dcast_cumulative$t_dur <- NA

ebola_daily_total <- rbind(ebola_daily_confirmed, ebola_daily_suspected)

for (row in 1:nrow(ebola_daily_total_dcast_cumulative)){
  chief <- ebola_daily_total_dcast_cumulative[row, "CHCODE"]
  days <- ebola_daily_total[ebola_daily_total$CHCODE == chief, "date_num"]
  ebola_daily_total_dcast_cumulative[row, "First_Case"] <- days[1]
  ebola_daily_total_dcast_cumulative[row, "Last_Case"] <- days[length(days)]
  ebola_daily_total_dcast_cumulative[row, "t_dur"] <- ebola_daily_total_dcast_cumulative[row, "Last_Case"] - ebola_daily_total_dcast_cumulative[row, "First_Case"]
}


#### Make sure we have all the cholera cases ####
sum(cholera_daily$cases) == 22691

#### Manipulate cholera data ####
cholera_cumulative$First_Case <- NA
cholera_cumulative$Last_Case <- NA
cholera_cumulative$t_dur <- NA

for (row in 1:nrow(cholera_cumulative)){
  chief <- cholera_cumulative[row, "CHCODE"]
  days <- cholera_daily[cholera_daily$cases > 0 & cholera_daily$CHCODE == chief, "date_num"]
  cholera_cumulative[row, "First_Case"] <- days[1]
  cholera_cumulative[row, "Last_Case"] <- days[length(days)]
  cholera_cumulative[row, "t_dur"] <- cholera_cumulative[row, "Last_Case"] - cholera_cumulative[row, "First_Case"]
}

#### Combine daily datasets ####
CHCODEs <- population$CHCODE
date_start <- min(cholera_daily$date)
date_end <- max(ebola_daily_suspected$date)
dates <- seq(from = date_start, to = date_end, by = 1)
date_num <- seq(from = 0, to = length(dates)-1, by = 1)

df_daily <- data.frame("CHCODE" = rep(CHCODEs, each = length(dates)), "date" = rep(dates, times = length(CHCODEs)), "date_num" = rep(date_num, times = length(CHCODEs)))
df_daily$confirmed_ebola <- NA
df_daily$suspected_ebola <- NA
df_daily$cholera <- NA
df_daily$rain <- NA

for (i in unique(df_daily$date)){
 
  # Confirmed Ebola
  if (i %in% ebola_daily_confirmed_dcast_daily$date){
  df_daily[df_daily$date == i,"confirmed_ebola"] <- left_join(df_daily[df_daily$date == i,], 
                                                               ebola_daily_confirmed_dcast_daily[ebola_daily_confirmed_dcast_daily$date == i, c("CHCODE", "cases")], 
                                                               by = c("CHCODE"))$cases
  }
  #Suspected Ebola
  if (i %in% ebola_daily_suspected_dcast_daily$date){
    df_daily[df_daily$date == i,"suspected_ebola"] <- left_join(df_daily[df_daily$date == i,], 
                                                                ebola_daily_suspected_dcast_daily[ebola_daily_suspected_dcast_daily$date == i, c("CHCODE", "cases")], 
                                                                by = c("CHCODE"))$cases
  }
  #Cholera
  if (i %in% cholera_daily$date){
    df_daily[df_daily$date == i,"cholera"] <- left_join(df_daily[df_daily$date == i,], 
                                                        cholera_daily[cholera_daily$date == i, c("CHCODE", "cases")], 
                                                                by = c("CHCODE"))$cases
  }
  cat(".")
}

# Replacae NA with 0
df_daily[is.na(df_daily$confirmed_ebola),"confirmed_ebola"] <- 0
df_daily[is.na(df_daily$suspected_ebola),"suspected_ebola"] <- 0
df_daily[is.na(df_daily$cholera),"cholera"] <- 0

# Add "total_ebola" field
df_daily$total_ebola <- df_daily$suspected_ebola + df_daily$confirmed_ebola

# Check if we have th expected number of cases
sum(df_daily$confirmed_ebola) == 8358
sum(df_daily$suspected_ebola) == 3545
sum(df_daily$cholera) == 22691

# Add days since the start of each diseasea class
days_since_start_cholera <- seq(from = 0, to = length(dates)-1, by = 1)

days_since_start_ebola_suspected <- c(rep(0, which(dates == min(ebola_daily_suspected$date)) - 1), seq(0, length(dates) - which(dates == min(ebola_daily_suspected$date))))

days_since_start_ebola_confirmed <- c(rep(0, which(dates == min(ebola_daily_confirmed$date)) - 1), seq(0, length(dates) - which(dates == min(ebola_daily_confirmed$date))))

days_since_start_ebola_total <- c(rep(0, min(which(dates == min(ebola_daily_suspected$date)), which(dates == min(ebola_daily_confirmed$date))) - 1),
                                  seq(0, length(dates) - min(which(dates == min(ebola_daily_suspected$date)), which(dates == min(ebola_daily_confirmed$date)))))

# Check to see number of dates for each CHCODE
for (chief in unique(df_daily$CHCODE)){
  if (length(df_daily[df_daily$CHCODE == chief, "date"]) == 1345){cat(".")
    } else {cat(chief)}
}

df_daily <- cbind(df_daily, "days_since_start_cholera" = rep(days_since_start_cholera, length(unique(CHCODEs))), 
                  "days_since_start_ebola_suspected" = rep(days_since_start_ebola_suspected, length(unique(CHCODEs))),
                  "days_since_start_ebola_confirmed" = rep(days_since_start_ebola_confirmed, length(unique(CHCODEs))), 
                  "days_since_start_ebola_total" = rep(days_since_start_ebola_total, length(unique(CHCODEs))))

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
  cat(".")
}

#### Add Climate Prediction Center Rainfall data (New way) ####
setwd("/Users/peakcm/Documents/2014 Cholera OCV/Original Work/Epidemics/Data_Files/Climate Prediction Center/dbf")
files <- list.files(pattern = "dbf")

df_daily <- df_daily[order(df_daily$CHCODE, df_daily$date),]
df_daily$rain_avg <- NA
df_daily$area <- NA

for (file in files){
  date <- as.Date(strsplit(file, "\\.")[[1]][1], format = "%Y%m%d")
  if (date %in% df_daily$date){
    rain_data <- read.dbf(file)
    if (nrow(rain_data) == 151){
      if (sum(df_daily[df_daily$date == date,"CHCODE"] == rain_data$CHCODE)==151){
        df_daily[df_daily$date == date, "rain_avg"] <- rain_data$MEAN
        df_daily[df_daily$date == date, "area"] <- rain_data$AREA
        cat(".")
      } else{cat("Error 2 with", file, "\n")}
    } else{cat("Error 1 with", file, "\n")}
  } else{cat("x")}
}

ggplot(df_daily, aes(x = date, group = CHCODE, y = rain_avg, color = CHCODE)) + geom_point()

#### Confirm combined daily dataset is correct ####
chief <- 3101
plot(df_daily[df_daily$CHCODE == chief,]$date, df_daily[df_daily$CHCODE == chief,]$confirmed_ebola, type = "l")
points(ebola_daily_confirmed_dcast_daily[ebola_daily_confirmed_dcast_daily$CHCODE == chief,]$date, ebola_daily_confirmed_dcast_daily[ebola_daily_confirmed_dcast_daily$CHCODE == chief,]$cases)

chief <- 4299
plot(df_daily[df_daily$CHCODE == chief,]$date, df_daily[df_daily$CHCODE == chief,]$suspected_ebola, type = "l")
points(ebola_daily_suspected_dcast_daily[ebola_daily_suspected_dcast_daily$CHCODE == chief,]$date, ebola_daily_suspected_dcast_daily[ebola_daily_suspected_dcast_daily$CHCODE == chief,]$cases)

chief <- 3410
plot(df_daily[df_daily$CHCODE == chief,]$date, df_daily[df_daily$CHCODE == chief,]$cholera, type = "l")
points(cholera_daily[cholera_daily$CHCODE == chief,]$date, cholera_daily[cholera_daily$CHCODE == chief,]$cases)

#### Combine cumulative datasets ####
df_cumulative <- data.frame("CHCODE" = CHCODEs)
df_cumulative$Chiefdom <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population[,c("CHCODE")], value_column = population[,c("Chiefdom")], transformation = "as.character"))
df_cumulative$District <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population[,c("CHCODE")], value_column = population[,c("District")], transformation = "as.character"))
df_cumulative$Pop2014 <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population[,c("CHCODE")], value_column = population[,c("Total_2014")]))
df_cumulative$Pop2012 <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = population[,c("CHCODE")], value_column = population[,c("Total_2012")]))

df_cumulative$confirmed_ebola <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = ebola_daily_confirmed_dcast_cumulative[,c("CHCODE")], value_column = ebola_daily_confirmed_dcast_cumulative[,c("cases")]))
df_cumulative$suspected_ebola <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = ebola_daily_suspected_dcast_cumulative[,c("CHCODE")], value_column = ebola_daily_suspected_dcast_cumulative[,c("cases")]))
df_cumulative$total_ebola <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = ebola_daily_total_dcast_cumulative[,c("CHCODE")], value_column = ebola_daily_total_dcast_cumulative[,c("cases")]))
df_cumulative$cholera <- sapply(df_cumulative$CHCODE, function(x) fcn_lookup(query_1 = x, reference = cholera_cumulative[,c("CHCODE")], value_column = cholera_cumulative[,c("cases")]))

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
  cat(".")
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

hist(df_cumulative$total_ebola_onset, breaks = 20)
hist(df_cumulative$cholera_onset, breaks = 20)
plot(df_cumulative$cholera_onset, df_cumulative$total_ebola_onset)
ggplot(df_cumulative, aes(x = cholera_onset, y = total_ebola_onset, color = log(Pop2014+1))) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")
ggplot(df_cumulative, aes(x = cholera_onset, y = total_ebola_onset, color = log(Pop2014+1), shape = Region)) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")

ggplot(df_cumulative[df_cumulative$cholera_onset > 100,], aes(x = cholera_onset, y = total_ebola_onset, color = log(Pop2014+1))) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")
ggplot(df_cumulative[df_cumulative$cholera_onset > 100,], aes(x = cholera_onset, y = total_ebola_onset, color = log(Pop2014+1), shape = Region)) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")

hist(df_cumulative$total_ebola_duration, breaks = 20)
hist(df_cumulative$cholera_duration, breaks = 20)
plot(df_cumulative$cholera_duration, df_cumulative$total_ebola_duration)
ggplot(df_cumulative, aes(x = cholera_duration, y = total_ebola_duration, color = log(Pop2014+1))) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")
ggplot(df_cumulative, aes(x = cholera_duration, y = total_ebola_duration, color = log(Pop2014+1), shape = Region)) + geom_point(size = 4) + theme_classic() + stat_smooth(method = "lm")

summary(lm(formula = df_cumulative$total_ebola_duration ~ df_cumulative$cholera_duration + df_cumulative$Pop2014))
summary(lm(formula = df_cumulative$total_ebola_duration ~ df_cumulative$cholera_duration + df_cumulative$Pop2014 + factor(df_cumulative$Region)))

# Which chiefdoms were affected
cholera_affected <- df_cumulative[df_cumulative$cholera > 0,"CHCODE"]
ebola_affected <- df_cumulative[df_cumulative$total_ebola > 0,"CHCODE"]

length(cholera_affected)
length(ebola_affected)

# Risk ratio for being affected by ebola given that you were affected by cholera
a <- intersect(cholera_affected, ebola_affected) # Affected by both cholera and ebola
b <- setdiff(cholera_affected, ebola_affected) # Affected by cholera but not ebola
c <- setdiff(ebola_affected, cholera_affected) # Affected by ebola but not cholera
d <- setdiff(df_cumulative$CHCODE, union(cholera_affected, ebola_affected)) # Not affected by cholera or ebola
table <- epitable(length(a), length(b), length(c), length(d))
epitab(table, method = "riskratio")

cat("Chiefdoms affected by cholera were not any more likely to report ebola (RR = 0.83)")

#### Compare epi curves by chiefdom ####
ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$total_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = days_since_start_ebola_total, y = total_ebola)) +
  geom_line(aes(x = days_since_start_cholera, y = cholera), color = "red") +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$total_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = total_ebola), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$total_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera_cumulative), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = total_ebola_cumulative), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$total_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = days_since_start_cholera, y = cholera_cumulative_proportion), color = "red", alpha = 0.8) +
  geom_line(aes(x = days_since_start_ebola_total, y = total_ebola_cumulative_proportion), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$total_ebola > 10,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera_cumulative_proportion), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = total_ebola_cumulative_proportion), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

ggplot(df_daily[df_daily$CHCODE %in% df_cumulative[df_cumulative$total_ebola > 100,"CHCODE"],]) +
  geom_line(aes(x = date_num, y = cholera_cumulative_proportion), color = "red", alpha = 0.8) +
  geom_line(aes(x = date_num - 365*2, y = total_ebola_cumulative_proportion), alpha = 0.8) +
  xlim(0, 500) +
  facet_wrap(~CHCODE, scales = "free_y")

#### Compare epi curves by country ####
df_daily_country_confirmed_ebola <- df_daily %>% group_by(date) %>% summarize(confirmed_ebola = sum(confirmed_ebola)) 
df_daily_country_confirmed_ebola$confirmed_ebola_7dayMA <- ma(df_daily_country_confirmed_ebola$confirmed_ebola,order = 7 )

df_daily_country_suspected_ebola <- df_daily %>% group_by(date) %>% summarize(suspected_ebola = sum(suspected_ebola))
df_daily_country_suspected_ebola$suspected_ebola_7dayMA <- ma(df_daily_country_suspected_ebola$suspected_ebola,order = 7 )

df_daily_country_total_ebola <- df_daily %>% group_by(date) %>% summarize(total_ebola = sum(total_ebola))
df_daily_country_total_ebola$total_ebola_7dayMA <- ma(df_daily_country_total_ebola$total_ebola,order = 7 )

df_daily_country_cholera <- df_daily %>% group_by(date) %>% summarize(cholera = sum(cholera))
df_daily_country_cholera$cholera_7dayMA <- ma(df_daily_country_cholera$cholera,order = 7 )

# Confirmed ebola started on "2014-05-19". Freetown 38 days later on "2014-06-26"
# cholera on 2012-01-01, Freetown 174 days later on "2012-06-23"

ggplot() +
  theme_bw() +
  geom_bar(data = df_daily_country_cholera, aes(x = date, y = cholera), color = "pink", stat = "identity", alpha = 0.5) +
  geom_bar(data = df_daily_country_total_ebola, aes(x=date-365*2, y=total_ebola), color = "skyblue", stat = "identity", alpha = 0.5) +
  geom_line(data = df_daily_country_cholera, aes(x = date, y = cholera_7dayMA), color = "red") +
  geom_line(data = df_daily_country_total_ebola, aes(x=date-365*2, y=total_ebola_7dayMA), color = "blue") +
  scale_x_date(limits = c(as.Date("2012-01-01"), as.Date("2014-01-01")), name = "Month", date_breaks = "2 month", date_labels = "%b") +
  ylab("Daily Cases") +
  theme(text = element_text(size = 8)) +
  ggtitle("Cholera (red) and Total Ebola (blue) Cases\n{Ebola shifted by -2 years}")
  
ggsave(filename = "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/20160707_epicurves.pdf", height = 4, width = 4, units = "in")

#### Save workspace ####
save.image("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")

#### Load workspace ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")

#### Export data as csv files for ArcGIS ####
names(ebola_daily_total_dcast_cumulative)
names(ebola_daily_total_dcast_cumulative) <- c("CHCODE", "Ebola_Cases", "Last_Ebola_Case", "First_Ebola_Case", "t_dur_ebola")
write.csv(ebola_daily_total_dcast_cumulative, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/total_ebola_cumulative.csv")

cholera_cumulative_GIS <- cholera_cumulative[,c("CHCODE","cases", "First_Case", "Last_Case", "t_dur")]
names(cholera_cumulative_GIS) <- c("CHCODE", "Cholera_Cases", "Cholera_Onset", "Cholera_Last_Case", "Cholera_Duration")
write.csv(cholera_cumulative_GIS, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/cholera_cumulative_GIS.csv")

#### Export df_daily ####
write.csv(df_daily, "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_cholera_daily.csv")

#### Read df_daily ####
df_daily <- read.csv("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Data Files/ebola_cholera_daily.csv")

