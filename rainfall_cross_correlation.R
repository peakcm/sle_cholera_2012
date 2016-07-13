#### Header ####
# cross-correlation between case metrics and rain

#### Load this workspace ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/rainfall_cross_correlation.RData")

#### Save this workspace ####
# save.image("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/rainfall_cross_correlation.RData")

#### Load libraries ####
library(ggplot2)
library(dplyr)

#### Function ggplot colors ####
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color = ggplotColours(2)

#### Load Data workspace ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")
View(df_daily)

#### auto-correlation rainfall ####
lag.max = 100
acf.data.rain <- data.frame(lag = NA,
                       acf = NA,
                       CHCODE = NA)

for (chief in unique(df_daily$CHCODE)){
  acf <- acf(df_daily[df_daily$CHCODE == chief, "rain_avg"], lag.max = lag.max)
  acf.data.rain <- rbind(acf.data.rain, 
                    data.frame(lag = 1:(lag.max+1), 
                         acf = acf$acf,
                         CHCODE = rep(chief, lag.max+1)))
  cat(".")
}

ggplot(acf.data.rain, aes(x = jitter(lag), y = acf, group = CHCODE)) +
  geom_point(alpha = 0.2)

#### auto-correlation cholera cases ####
lag.max = 100
acf.data.cholera <- data.frame(lag = NA,
                       acf = NA,
                       CHCODE = NA)

for (chief in unique(df_daily$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "cholera"] > 0){
    first <- df_cumulative[df_cumulative$CHCODE == chief, "cholera_onset"] 
    last <- first + df_cumulative[df_cumulative$CHCODE == chief, "cholera_duration"] 
    acf <- acf(df_daily[df_daily$CHCODE == chief, "cholera"][first:last], lag.max = lag.max)
    acf.data.cholera <- rbind(acf.data.cholera, 
                              data.frame(lag = 1:length(acf$acf), 
                                         acf = acf$acf,
                                         CHCODE = rep(chief, length(acf$acf))))
    cat(".")
  }
}

ggplot(acf.data.cholera, aes(x = jitter(lag), y = acf, group = CHCODE, color = CHCODE)) + geom_line(alpha = 0.2)

#### auto-correlation ebola cases ####
lag.max = 100
acf.data.ebola <- data.frame(lag = NA,
                               acf = NA,
                               CHCODE = NA)

for (chief in unique(df_daily$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "total_ebola"] > 0){
    first <- df_cumulative[df_cumulative$CHCODE == chief, "total_ebola_onset"] + 863 # shift to first day of ebola
    last <- first + df_cumulative[df_cumulative$CHCODE == chief, "total_ebola_duration"] 
    acf <- acf(df_daily[df_daily$CHCODE == chief, "total_ebola"][first:last], lag.max = lag.max)
    acf.data.ebola <- rbind(acf.data.ebola, 
                              data.frame(lag = 1:length(acf$acf), 
                                         acf = acf$acf,
                                         CHCODE = rep(chief, length(acf$acf))))
    cat(".")
  }
}

ggplot(acf.data.ebola, aes(x = jitter(lag), y = acf, group = CHCODE, color = CHCODE)) + geom_line(alpha = 0.2)

#### cross-correlation rainfall-cholera ####
lag.max = 100
ccf.data.rain.cholera <- data.frame(lag = NA,
                               ccf = NA,
                               CHCODE = NA)

for (chief in unique(df_daily$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "cholera"] > 0){
    ccf <- ccf(df_daily[df_daily$CHCODE == chief, "rain_avg"], df_daily[df_daily$CHCODE == chief, "cholera"], lag.max = lag.max)
    ccf.data.rain.cholera <- rbind(ccf.data.rain.cholera, 
                              data.frame(lag = (-lag.max):(lag.max), 
                                         ccf = ccf$acf,
                                         CHCODE = rep(chief, 2*(lag.max)+1)))
    cat(".")
  }
}

ggplot(ccf.data.rain.cholera, aes(x = jitter(lag), y = ccf, group = CHCODE, color = CHCODE)) + geom_point(alpha = 0.2)

ccf.data.rain.cholera_summary <- ccf.data.rain.cholera %>% group_by(CHCODE) %>%
  summarize(max_lag = max(ccf, na.rm = TRUE))
ccf.data.rain.cholera_summary <- ccf.data.rain.cholera_summary[is.na(ccf.data.rain.cholera_summary$CHCODE)==0,]

ccf.data.rain.cholera_summary$max_lag_time <- NA
for (row in 1:nrow(ccf.data.rain.cholera_summary)){
  chief <- as.numeric(ccf.data.rain.cholera_summary[row, "CHCODE"])
  ccf.data.rain.cholera_summary[row, "max_lag_time"] <- ccf.data.rain.cholera[ccf.data.rain.cholera$CHCODE == chief & ccf.data.rain.cholera$ccf %in% ccf.data.rain.cholera_summary[row, "max_lag"], "lag"]
  cat(".")
}
ggplot(ccf.data.rain.cholera_summary, aes(x = max_lag_time)) + geom_bar(width = 7)

#### cross-correlation rainfall-ebola ####
lag.max = 100
ccf.data.rain.ebola <- data.frame(lag = NA,
                                    ccf = NA,
                                    CHCODE = NA)

for (chief in unique(df_daily$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "total_ebola"] > 0){
    ccf <- ccf(df_daily[df_daily$CHCODE == chief, "rain_avg"], df_daily[df_daily$CHCODE == chief, "total_ebola"], lag.max = lag.max)
    ccf.data.rain.ebola <- rbind(ccf.data.rain.ebola, 
                                   data.frame(lag = (-lag.max):(lag.max), 
                                              ccf = ccf$acf,
                                              CHCODE = rep(chief, 2*(lag.max)+1)))
    cat(".")
  }
}

ggplot(ccf.data.rain.ebola, aes(x = jitter(lag), y = ccf, group = CHCODE, color = CHCODE)) + geom_point(alpha = 0.2)

ccf.data.rain.ebola_summary <- ccf.data.rain.ebola %>% group_by(CHCODE) %>%
  summarize(max_lag = max(ccf, na.rm = TRUE))
ccf.data.rain.ebola_summary <- ccf.data.rain.ebola_summary[is.na(ccf.data.rain.ebola_summary$CHCODE)==0,]

ccf.data.rain.ebola_summary$max_lag_time <- NA
for (row in 1:nrow(ccf.data.rain.ebola_summary)){
  chief <- as.numeric(ccf.data.rain.ebola_summary[row, "CHCODE"])
  ccf.data.rain.ebola_summary[row, "max_lag_time"] <- ccf.data.rain.ebola[ccf.data.rain.ebola$CHCODE == chief & ccf.data.rain.ebola$ccf %in% ccf.data.rain.ebola_summary[row, "max_lag"], "lag"]
  cat(".")
}
ggplot(ccf.data.rain.ebola_summary, aes(x = max_lag_time)) + geom_bar(width = 7)

#### cross-correlation rainfall-cholera Reff ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial_Cholera.RData")

lag.max = 100
ccf.data.rain.choleraRe <- data.frame(lag = NA,
                                    ccf = NA,
                                    CHCODE = NA)

for (chief in unique(cholera_df_use$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "cholera"] > 0){
    ccf <- ccf(cholera_df_use[cholera_df_use$CHCODE == chief, "rain_avg"], cholera_df_use[cholera_df_use$CHCODE == chief, "Reff"], lag.max = lag.max)
    ccf.data.rain.choleraRe <- rbind(ccf.data.rain.choleraRe, 
                                   data.frame(lag = (-lag.max):(lag.max), 
                                              ccf = ccf$acf,
                                              CHCODE = rep(chief, 2*(lag.max)+1)))
    cat(".")
  }
}

ggplot(ccf.data.rain.choleraRe, aes(x = jitter(lag), y = ccf, group = CHCODE, color = CHCODE)) + geom_point(alpha = 0.2)

ccf.data.rain.choleraRe_summary <- ccf.data.rain.choleraRe %>% group_by(CHCODE) %>%
  summarize(max_lag = max(ccf, na.rm = TRUE))
ccf.data.rain.choleraRe_summary <- ccf.data.rain.choleraRe_summary[is.na(ccf.data.rain.choleraRe_summary$CHCODE)==0,]

ccf.data.rain.choleraRe_summary$max_lag_time <- NA
for (row in 1:nrow(ccf.data.rain.choleraRe_summary)){
  chief <- as.numeric(ccf.data.rain.choleraRe_summary[row, "CHCODE"])
  ccf.data.rain.choleraRe_summary[row, "max_lag_time"] <- ccf.data.rain.choleraRe[ccf.data.rain.choleraRe$CHCODE == chief & ccf.data.rain.choleraRe$ccf %in% ccf.data.rain.choleraRe_summary[row, "max_lag"], "lag"]
  cat(".")
}
ggplot(ccf.data.rain.choleraRe_summary, aes(x = max_lag_time)) + geom_bar(width = 7)

#### cross-correlation rainfall-ebola Reff ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/sle_cholera_2012/WT_Spatial_Ebola.RData")

lag.max = 100
ccf.data.rain.ebolaRe <- data.frame(lag = NA,
                                      ccf = NA,
                                      CHCODE = NA)

for (chief in unique(ebola_df_use$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "confirmed_ebola"] > 0){
    ccf <- ccf(ebola_df_use[ebola_df_use$CHCODE == chief, "rain_avg"], ebola_df_use[ebola_df_use$CHCODE == chief, "Reff"], lag.max = lag.max)
    ccf.data.rain.ebolaRe <- rbind(ccf.data.rain.ebolaRe, 
                                     data.frame(lag = (-lag.max):(lag.max), 
                                                ccf = ccf$acf,
                                                CHCODE = rep(chief, 2*(lag.max)+1)))
    cat(".")
  }
}

ggplot(ccf.data.rain.ebolaRe, aes(x = jitter(lag), y = ccf, group = CHCODE, color = CHCODE)) + geom_point(alpha = 0.2)

ccf.data.rain.ebolaRe_summary <- ccf.data.rain.ebolaRe %>% group_by(CHCODE) %>%
  summarize(max_lag = max(ccf, na.rm = TRUE))
ccf.data.rain.ebolaRe_summary <- ccf.data.rain.ebolaRe_summary[is.na(ccf.data.rain.ebolaRe_summary$CHCODE)==0,]

ccf.data.rain.ebolaRe_summary$max_lag_time <- NA
for (row in 1:nrow(ccf.data.rain.ebolaRe_summary)){
  chief <- as.numeric(ccf.data.rain.ebolaRe_summary[row, "CHCODE"])
  ccf.data.rain.ebolaRe_summary[row, "max_lag_time"] <- ccf.data.rain.ebolaRe[ccf.data.rain.ebolaRe$CHCODE == chief & ccf.data.rain.ebolaRe$ccf %in% ccf.data.rain.ebolaRe_summary[row, "max_lag"], "lag"]
  cat(".")
}
ggplot(ccf.data.rain.ebolaRe_summary, aes(x = max_lag_time)) + geom_bar(width = 7)

#### Cholera: Restructure data so that it's all in one time series with dead space in between ####
# http://stats.stackexchange.com/questions/23036/estimating-same-model-over-multiple-time-series
cholera.ts <- data.frame(CHCODE = NA, dates = NA, cholera = NA, rain = NA, cholera_Re = NA)
  
for (chief in unique(df_cumulative$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "cholera"] > 0){
    first <- max(1,df_cumulative[df_cumulative$CHCODE == chief, "cholera_onset"] )
    last <- first + df_cumulative[df_cumulative$CHCODE == chief, "cholera_duration"] 
    
    start <- max(1, first - 50)
    end <- last + 50
    
    dates <- df_daily[start:end, "date"]
    cholera <- c(rep(NA, times = first - start),
                 df_daily[df_daily$CHCODE == chief,"cholera"][first:last], 
                 rep(NA, times = end - last))
    cholera_Re <- c(rep(NA, times = first - start),
                 cholera_df_use[cholera_df_use$CHCODE == chief,]$Reff[(first+6):(last+6)],  #adjust by 7 because Reff goes to Jan 1
                 rep(NA, times = end - last))
    CHCODE <- rep(chief, length(dates))
    rain <- df_daily[df_daily$CHCODE == chief, "rain_avg"][start:end]
    
    cholera.ts <- rbind(cholera.ts, data.frame(CHCODE, dates, cholera, rain, cholera_Re))
    
  }
}

cholera.ts.ccf <- ccf(cholera.ts$cholera,cholera.ts$rain,  na.action = na.pass, lag.max = lag.max)
cholera_Re.ts.ccf <- ccf(cholera.ts$cholera_Re, cholera.ts$rain, na.action = na.pass, lag.max = lag.max)

cholera.ts.df <- data.frame(lag = cholera.ts.ccf$lag, 
                            cholera = cholera.ts.ccf$acf, 
                            cholera_Re = cholera_Re.ts.ccf$acf)
# http://stats.stackexchange.com/questions/211628/how-is-the-confidence-interval-calculated-for-the-acf-function
ci = 0.95
cholera.ts.df_ci <- qnorm((1 + ci)/2)/sqrt(cholera.ts.ccf$n.used)
cholera.ts.df_ci_Re <- qnorm((1 + ci)/2)/sqrt(cholera_Re.ts.ccf$n.used)

ggplot() +
  geom_hline(yintercept = cholera.ts.df_ci, lty = "dashed", col = "blue") +
  geom_hline(yintercept = -cholera.ts.df_ci, lty = "dashed", col = "blue") +
  geom_bar(data = cholera.ts.df, aes(x = lag, y = cholera),stat = "identity", fill = "red", alpha = 0.2) +
  geom_line(data = cholera.ts.df, aes(x = lag, y = ma(cholera,14)), color = "red", alpha = 0.5)

ggplot() +
  geom_hline(yintercept = cholera.ts.df_ci_Re, lty = "dashed", col = "blue") +
  geom_hline(yintercept = -cholera.ts.df_ci_Re, lty = "dashed", col = "blue") +
  geom_bar(data = cholera.ts.df, aes(x = lag, y = cholera_Re),stat = "identity", fill = "red", alpha = 0.2) +
  geom_line(data = cholera.ts.df, aes(x = lag, y = ma(cholera_Re,14)), color = "red", alpha = 0.5)

#### Ebola: Restructure data so that it's all in one time series with dead space in between ####
# http://stats.stackexchange.com/questions/23036/estimating-same-model-over-multiple-time-series
ebola.ts <- data.frame(CHCODE = NA, dates = NA, ebola = NA, rain = NA, ebola_Re = NA)

first_ebola <- 863 # shift to first day of ebola
# 2015-09-12 # last day
# 482 #days of ebola

for (chief in unique(df_cumulative$CHCODE)){
  if (df_cumulative[df_cumulative$CHCODE == chief, "total_ebola"] > 0){
    first <- max(1,df_cumulative[df_cumulative$CHCODE == chief, "total_ebola_onset"] )
    last <- first + df_cumulative[df_cumulative$CHCODE == chief, "total_ebola_duration"] 
    
    start <- first - 50
    end <- min(482, last + 50)
    
    dates <- df_daily[(first_ebola+start):(first_ebola+end), "date"]
    ebola <- c(rep(NA, times = first - start),
                 df_daily[df_daily$CHCODE == chief,"total_ebola"][(first_ebola+first):(first_ebola+last)], 
                 rep(NA, times = end - last))
    ebola_Re <- c(rep(NA, times = first - start),
                    ebola_df_use[ebola_df_use$CHCODE == chief,]$Reff[(first):(last)],  #adjust by 7 because Reff goes to Jan 1
                    rep(NA, times = end - last))
    CHCODE <- rep(chief, length(dates))
    rain <- df_daily[df_daily$CHCODE == chief, "rain_avg"][(first_ebola+start):(first_ebola+end)]
    
    ebola.ts <- rbind(ebola.ts, data.frame(CHCODE, dates, ebola, rain, ebola_Re))
    
  }
}

ebola.ts.ccf <- ccf( ebola.ts$ebola, ebola.ts$rain, na.action = na.pass, lag.max = lag.max)
ebola_Re.ts.ccf <- ccf(ebola.ts$ebola_Re, ebola.ts$rain,  na.action = na.pass, lag.max = lag.max)

ebola.ts.df <- data.frame(lag = ebola.ts.ccf$lag, 
                            ebola = ebola.ts.ccf$acf, 
                            ebola_Re = ebola_Re.ts.ccf$acf)
# http://stats.stackexchange.com/questions/211628/how-is-the-confidence-interval-calculated-for-the-acf-function
ci = 0.95
ebola.ts.df_ci <- qnorm((1 + ci)/2)/sqrt(ebola.ts.ccf$n.used)
ebola.ts.df_ci_Re <- qnorm((1 + ci)/2)/sqrt(ebola_Re.ts.ccf$n.used)

ggplot() +
  geom_hline(yintercept = ebola.ts.df_ci, lty = "dashed", col = "blue") +
  geom_hline(yintercept = -ebola.ts.df_ci, lty = "dashed", col = "blue") +
  geom_bar(data = ebola.ts.df, aes(x = lag, y = ebola),stat = "identity", fill = "blue", alpha = 0.2) +
  geom_line(data = ebola.ts.df, aes(x = lag, y = ma(ebola,14)), color = "blue", alpha = 0.5)

ggplot() +
  geom_hline(yintercept = ebola.ts.df_ci_Re, lty = "dashed", col = "blue") +
  geom_hline(yintercept = -ebola.ts.df_ci_Re, lty = "dashed", col = "blue") +
  geom_bar(data = ebola.ts.df, aes(x = lag, y = ebola_Re),stat = "identity", fill = "blue", alpha = 0.2) +
  geom_line(data = ebola.ts.df, aes(x = lag, y = ma(ebola_Re,14)), color = "blue", alpha = 0.5)

#### Plot both time series ####

ggplot() +
  theme_bw() +
  
  geom_vline(xintercept = 0, col = "darkgrey") +
  geom_hline(yintercept = 0, col = "darkgrey") +
  
  geom_hline(yintercept = cholera.ts.df_ci, lty = "dashed", col = color[1]) +
  geom_hline(yintercept = -cholera.ts.df_ci, lty = "dashed", col = color[1]) +
  
  geom_hline(yintercept = ebola.ts.df_ci, lty = "dashed", col = color[2]) +
  geom_hline(yintercept = -ebola.ts.df_ci, lty = "dashed", col = color[2]) +
  
  geom_bar(data = cholera.ts.df, aes(x = lag, y = cholera),stat = "identity", fill = color[1], alpha = 0.5) +
  geom_line(data = cholera.ts.df, aes(x = lag, y = ma(cholera,14)), color = color[1], alpha = 0.8) +
  
  geom_bar(data = ebola.ts.df, aes(x = lag, y = ebola),stat = "identity", fill = color[2], alpha = 0.5) +
  geom_line(data = ebola.ts.df, aes(x = lag, y = ma(ebola,14)), color = color[2], alpha = 0.8) +
  
  theme(text = element_text(size = 8)) +
  
  scale_x_continuous(name = "Lag (Weeks)", breaks = seq((7*-15),(7*15), by=14), labels = seq(-15, 15, by = 2)) +
  scale_y_continuous(name = "Cross-Correlation\nRainfall and Case Count", limits = c(-0.06, 0.11))

ggsave(filename = "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/20160712_ccf_rain_cases.pdf", height = 3, width = 5, units = "in")

ggplot() +
  theme_bw() +
  
  geom_vline(xintercept = 0, col = "darkgrey") +
  geom_hline(yintercept = 0, col = "darkgrey") +
  
  geom_hline(yintercept = cholera.ts.df_ci_Re, lty = "dashed", col = color[1]) +
  geom_hline(yintercept = -cholera.ts.df_ci_Re, lty = "dashed", col = color[1]) +
  
  geom_hline(yintercept = ebola.ts.df_ci_Re, lty = "dashed", col = color[2]) +
  geom_hline(yintercept = -ebola.ts.df_ci_Re, lty = "dashed", col = color[2]) +
  
  geom_bar(data = cholera.ts.df, aes(x = lag, y = cholera_Re),stat = "identity", fill = color[1], alpha = 0.5) +
  geom_line(data = cholera.ts.df, aes(x = lag, y = ma(cholera_Re,14)), color = color[1], alpha = 0.8) +
  
  geom_bar(data = ebola.ts.df, aes(x = lag, y = ebola_Re),stat = "identity", fill = color[2], alpha = 0.5) +
  geom_line(data = ebola.ts.df, aes(x = lag, y = ma(ebola_Re,14)), color = color[2], alpha = 0.8) +
  
  theme(text = element_text(size = 8)) +
  
  scale_x_continuous(name = "Lag (Weeks)", breaks = seq((7*-15),(7*15), by=14), labels = seq(-15, 15, by = 2)) +
  scale_y_continuous(name = "Cross-Correlation\nRainfall and Effective Reproductive Number", limits = c(-0.06, 0.11))

ggsave(filename = "/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/Figures/20160712_ccf_rain_Re.pdf", height = 3, width = 5, units = "in")

