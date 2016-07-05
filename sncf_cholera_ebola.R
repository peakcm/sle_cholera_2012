#Create a figure like S2 in WHO Ebola Response NEJM Supplement
# install.packages("ncf")
library(ncf)
library(ggplot2)

help(Sncf)

#### Pull strictly centroid geodata ####
# Create a dataset x with the longitude coordinates. Sort by chiefdom number
setwd("~/Documents/2014 Cholera OCV/Data - Raw/GIS Shapefiles/Shapefiles/Admin 3")
centroids <- read.csv("WAMerged_centroids_XY.csv")
# Create a dataset y with latitutde coordinates. Sort by chiefdom number
centroids <- centroids[ order(centroids[,"CHCODE"]), ]
x <- centroids$POINT_X
y <- centroids$POINT_Y

#### Load cumulative disease data ####
load("/Users/peakcm/Documents/2014 Cholera OCV/Data - Analysis/R codes/Data.RData")
head(df_cumulative)

df_sncf <- df_cumulative
df_sncf <- df_sncf[order(df_sncf$CHCODE),]

#### Add variables ####
df_sncf$cholera_any_case <- 0
df_sncf[df_sncf$cholera > 0,]$cholera_any_case <- 1

df_sncf$confirmed_ebola_any_case <- 0
df_sncf[df_sncf$confirmed_ebola > 0,]$confirmed_ebola_any_case <- 1

df_sncf$cholera_CAR <- df_sncf$cholera/df_sncf$Pop2012
df_sncf$confirmed_ebola_CAR <- df_sncf$confirmed_ebola/df_sncf$Pop2014

#### Function to choose x,y,z ####
fcn_choose_xyz <- function(x, y, df_z, column){
  rows_to_keep <- which(is.na(df_z[,column]) == 0)
  x <- x[rows_to_keep]
  y <- y[rows_to_keep]
  z <- df_z[rows_to_keep, column]
  
  return(list(x=x,y=y,z=z))
}

#### Function to run spline.correlog ####
fcn_run_spline.correlog <- function(x, y, df_z, column, resamp=1000, plot = FALSE){
  temp <- fcn_choose_xyz(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column)
  
  fit <- spline.correlog(temp$x,temp$y,temp$z, na.rm=TRUE, resamp=resamp, latlon=TRUE, df = 10, xmax = 150)
  if (plot){plot.spline.correlog(fit)}
  return(fit)
}

#### Function to transform spline.correlog output to useful ggplot format ####
fcn_transform_spline.correlog_output <- function(input){
  data.frame(x = input$real$predicted$x,
             y = input$real$predicted$y,
             ymin = input$boot$boot.summary$predicted$y["0.025",],
             ymax = input$boot$boot.summary$predicted$y["0.975",])
}

#### Function to plot ####
fcn_ggplot_spline.correlog <- function(fit1, fit2 = NA, title = NA, limits = NA, color = NA){
  fit1 <- fcn_transform_spline.correlog_output(fit1)
  
  plot <-
    ggplot() +
    
    theme_bw() +
    
    geom_hline(yintercept = 0, color = "black", lty = "dashed") +
    
    geom_ribbon(data = fit1, aes(x,ymin=ymin, ymax=ymax), fill = color[1], alpha = 0.2) +
    geom_line(data = fit1, aes(x,y), col = colors[1]) +
    
    ylab("Correlation") + xlab("Distance (km)") +
    ggtitle(title)
  
  if (is.na(fit2)[1] !=1){
    fit2 <- fcn_transform_spline.correlog_output(fit2)
    
    plot <- plot +  
      geom_ribbon(data = fit2, aes(x,ymin=ymin, ymax=ymax), fill = color[2], alpha = 0.2) +
      geom_line(data = fit2, aes(x,y), col = colors[2])
  }
  
  if (is.na(fit2)[1] !=1){ plot <- plot + coord_cartesian(ylim=limits)
  }
  
  return(plot)
}

#### Outcome: Day of first case ####
# Cholera
column = "cholera_onset"
fit1c <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

# Ebola
column = "confirmed_ebola_onset"
fit1e <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

#### Outcome: Cumulative cases ####
# Cholera
column = "cholera"
fit2c <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

# Ebola
column = "confirmed_ebola"
fit2e <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

#### Outcome: Cumulative Attack Rate ####
# Cholera
column = "cholera_CAR"
fit3c <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

# Ebola
column = "confirmed_ebola_CAR"
fit3e <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

#### Outcome: Any Case ####
# Cholera
column = "cholera_any_case"
fit4c <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

# Ebola
column = "confirmed_ebola_any_case"
fit4e <- fcn_run_spline.correlog(x = centroids$POINT_X, y = centroids$POINT_Y, df_z = df_sncf, column = column, resamp = 1000, plot = TRUE)

#### GGplot both diseases ####
limits = c(-1,1)
color = c("red", "blue")

fcn_ggplot_spline.correlog(fit1 = fit1c, fit2 = fit1e, title = "Outbreak Onset", limits = limits, color = color)
fcn_ggplot_spline.correlog(fit1 = fit2c, fit2 = fit2e, title = "Cumulative Cases", limits = limits, color = color)
fcn_ggplot_spline.correlog(fit1 = fit3c, fit2 = fit3e, title = "Cumulative Attack Rate", limits = limits, color = color)
fcn_ggplot_spline.correlog(fit1 = fit4c, fit2 = fit4e, title = "Disease Presence", limits = limits, color = color)

#### Experimental ####
# Try Sncf function
setwd("~/Documents/2014 Cholera OCV/Data - Analysis/Data Files/Admin 3/Admin 3 Weekly")
weekly <- read.csv("Weekly_Data.csv")
weekly <- weekly[order(weekly[,"CHCODE"]),]
weekly$CAR_per_10000
max(weekly$Day)
# Create a new dataset with each row a chiefdom and each column a week. Try with CAR  
z_CAR <- matrix(rep(NA, max(weekly$Day) * length(unique(weekly$CHCODE))), nrow = length(unique(weekly$CHCODE)))
row.names(z_CAR) <- sort(unique(weekly$CHCODE))
View(z_CAR)
for (i in 1:nrow(z_CAR)){
  chief <- as.numeric(row.names(z_CAR)[i])
  z_CAR[i,] <- weekly[weekly$CHCODE == chief,"CAR_per_10000"]
}
fit_5 <- Sncf(x,y, z_CAR, na.rm=TRUE, latlon=TRUE, resamp=1000, df = 10)
summary.Sncf(fit_5)
plot.Sncf(fit_5) #note this includes chiefdoms with 0 cases, so super highly correlated


#Try with just weekly attack rate
weekly$WAR <- weekly$Case / weekly$Population
z_WAR <- matrix(rep(NA, max(weekly$Day) * length(unique(weekly$CHCODE))), nrow = length(unique(weekly$CHCODE)))
row.names(z_WAR) <- sort(unique(weekly$CHCODE))
View(z_WAR)
for (i in 1:nrow(z_WAR)){
  chief <- as.numeric(row.names(z_WAR)[i])
  z_WAR[i,] <- weekly[weekly$CHCODE == chief,"WAR"]
}
z_WAR_nomiss <- z_WAR
for (i in  1: nrow(z_WAR_nomiss)){
  chief <- as.numeric(row.names(z_WAR_nomiss)[i])
  if (cdata[cdata$CHCODE == chief, "Any_Case"] == 0){
    z_WAR_nomiss[i,] <- NA
  }
}

fit_6 <- Sncf(x,y, z_WAR_nomiss, na.rm=TRUE, latlon=TRUE, resamp=1000, df = 10)
summary.Sncf(fit_6)
plot.Sncf(fit_6) #note this includes chiefdoms while zero cases were reported, so highly correlated. List as missing until first case?

#mark all "0" as NA until first case
z_WAR_nomiss_firstcase <- z_WAR_nomiss
for (i in 1:nrow(z_WAR_nomiss_firstcase)){
  chief <- as.numeric(row.names(z_WAR_nomiss_firstcase)[i])
  if (cdata[cdata$CHCODE == chief, "Any_Case"] == 1){
    counter = 0
    for (j in 1:ncol(z_WAR_nomiss_firstcase)){
      if (counter == 0){
        if (z_WAR_nomiss_firstcase[i,j] > 0){
          counter <- 1
        }
      }
      if (counter == 0){
        z_WAR_nomiss_firstcase[i,j] <- NA
      }
    }
  }
}
fit_7 <- Sncf(x,y, z_WAR_nomiss_firstcase, na.rm=TRUE, latlon=TRUE, resamp=10, df=10)
summary.Sncf(fit_7)
plot.Sncf(fit_7) # List as missing until first case. 

#mark all 0 as NA
z_WAR_nozero <- z_WAR
z_WAR_nozero[z_WAR_nozero == 0] <- NA
View(z_WAR_nozero)

fit_8 <- Sncf(x,y, z_WAR_nozero, na.rm=TRUE, latlon=TRUE, resamp=1000, df=10)
summary.Sncf(fit_8)
plot.Sncf(fit_8) # List all "0" as nA
