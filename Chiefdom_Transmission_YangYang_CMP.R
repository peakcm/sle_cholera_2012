library(gdata)
library(mgcv)

# weather.raw <- read.table("c:/research/ebola/data/weather/sl_weather_daily_20140101-20151101.txt", header=FALSE, skip=1, as.is=TRUE)
# colnames(weather.raw) <- c('station', 'WBAN', 'DATE', 'temp', 't1', 'dewp', 't2', 'SLP', 't3', 'STP', 't4', 'visib', 't5', 'wind', 't6', 'MXSPD', 'gust', 'max', 'min',  'prcp', 'SNDP', 'FRSHTT')
# weather <- weather.raw[,c('DATE', 'temp', 'dewp', 'wind')]
# weather$tmp <- (weather$temp - 32) * 5 / 9
# weather$dew <- (weather$dewp - 32) * 5 / 9
# a <- 17.625 ; b <- 243.04
# weather$rhx <- with(weather, 100*(exp((a*dew)/(b+dew) - (a*tmp)/(b+tmp))))
# weather$date <- as.Date(as.character(weather$DATE), '%Y%m%d')
# weather$temp<-NULL
# weather$dewp<-NULL
# weather$DATE<-NULL
# with(weather,cor(cbind(tmp, rhx, wind)))
# weather$week <- ceiling(as.numeric(weather$date - as.Date('20131229', '%Y%m%d'))/7)
# weekly.tmp <- aggregate(tmp~week, data=weather, FUN='mean')
# weekly.rhx <- aggregate(rhx~week, data=weather, FUN='mean')
# weekly.wind <- aggregate(wind~week, data=weather, FUN='mean')
# climate <- merge(weekly.tmp, weekly.rhx, by='week')
# climate <- merge(climate, weekly.wind, by='week', all.x=TRUE)
# 


cases.raw <- read.csv("/Users/peakcm/Dropbox/Ebola/Spatial Analysis SL PNAS/PNAS_Cases_Suspected.csv", header=TRUE, as.is=TRUE)

cases <- subset(cases.raw, select=c(code, district_name, chiefdom_name, poputotal, popudensity, week, 
                Local_case, NLag1_case, NLag2_case, NLag3_case, NLag4_case, LocalLag1_case, LocalLag2_case, LocalLag3_case, LocalLag4_case, 
                hospital_dis10, etc_dis10, primary_dis10, secondary_dis10, road_dis10,
                poorlevel, cropland10, forest10, shrub10, tribe_type))
colnames(cases)<-c('code', 'district', 'chiefdom', 'pop_size', 'pop_density', 'week', 
                'case', 'case_out1', 'case_out2', 'case_out3', 'case_out4', 'case_in1', 'case_in2', 'case_in3', 'case_in4', 
                'dis_hosp', 'dis_etc', 'dis_pri_rd', 'dis_sec_rd', 'dis_all_rd',
                'ecolevel', 'cropland', 'forest', 'shrub', 'tribe_type')
cases <- merge(cases, climate, by='week', all.x=TRUE)
save(cases, file="c:/research/ebola/data/R/panel_update.RData")                


load("c:/research/ebola/data/R/panel_update.RData")
#load("/home/yangyang/scratch/research/ebola/poisson_bootstrap/panel_update.RData")
#load("/home/yang/research/ebola/poisson/R_data/panel_update.RData")
cases$case_in  <- with(cases, 0.215 * case_in1  + 0.506 * case_in2  + 0.225 * case_in3  + 0.054 * case_in4)
cases$case_out <- with(cases, 0.215 * case_out1 + 0.506 * case_out2 + 0.225 * case_out3 + 0.054 * case_out4)  
cases$dis_clinic <- apply(subset(cases, select=c(dis_hosp, dis_etc)), 1, min)
cases$tribe4 <- cases$tribe3 <- cases$tribe2 <- cases$tribe1 <- 0 
cases$tribe1[cases$tribe_type == 2] <- 1 
cases$tribe2[cases$tribe_type == 3] <- 1 
cases$tribe3[cases$tribe_type == 4] <- 1 
cases$tribe4[cases$tribe_type == 5] <- 1 
# dichotomize poverty level, 0:extreme poverty <20%, 1:extreme poverty >= 20%
k<-which(cases$ecolevel <3)
if(length(k) > 0)  cases$ecolevel[k] <- 0
k<-which(cases$ecolevel >= 3)
if(length(k) > 0)  cases$ecolevel[k] <- 1

k<-which(cases$week < 21)
if(length(k) > 0)  cases$week[k] <- cases$week[k] + 52
cases <- cases[order(cases$code, cases$week),]
n.week <- length(unique(cases$week))
k<-which(cases$week > 22)
cases$tmp_lag[k] <- cases$tmp[k-2] 
cases$rhx_lag[k] <- cases$rhx[k-2] 
cases$wind_lag[k] <- cases$wind[k-2] 
k<-which(cases$shrub==0) 
if(length(k) > 0)  cases$shrub[k] <- 0.001
k<-which(is.na(cases$cropland)) 
if(length(k) > 0)  cases$cropland[k] <- 0.01
k<-which(cases$cropland>10) 
if(length(k) > 0)  cases$cropland[k] <- cases$cropland[k]/10
 
#cases$phase1 <- as.numeric(cases$week >= 38) 
cases$phase1 <- as.numeric(cases$week >= 40 & cases$week < 53) 
#cases$phase1 <- as.numeric(cases$week >= 40) 
cases$phase2 <- as.numeric(cases$week >= 53) 

cases$log_pop_den <- log(cases$pop_density)
#with(cases, cor(cbind(tmp, rhx, wind)))
#with(cases, cor(cbind(dis_all_rd, dis_hosp, dis_etc)))
#with(cases, cor(cbind(cropland, shrub, forest)))



cases.use <- subset(cases, week > 22)
cases.std <- cases.use
if(TRUE)
{
   x.mu <- list(tmp_lag=NA, rhx_lag=NA, wind_lag=NA, dis_pri_rd=NA, dis_sec_rd=NA, dis_all_rd=NA, 
                dis_hosp=NA, dis_etc=NA, dis_clinic=NA, cropland=NA, forest=NA, shrub=NA, log_pop_den=NA)
   x.sd <- list(tmp_lag=NA, rhx_lag=NA, wind_lag=NA, dis_pri_rd=NA, dis_sec_rd=NA, dis_all_rd=NA, 
                dis_hosp=NA, dis_etc=NA, dis_clinic=NA, cropland=NA, forest=NA, shrub=NA, log_pop_den=NA)
   x.range <- list(tmp_lag=c(NA, NA), rhx_lag=c(NA, NA), wind_lag=c(NA, NA), 
                   dis_pri_rd=c(NA, NA), dis_sec_rd=c(NA, NA), dis_all_rd=c(NA, NA), 
                   dis_hosp=c(NA, NA), dis_etc=c(NA, NA), dis_clinic=c(NA, NA), 
                   cropland=c(NA, NA), forest=c(NA, NA), shrub=c(NA, NA), log_pop_den=c(NA, NA))
                  
   x<- cases.use$tmp_lag
   x.mu$tmp_lag<-mean(x);   x.sd$tmp_lag<-sd(x);  x.range$tmp_lag <- (range(x) - x.mu$tmp_lag) / x.sd$tmp_lag
   x<- cases.use$rhx_lag
   x.mu$rhx_lag<-mean(x);   x.sd$rhx_lag<-sd(x);  x.range$rhx_lag <- (range(x) - x.mu$rhx_lag) / x.sd$rhx_lag
   x<- cases.use$wind_lag
   x.mu$wind_lag<-mean(x);  x.sd$wind_lag<-sd(x);  x.range$wind_lag <- (range(x) - x.mu$wind_lag) / x.sd$wind_lag
   x<- cases.use$dis_pri_rd
   x.mu$dis_pri_rd<-mean(x);  x.sd$dis_pri_rd<-sd(x);  x.range$dis_pri_rd <- (range(x) - x.mu$dis_pri_rd) / x.sd$dis_pri_rd
   x<- cases.use$dis_sec_rd
   x.mu$dis_sec_rd<-mean(x);  x.sd$dis_sec_rd<-sd(x);  x.range$dis_sec_rd <- (range(x) - x.mu$dis_sec_rd) / x.sd$dis_sec_rd
   x<- cases.use$dis_all_rd
   x.mu$dis_all_rd<-mean(x);  x.sd$dis_all_rd<-sd(x);  x.range$dis_all_rd <- (range(x) - x.mu$dis_all_rd) / x.sd$dis_all_rd
   x<- cases.use$dis_hosp
   x.mu$dis_hosp<-mean(x);  x.sd$dis_hosp<-sd(x);  x.range$dis_hosp <- (range(x) - x.mu$dis_hosp) / x.sd$dis_hosp
   x<- cases.use$dis_etc
   x.mu$dis_etc<-mean(x);   x.sd$dis_etc<-sd(x);  x.range$dis_etc <- (range(x) - x.mu$dis_etc) / x.sd$dis_etc
   x<- cases.use$dis_clinic
   x.mu$dis_clinic<-mean(x);   x.sd$dis_clinic<-sd(x);  x.range$dis_clinic <- (range(x) - x.mu$dis_clinic) / x.sd$dis_clinic
   x<- cases.use$cropland
   x.mu$cropland<-mean(x);   x.sd$cropland<-sd(x);  x.range$cropland <- (range(x) - x.mu$cropland) / x.sd$cropland
   x<- cases.use$forest
   x.mu$forest<-mean(x);   x.sd$forest<-sd(x);  x.range$forest <- (range(x) - x.mu$forest) / x.sd$forest
   x<- cases.use$shrub
   x.mu$shrub<-mean(x);   x.sd$shrub<-sd(x);  x.range$shrub <- (range(x) - x.mu$shrub) / x.sd$shrub
   x<- cases.use$log_pop_den
   x.mu$log_pop_den<-mean(x);   x.sd$log_pop_den<-sd(x);  x.range$log_pop_den <- (range(x) - x.mu$log_pop_den) / x.sd$log_pop_den
   
   cases.std$tmp_lag <- (cases.use$tmp_lag - x.mu$tmp_lag) / x.sd$tmp_lag
   cases.std$rhx_lag <- (cases.use$rhx_lag - x.mu$rhx_lag) / x.sd$rhx_lag
   cases.std$wind_lag <- (cases.use$wind_lag - x.mu$wind_lag) / x.sd$wind_lag
   cases.std$dis_pri_rd <- (cases.use$dis_pri_rd - x.mu$dis_pri_rd) / x.sd$dis_pri_rd
   cases.std$dis_sec_rd <- (cases.use$dis_sec_rd - x.mu$dis_sec_rd) / x.sd$dis_sec_rd
   cases.std$dis_all_rd <- (cases.use$dis_all_rd - x.mu$dis_all_rd) / x.sd$dis_all_rd
   cases.std$dis_hosp <- (cases.use$dis_hosp - x.mu$dis_hosp) / x.sd$dis_hosp
   cases.std$dis_etc <- (cases.use$dis_etc - x.mu$dis_etc) / x.sd$dis_etc
   cases.std$dis_clinic <- (cases.use$dis_clinic - x.mu$dis_clinic) / x.sd$dis_clinic
   cases.std$cropland <- (cases.use$cropland - x.mu$cropland) / x.sd$cropland
   cases.std$forest <- (cases.use$forest - x.mu$forest) / x.sd$forest
   cases.std$shrub <- (cases.use$shrub - x.mu$shrub) / x.sd$shrub
   cases.std$log_pop_den <- (cases.use$log_pop_den - x.mu$log_pop_den) / x.sd$log_pop_den
}

add.relative.humidity <- FALSE

if(add.relative.humidity)
{
   #final.model <- case ~ poly(tmp_lag, 2, raw=TRUE) + poly(rhx_lag, 2, raw=TRUE) +
   #                      poly(dis_pri_rd, 1, raw=TRUE) + poly(dis_etc, 2, raw=TRUE) + 
   #                      poly(cropland, 3, raw=TRUE)  + poly(forest, 3, raw=TRUE) + 
   #                      phase1 + phase2 + tribe1 + tribe2 + tribe3 + tribe4 +
   #                      offset(log(pop_size)) - 1
   final.model <- case ~ poly(tmp_lag, 2, raw=TRUE) + poly(rhx_lag, 2, raw=TRUE) +
                         poly(dis_etc, 2, raw=TRUE) + 
                         poly(cropland, 3, raw=TRUE)  + 
                         poly(log_pop_den, 1, raw=TRUE)  + 
                         phase1 + phase2 + tribe1 + tribe2 + tribe3 + tribe4 -1 
}  else
{
   final.model <- case ~ poly(tmp_lag, 2, raw=TRUE) + 
                         poly(dis_etc, 2, raw=TRUE) + 
                         poly(cropland, 3, raw=TRUE)  + 
                         poly(log_pop_den, 1, raw=TRUE)  + 
                         phase1 + phase2 + tribe1 + tribe2 + tribe3 + tribe4 -1 
}

# RR is the risk ratio estimated from GAM
log.lik <- function(pars, cases, RR)
{
    lambda0 <- exp(pars[1])
    lambda <- exp(pars[2])
    theta <- exp(pars[3])
    rate <- with(cases, (lambda0 + lambda * (case_in + theta * case_out)) * RR)
    sum(cases$case * log(rate) - rate)
}
       
lambda0 <- 0.01
lambda <- 0.5 
theta <- 0.05
for (i in 1:500)
{
   my.offset <- log(lambda0 + lambda*(cases.std$case_in + theta * cases.std$case_out))
   #fit.gam <- gam(final.model, offset=my.offset, method='REML', optimizer=c("outer","bfgs"), 
   #               control=list(maxit=200, epsilon = 1e-7,mgcv.tol=1e-7), family = poisson, data = cases.std)
   fit.gam <- gam(final.model, offset=my.offset, method='GCV.Cp', optimizer=c("outer","newton"),
                  control=list(maxit=200, epsilon = 1e-5,mgcv.tol=1e-5), family = poisson, data = cases.std)
   RR <- exp(fit.gam$linear.predictors - my.offset)
   ini <- log(c(lambda0, lambda, theta))
   est <- optim(ini, log.lik, cases=cases.std, RR= RR, method = "BFGS", control=list(fnscale=-1,maxit=200,trace=0),hessian = FALSE)
   lambda0 <- exp(est$par[1])
   lambda <- exp(est$par[2])
   theta <- exp(est$par[3])

   coeff<-coef(fit.gam)
   names(coeff) <- NULL
   print(c(round(i), lambda0, lambda, theta, coeff[1:5]))
}

if(add.relative.humidity)
{
   final.results.wk40.rhx <- list(lambda0=lambda0, lambda=lambda, theta=theta, fit=fit.gam, offset=my.offset)
   save(final.results.wk40.rhx, file="c:/research/ebola/data/R/fit_gam_wk40_rhx.RData")                
}  else
{
   final.results.wk40 <- list(lambda0=lambda0, lambda=lambda, theta=theta, fit=fit.gam, offset=my.offset)
   save(final.results.wk40, file="c:/research/ebola/data/R/fit_gam_wk40.RData")                
}

if(add.relative.humidity)
{
   load("c:/research/ebola/data/R/fit_gam_wk40_rhx.RData")
   final.results <- final.results.wk40.rhx
}  else
{
   load("c:/research/ebola/data/R/fit_gam_wk40.RData")
   final.results <- final.results.wk40
}
fit.gam <- final.results$fit

summary(final.results$fit)
exp(final.results$fit$coefficients)

