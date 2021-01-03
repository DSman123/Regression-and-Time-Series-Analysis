# Joel Cabrera
# Regression and Time Series Analysis (16:954:596)
# Professor Robin
# December 5, 2020

### Preliminaries
library(TSA)
library(forecast)
library(tseries)
library(RColorBrewer)
library(ggplot2)
### 1. Loading & Checking NJ employment rate data
# Note: both data are already seasonally adjusted
# Loading data
data1 = read.csv("NJUR.csv") # NJ UR data
data2 = read.csv("LBSSA34.csv") # NJ LFPR data
# Checking UR data
dim(data1) # 538 x 2
head(data1)
summary(data1)
sum(is.na(data1)) # check for missing values
# Checking LFPR data
dim(data2)
head(data2)
summary(data2)
sum(is.na(data2)) # check for missing values
## 2a. Setting datasets as time series data
ur = ts(data1$NJUR, frequency = 12, start = c(1976, 1)) # monthly data, starting from January 1976
lfpr = ts(data2$LBSSA34, frequency = 12, start = c(1976, 1)) # monthly data, starting from January 1976
frequency(ur) # check if frequency = 12, instead of 1
frequency(lfpr)
## 2b. Plotting data
par(mfrow = c(2, 2))
# For UR
plot(ur, 
     xlab = "Year",
     ylab = "Monthly UR",
     main = "NJ Unemployment Rate (UR) from 1976-2020")
hist(ur, 
     xlab = "Monthly UR",
     ylab = "Frequency",
     main = "NJ Unemployment Rate (UR) from 1976-2020")
# For LFPR
plot(lfpr, 
     xlab = "Year",
     ylab = "Monthly LFPR",
     main = "NJ Labor Force Participation Rate (LFPR) from 1976-2020")
hist(lfpr, 
     xlab = "Monthly LFPR",
     ylab = "Frequency",
     main = "NJ Unemployment Rate (UR) from 1976-2020")
## 2c. Seasonal-means trend OLS
temp.color = c(rev(brewer.pal(6, 'RdYlBu')), brewer.pal(6, 'RdYlBu'))
par(mfrow = c(1, 2))
# For UR
month1. = season(ur)
ur.lm = lm(ur ~ month1.)
summary(ur.lm)
plot(y = rstudent(ur.lm), 
     x = as.vector(time(ur)),
     type = "l",
     ylab = "Standardized Residuals",
     main = "Residual Plot of Seasonal-Means Trend Model",
     xlab = "time")
points(y=rstudent(ur.lm),x=as.vector(time(ur)), col=temp.color, pch = 15)
hist(rstudent(ur.lm), xlab = "Residuals")
# ACF, runs, QQ-Plot, & Shapiro tests
acf(rstudent(ur.lm))
runs(rstudent(ur.lm))
qqnorm(rstudent(ur.lm))
qqline(rstudent(ur.lm))
shapiro.test(rstudent(ur.lm))
# For LFPR
month2. = season(lfpr)
lfpr.lm = lm(lfpr ~ month2.)
summary(lfpr.lm)
plot(y = rstudent(lfpr.lm), 
     x = as.vector(time(lfpr)),
     type = "l",
     ylab = "Standardized Residuals",
     main = "Residual Plot of Seasonal-Means Trend Model",
     xlab = "time")
points(y=rstudent(lfpr.lm),x=as.vector(time(lfpr)), col=temp.color, pch = 15)
hist(rstudent(lfpr.lm), xlab = "Residuals")
# ACF, runs, QQ-Plot, & Shapiro tests
acf(rstudent(lfpr.lm))
runs(rstudent(lfpr.lm))
qqnorm(rstudent(lfpr.lm))
qqline(rstudent(lfpr.lm))
shapiro.test(rstudent(lfpr.lm))
### 3. (More) Model Specification and Estimation
## 3a. Testing for stationarity
# ADF test: H0 hypothesis = observed time series is non-stationary (unit root) (HA hypothesis = stationary)
adf.test(ur) # p-value > 0.05, so non-stationary
adf.test(lfpr) # p-value > 0.05, so non-stationary
## 3b. ACFs and PACFs of Data
par(mfrow = c(2, 2))
# For UR
acf(diff(ur))
pacf(diff(ur))
# For LFPR
acf(diff(lfpr))
pacf(diff(lfpr))
## 3c. Allowing R to choose ARIMA orders/estimates
# For UR
ur_arima = auto.arima(ur)
ur_arima
# For LFPR
lfpr_arima = auto.arima(lfpr)
lfpr_arima
### 4. Model Diagnostics
## Residual Analysis
## For UR
ur_arima_res = rstandard(ur_arima) # different from ur_arima$residuals
# Residual plot & histogram
plot(ur_arima_res, 
     xlab = "Year",
     main = "Residual Plot",
     ylab = "ARIMA(0, 1, 0) Residuals")
abline(h = 0)
# QQ-plot & Shapiro test (for testing normality assumption)
qqnorm(ur_arima_res)
qqline(ur_arima_res)
shapiro.test(ur_arima_res)
# ACF, Ljung-Box, & runs tests
acf(ur_arima_res, main = "ARIMA(0, 1, 0) Model", ylab = "ACF of Residuals")
pacf(ur_arima_res, main = "ARIMA(0, 1, 0) Model", ylab = "PACF of Residuals")
LB.test(ur_arima) # model only, not residuals of model
runs(ur_arima_res)
## For LFPR
lfpr_arima_res = rstandard(lfpr_arima) # different from lfpr_arima$residuals
# Residual plot & histogram
plot(lfpr_arima_res, 
     xlab = "Year",
     main = "Residual Plot",
     ylab = "SARIMA(2, 2, 4)(0, 0, 1)[12] Residuals")
abline(h = 0)
# QQ-plot & Shapiro test (for testing normality assumption)
qqnorm(lfpr_arima_res)
qqline(lfpr_arima_res)
shapiro.test(lfpr_arima_res)
# ACF, Ljung-Box, & runs tests
acf(lfpr_arima_res, main = "SARIMA(2, 2, 4)(0, 0, 1)[12] Model", ylab = "ACF of Residuals")
pacf(lfpr_arima_res, main = "SARIMA(2, 2, 4)(0, 0, 1)[12] Model", ylab = "PACF of Residuals")
LB.test(lfpr_arima) # model only, not residuals of model
runs(lfpr_arima_res)
### 5. Model Forecasting
par(mfrow = c(1, 1))
## For UR
ur_arima_forecast = forecast(ur_arima, level = c(95), h = 12) # h = # of months ahead, 60 = 5 years; c(95) = 95% confidence interval
ur_arima_forecast # Note: refer to NJ news article for 8.2% increase
autoplot(ur_arima_forecast, 
         xlab = "Year",
         ylab = "Monthly LFPR",
         pch = 19)
## For LFPR
lfpr_arima_forecast = forecast(lfpr_arima, level = c(95), h = 12)
lfpr_arima_forecast
autoplot(lfpr_arima_forecast, 
         xlab = "Year",
         ylab = "Monthly LFPR",
         pch = 19)