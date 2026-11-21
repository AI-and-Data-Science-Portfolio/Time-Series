###############################################################################################
## EXPLORATORY ANALYSIS ##

# Import libraries
library(fpp2)
library(fpp3)
require(ggpubr)
require(tseries)
library(zoo)
library(ICSS)
library(forecast)

# Import data
project = read.csv("unemployment_ts.csv", header = T)
head(project)
str(project) 
class(project)
dim(project)
summary(project$VALUE)

# Create initial time series
ts1 = ts(project$VALUE, frequency = 12, start = c(1998, 1))

# Check for missing value
project[is.na(project$VALUE), ]

# Replace missing value using data from CSO website
missing = read.csv("MUM01.20260219102834.csv")[, c(4, 6, 8, 10)] %>%
  subset(Sex == "Both sexes" & Age.Group == "15 - 74 years") %>%
  .[, c(1, 4)]
missing = missing[c(14,146,162,179,182,234,245,253),]
project1 = project
project1[c(14,146,162,179,182,234,245,253),] = missing
table(complete.cases(project1))

project1[c(14,146,162,179,182,234,245,253),]

# Create time series with no missing values
ts2 = ts(project1$VALUE, frequency = 12, start = c(1998, 1))

# Check for outliers
outliers = tsoutliers(ts2)
outliers

# Replace outliers and create adjusted time series
ts3 = tsclean(ts2)
tsoutliers(ts3)
ts3

# Table of outliers and replacements 
months = project1$Month[outliers$index]
values1 = ts2[outliers$index]
values2 = ts3[outliers$index] 
outlier_details = data.frame( 
  Index = outliers$index,     
  Month = months,   
  Outlier = values1,
  Replacement = round(values2, digits = 1))
print(outlier_details)

# Summary statistics
summary(ts3)
format(as.yearmon(time(ts3)[which.max(ts3)]), "%b %Y")
format(as.yearmon(time(ts3)[which.min(ts3)]), "%b %Y")

# Plot adjusted time series vs. initial time series
plot1 = autoplot(ts3, series = "Adjusted") + 
  autolayer(ts1, series = "Original") + 
  xlab("Year") + ylab("Unemployment (000s)") + 
  ggtitle("Unemployment in Ireland 1998-2024") +
  theme(
    axis.title.x=element_text(angle = 0, color='black', face='bold', size = 12),
    axis.title.y=element_text(color='black', face='bold', size = 12),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.text.y = element_text(color = 'black'))
  scale_color_manual(values = c(`Adjusted` = "red", `Original` = "black"))
plot1

# Labour force data
library(readxl)
labour_force = read_excel("Labour Force 1998-2024.xlsx")
labour_force
ts(labour_force$LF, frequency = 4, start = c(1998, 1)) %>% autoplot() + xlab("Year") + ylab("Labour Force (000s)") +
  ggtitle("Labour Force 1998-2024")

# Time series, ACF and PACF plots
ggtsdisplay(ts3, main = "Time series, ACF and PACF plots of unemployment in Ireland 1998-2024")

# Lag plot
gglagplot(ts3) + xlab("Unemployment (000s)") + ylab("Unemployment (000s)") + 
  ggtitle("Lag plots of unemployment in Ireland 1998-2024") 

# Seasonal plot
ggseasonplot(ts3, year.labels = T, year.labels.left = T) + xlab("Month") + 
  ylab("Unemployment (000s)") + ggtitle("Seasonality of unemployment in Ireland 1998-2024")

###############################################################################################
### FORECASTING APPROACH ###

# Split data into training and test sets
ts_training = window(ts3, start=c(1998,1), end=c(2018,12))
ts_test = window(ts3, start = c(2019,1))

test_cases1 = length(ts3) - length(ts_training)

###################### ETS ###############################

# time series decomposition
comp3 = stl(ts3, s.window = "periodic")
autoplot(comp3) + xlab("Year") + ggtitle("Components of unemployment time series")

# ETS model
ETS = ets(ts_training, model="MAN", damped=T)
ETS
ETS_forecast = forecast(ETS, h = test_cases1)
summary(ETS_forecast)

# Holt model
Holt = holt(ts_training, initial = c("optimal"), h = test_cases1)
summary(Holt)
Holt

# damped Holt model
Holt_D = holt(ts_training, damped=T, h = test_cases1) 
summary(Holt_D)
Holt_D

# check accuracy of test set
round(accuracy(ETS_forecast, ts_test),2)
round(accuracy(Holt_D, ts_test),2) # best fit
round(accuracy(Holt, ts_test), 2)

# visualise point forecasts for three models 
ts_final = window(ts3, start = c(2019,1))
autoplot(ts_final, series = "Data", size=0.4) +
  autolayer(Holt_D$mean, series="Damped Holt", size=0.4) +
  autolayer(ETS_forecast$mean, series="Multiplicative", size=0.4) +
  autolayer(Holt$mean, series="Holt", size=0.4) +
  scale_colour_brewer(palette = "Dark2") + ylab("Unemployment") + 
  ggtitle("Point forecasts on test set") +
theme(
  axis.title.x=element_text(angle = 0, color='black', face='bold', size = 10),
  axis.title.y=element_text(color='black', face='bold', size = 10),
  axis.text.x = element_text(color = 'black', size = 10),
  axis.text.y = element_text(color = 'black'))

# check residuals
checkresiduals(Holt_D, lag = 12) 
checkresiduals(ETS_forecast, lag = 12)

# visualise with prediction intervals
autoplot(Holt_D) +
xlab("Year") + ylab("Unemployment (000s)") + 
  ggtitle("Damped Holt model - Forecast for test set") +
  theme(
    axis.title.x=element_text(angle = 0, color='black', face='bold', size = 12),
    axis.title.y=element_text(color='black', face='bold', size = 12),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.text.y = element_text(color = 'black'))
scale_color_manual(values = c(`Adjusted` = "red", `Original` = "black"))

###################### ARIMA ###############################

# check for constant variance
ICSS(ts_training, demean = FALSE)

# check time series is level or trend stationary
ggtsdisplay(ts_training, main = "Time series, ACF and PACF plots for training data") 
kpss.test(ts_training) # not stationary

# check number of difference required 
ndiffs(ts_training)
nsdiffs(ts_training)

# difference and second difference
ts_diff = diff(ts_training)
kpss.test(ts_diff)  # not stationary
adf.test(ts_diff)   # not stationary
ggtsdisplay(ts_diff, main = "Plot of time series with one level of differencing")  

ts_diff2 = diff(ts_diff)
kpss.test(ts_diff2)  # stationary
ggtsdisplay(ts_diff2, main = "Twice-differenced time series, ACF and PACF plots")

# auto ARIMA
ARIMA = auto.arima(ts_training, stepwise = F,  approximation= F, seasonal = F)
summary(ARIMA) # for example; write out actual model in project. This gives the constant form. Arima() gives the mean form

# forecast with ARIMA(2,2,3)
ARIMA_forecast = forecast(ARIMA, h = test_cases1)
autoplot(ARIMA_forecast)

# check accuracy and residuals
round(accuracy(ARIMA_forecast, ts_test), 2)
checkresiduals(ARIMA_forecast, lag =12)

# check other ARIMA(p,2,q) models
fit1 = Arima(ts_training, order=c(1,2,2))
ARIMA_forecast1 = forecast(fit1, h = test_cases1)
fit2 = Arima(ts_training, order=c(2,2,2))
ARIMA_forecast2 = forecast(fit2, h = test_cases1)
fit3 = Arima(ts_training, order=c(1,2,3))
ARIMA_forecast3 = forecast(fit3, h = test_cases1)
fit4 = Arima(ts_training, order=c(3,2,2))
ARIMA_forecast4 = forecast(fit4, h = test_cases1)
fit5 = Arima(ts_training, order=c(3,2,3))
ARIMA_forecast5 = forecast(fit5, h = test_cases1)
round(accuracy(ARIMA_forecast1, ts_test), 2)
round(accuracy(ARIMA_forecast2, ts_test), 2)
round(accuracy(ARIMA_forecast3, ts_test), 2)
round(accuracy(ARIMA_forecast4, ts_test), 2)
round(accuracy(ARIMA_forecast5, ts_test), 2)

# ARIMA with one level of differencing
ARIMA2 = auto.arima(ts_training, stepwise = F,  approximation= F, seasonal = F, max.d = 1, max.q = 3, max.p = 3)
ARIMA2

# forecast with ARIMA(1,1,3)
ARIMA_forecast2 = forecast(ARIMA2, h = test_cases1)
autoplot(ARIMA_forecast2)

# check accuracy and residuals
round(accuracy(ARIMA_forecast2, ts_test), 2)
checkresiduals(ARIMA_forecast2, lag = 12)

# visualise point forecasts for two models 
ts_final = window(ts3, start = c(2019,1))
autoplot(ts_final, series = "Data", size=0.4) +
  autolayer(Holt_D$mean, series="Damped Holt", size=0.4) +
  autolayer(ARIMA_forecast2$mean, series="ARIMA", size=0.4) +
  scale_colour_brewer(palette = "Dark2") + ylab("Unemployment") + 
  ggtitle("Model comparison on test set")

#################### CROSS VALIDATION ####################

# create cross validation objects
e1 = tsCV(ts3, holt, damped=T, h=12) # Damped Holt with additive errors
e1 

e2 = tsCV(ts3, function(x, h)        # Damped Holt with multiplicative errors
  {forecast(ets(x, model="MAN", damped=TRUE), h=h)}, h=12)
e2

e3 = tsCV(ts3, function(x, h)        # ARIMA(1,1,2)
{forecast(auto.arima(x, stepwise = F,  approximation= F, seasonal = F), h=h)}, h=12)
e3

e4 = tsCV(ts3, function(x, h)        # Random Walk with drift
  {rwf(x, drift = TRUE, h = h)}, h = 12)
e4

# compare MSE:
MSE1 = round(mean(e1^2, na.rm=T),2)
MSE2 = round(mean(e2^2, na.rm=T),2)
MSE3 = round(mean(e3^2, na.rm=T),2)
MSE4 = round(mean(e4^2, na.rm=T),2)

# compare MAE:
mae1 = round(mean(abs(e1), na.rm=T),2)
mae2 = round(mean(abs(e2), na.rm=T),2)
mae3 = round(mean(abs(e3), na.rm=T),2)
mae4 = round(mean(abs(e4), na.rm=T),2)

# fit naive method to the FULL time series to get in-sample errors
naive_fit = naive(ts3)
naive_insample_mae = mean(abs(naive_fit$residuals), na.rm = TRUE)

# calculate MASE for each model
mase1 = mae1 / naive_insample_mae
mase2 = mae2 / naive_insample_mae
mase3 = mae3 / naive_insample_mae
mase4 = mae4 / naive_insample_mae

# View results
results = data.frame(
  Model = c("Damped Holt", "Multiplicative", "ARIMA(1,1,2)", "Random Walk w/ Drift"),
  MSE = c(MSE1, MSE2, MSE3, MSE4),
  MAE = c(mae1, mae2, mae3, mae4),
  MASE = c(mase1, mase2, mase3, mase4))

print(results)

# check residuals
MULT = forecast(ets(ts3, model="ZZZ", damped=TRUE), h=12)
ARIMA112 = forecast(auto.arima(ts3, stepwise = F,  approximation= F, seasonal = F), h=12)
RWF = rwf(ts3, drift = TRUE, h = 12)
checkresiduals(MULT, lag = 12)
checkresiduals(ARIMA112, lag = 12)
checkresiduals(RWF, lag = 12)

###############################################################################################
## CONCLUSION ##

Holt_Final = holt(ts3, damped=T, h = 12) 
summary(Holt_Final)
Holt_Final$mean

# plot all methods together - add other methods later?
autoplot(ts3) +    
  autolayer(Holt_Final, PI=T, series = "Holt linear damped") +  
  ggtitle("Forecast Unemployment in Ireland 2025") +  
  xlab("Year") + ylab("Unemployment (000s)") +  
  guides(colour = guide_legend(title = "Forecast"))  

autoplot(ts_final) +    
  autolayer(Holt_Final, PI=T, series = "Holt linear damped") +  # add seasonal naive method for 11 periods
  ggtitle("Forecast Unemployment in Ireland 2025") +  
  xlab("Year") + ylab("Unemployment (000s)") +  
  guides(colour = guide_legend(title = "Forecast")) +
  theme(
    axis.title.x=element_text(angle = 0, color='black', face='bold', size = 12),
    axis.title.y=element_text(color='black', face='bold', size = 12),
    axis.text.x = element_text(color = 'black', size = 12),
    axis.text.y = element_text(color = 'black'))
scale_color_manual(values = c(`Adjusted` = "red", `Original` = "black"))

###############################################################################################
