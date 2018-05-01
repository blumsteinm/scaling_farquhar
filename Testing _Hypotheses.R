
## Load necessary packages
require(plyr)
require(dplyr)

#### Now look at real climate_dataa to compare ####
dat <- read.csv("/Users/Meghs/Documents/GitHub/scaling_farquhar/Aggregated_Climate_Data.csv", header = T)
lai <- read.csv('/Users/Meghs/Documents/GitHub/scaling_farquhar/hf150-01-hem-lai.csv', header=TRUE)

## Assign variables new names and configurations for model
dat$Tv <- dat$Air_Temp_K
dat$Tv<-ave(dat$Tv, dat$month, dat$year, dat$day, dat$hour)
dat$relHum <- dat$Relative_Humidity_Percent
dat$relHum<-ave(dat$relHum, dat$month, dat$year, dat$day, dat$hour)
dat$Ca <- as.numeric(dat$Atmospheric_CO2)
dat$Ca <- ave(dat$Ca, dat$month, dat$year, dat$day, dat$hour)
dat$comp <- comp(Tv = dat$Tv)
dat$comp<-ave(dat$comp, dat$month, dat$year, dat$day, dat$hour)
dat$Vcmax <- Vcmax(Vmo = 35, Tv = dat$Tv, Tvlo=277.85)
dat$Vcmax<-ave(dat$Vcmax, dat$month, dat$year, dat$day, dat$hour)
dat$el <- el(Tv = dat$Tv)
dat$el<-ave(dat$el, dat$month, dat$year, dat$day, dat$hour)
dat$ea <- ea(relHum = dat$relHum, Tv = dat$Tv)
dat$ea<-ave(dat$ea, dat$month, dat$year, dat$day, dat$hour)
dat$PAR <- dat$Par_moles_m2_s * 1000000 ### only keep this until push updated to git
dat$PAR<-ave(dat$PAR, dat$month, dat$year, dat$day, dat$hour)
dat$time<-dat$hour

## Subset LAI and HF climate climate_dataa down to August  13th, 2013 during the day to test co2-limited funcitons
lai<-lai%>%filter(date=="2013-08-13")
lai$lai.sum<-ave(lai$lai.masked, FUN=sum) ## 54.28

dat.check<-dat%>%filter(year==2013)%>%filter(month==8)%>%filter(day==13)%>% ### choose this date due to LAI information
  dplyr::select(hour, Tv, relHum, Ca, comp, Vcmax, el, ea, PAR, time)
dat.check<-dat.check[!duplicated(dat.check),]

dat.check<-dat.check[(dat.check$time<24),]

dat.check$VPD<- dat.check$el-dat.check$ea
dat.check$Tair<-dat.check$Tv-273.15

FARAO(Ca=dat.check$Ca, VPD=dat.check$VPD, Tair=dat.check$Tair)


# climate_data_subset <- climate_data[climate_data$year == 2013 & climate_data$day == 13 & climate_data$month == 8,]
# climate_data_subset <- climate_data[climate_data$hour > 6 & climate_data$hour < 18,]


## Run model on 
farquhar_solver(input.df = dat.check, stomata = 'open') 


