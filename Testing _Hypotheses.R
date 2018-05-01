
## Load necessary packages
require(plyr)
require(dplyr)

#### Now look at real climate_dataa to compare ####
dat <- read.csv("/Users/Meghs/Documents/GitHub/scaling_farquhar/Aggregated_Climate_Data.csv", header = T)
lai <- read.csv('/Users/Meghs/Documents/GitHub/scaling_farquhar/hf150-01-hem-lai.csv', header=TRUE)
plantecophys <- read.csv("/Users/Meghs/Documents/GitHub/scaling_farquhar/plantecophys_output.csv")

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

## Run model on subsetted data
assimilation <- farquhar_solver(input.df = dat.check, stomata = 'open') 

## Combine information 
assimilation_ours <- cbind(assimilation, dat.check)
assimilation_plant <- cbind(plantecophys, dat.check)

## Plot up initial results
plot(A.min ~ hour, data = assimilation_ours, pch = 16, family = "serif", col = "steelblue", ylab = "A.min", xlab = "Hour of Day", main = "Assimilation: August 13, 2013", type = "o", lwd = 2, ylim = c(0,10))

## Compare to plant ecophys
lines(ALEAF ~ hour, data = assimilation_plant, pch = 16, family = "serif", col = "forestgreen", type = 'o')

## Scale up and see what it looks like & plot against flux data - how different is it
plot(assimilation_ours$A.min * mean(lai$lai.sum), pch = 16, col = "coral", type = "o", ylab = "Assimilation per leaf * Leaf Area Index", xlab = "Hour of Day", main = "August 13, 2013 Scaled")






