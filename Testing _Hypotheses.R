
## Load necessary packages
require(plyr)
require(dplyr)

#### Now look at real climate_dataa to compare ####
dat <- read.csv("/Users/Meghs/Documents/GitHub/scaling_farquhar/Aggregated_Climate_Data.csv", header = T)
lai <- read.csv('/Users/Meghs/Documents/GitHub/scaling_farquhar/hf150-01-hem-lai.csv', header=TRUE)
plantecophys <- read.csv("/Users/Meghs/Documents/GitHub/scaling_farquhar/plantecophys_output.csv")
flux <- read.csv("/Users/Meghs/Documents/GitHub/scaling_farquhar/hf103-03-flux-2004-2013.csv")

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

## Subset LAI, flux,  and HF climate climate_dataa down to August  13th, 2013 during the day to test co2-limited funcitons
lai<-lai%>%filter(date=="2013-08-13")
lai$lai.sum<-ave(lai$lai.masked, FUN=sum) ## 54.28

dat.check<-dat%>%filter(year==2013)%>%filter(month==8)%>%filter(day==13)%>% ### choose this date due to LAI information
  dplyr::select(hour, Tv, relHum, Ca, comp, Vcmax, el, ea, PAR, time)
dat.check<-dat.check[!duplicated(dat.check),]

dat.check<-dat.check[(dat.check$time<24),]

dat.check$VPD<- dat.check$el-dat.check$ea
dat.check$Tair<-dat.check$Tv-273.15

flux <- flux[flux$year == 2013 & trunc(flux$doy) == 225,]
flux$hem.fco2 <- flux$hem.fco2 * -1
hour <- sapply(as.character(flux$datetime), function(x) strsplit( strsplit(x, split = "T")[[1]][2], split = ":")[[1]][1] )
flux$hour <- hour
flux_subset <- ddply(flux, .(hour, year), summarize, flux_assimilation = mean(hem.fco2, na.rm = T))

## Run model on subsetted data
assimilation <- farquhar_solver(input.df = dat.check, stomata = 'open') 

## Combine information 
assimilation_ours <- cbind(assimilation, dat.check)
assimilation_plant <- cbind(plantecophys, dat.check)

jpeg("Results_1.jpeg")
par(family = "serif")
## Plot up initial results
plot(A.min ~ hour, data = assimilation_ours, pch = 16, family = "serif", col = "steelblue", ylab = "A.min", xlab = "Hour of Day", main = "Assimilation: August 13, 2013", type = "o", lwd = 2, ylim = c(0,12))

## Compare to plant ecophys
lines(ALEAF ~ hour, data = assimilation_plant, pch = 16, family = "serif", col = "forestgreen", type = 'o')

legend("topleft", c("Our Prediction", "Plantecophys Package"), col = c("steelblue", "forestgreen"), lwd = 2, bty = 'n')
dev.off()

jpeg("Results_2.jpeg")
## Scale up and see what it looks like & plot against flux data - how different is it
plot(assimilation_ours$A.min * mean(lai$lai.masked), pch = 16, col = "coral", type = "o", ylab = "Assimilation per m2 * Leaf Area Index", xlab = "Hour of Day", main = "August 13, 2013 Scaled, (Stomata Open)", xlim = c(-2, 26), ylim = c(-20, 80))
lines(flux_assimilation * 10 ~ as.numeric(flux_subset$hour), data = flux_subset, col = "forestgreen", lwd = 2)
legend("topleft", c("Our Estimates * LAI", "Flux microMoles/m2/seconds * 10"), lwd = 2, col = c("coral", "forestgreen"), bty = "n")
dev.off()


