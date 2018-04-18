
## Load necessary packages
require(plyr)

## Read in Data
weather <- read.csv("/Users/Meghs/Downloads/hf001-10-15min-m.csv")
flux <- read.csv("/Users/Meghs/Downloads/hf103-03-flux-2004-2013.csv")

## Subset weather data into relevant parameters (PAR, AirTemp, & RH)
weather_subset <- weather[c("datetime", "jd", "airt", "parr", "rh")]

## Split out date for easier merging
dates <- weather_subset$datetime
year <- format(as.Date(dates), format = "%Y") ## Year
month <- format(as.Date(dates), format = "%m") ## month
day <- format(as.Date(dates), format = "%d") ## Day
hour <- sapply(as.character(dates), function(x) strsplit( strsplit(x, split = "T")[[1]][2], split = ":")[[1]][1] )
minute <- sapply(as.character(dates), function(x) strsplit( strsplit(x, split = "T")[[1]][2], split = ":")[[1]][2] )

## Add in Date/Time splits
weather_subset$year <- year
weather_subset$month <- month
weather_subset$day <- day
weather_subset$hour <- hour
weather_subset$minute <- minute

## Sumamrize to 30 minute intervals to match flux data; (convert Celcius to kelvins and Par micromoles/m2/s to moles/m2/s)
increment <- sapply(weather_subset$minute, function(x) ifelse(x == "30" | x == "00", 0, 1) )
weather_subset$summary_increment <- cumsum(increment)

weather <- ddply(weather_subset, .(summary_increment, year, month, day), summarize, 
           hour = hour[2],
           minute = minute[2],
           Julian_Day = mean(jd, na.rm = T), 
           Air_Temp_K = mean(airt, na.rm = T) + 273.15,
           Par_moles_m2_s = mean(parr, na.rm = T) * 1e-6, 
           Relative_Humidity_Percent = mean(rh, na.rm = T) 
           )

## merge with relevant flux parameters (co2)
(as.numeric(as.POSIXct(strptime(t, format = "%M:%OS"))) - 
    as.numeric(as.POSIXct(strptime("0", format = "%S"))))

# PAR (mol/m-2s-1) (micromolePerMeterSquaredPerSecond)
# Air temperature (K) (C)
# RH (unitless; percent, whichever) - preferred. If this is impossible, perhaps they report VPD instead? (%)
# atmospheric CO2 concentration (micromol/mol)













  