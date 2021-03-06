dummy <- data.frame(time = c(0:23))
dummy$Tv <- c(seq(280,300, length.out = 12), seq(300, 280, length.out = 12))
dummy$PAR <- c(rep(0, 5), seq(0, 2000, length.out = 7), seq(2000, 0, length.out = 7), rep(0, 5))
dummy$relHum <- c(rep(70, 3), 80, 80, 90, 100, 100, 90, 80, 70, 60, 50, 40, 50, 60, rep(70, 8))
dummy$relHum <- c(rep(70, 24))
dummy$Ca <- rep(400, 24)

### Dependent functions 
## Vcmax
Vcmax <- function (Vmo, Tvlo, Tv) {  #Vmo = 23, Tvlo = 277.15 are constants
   (Vmo * (exp(3000 * (1/288.15 -1/Tv)))) /
      ((1 + exp(.4 * (Tvlo - Tv))) * (1 + exp(0.4 * (Tv - 318.15))))
}

## compensation point
comp <- function(Tv) { 21.2 * exp(5000 * ((1/288.15) - (1/Tv))) } #units: umol/mol

## K1 and K2
K1 <- function(Tv) { 150 * exp(6000 * ((1/288.15) - (1/Tv))) } 
K2 <- function(Tv) {0.836 * exp(-1400 * ((1/288.15) - (1/Tv))) }

## el and ea
convertTv <- function(Tv) {
   TvC <- Tv -273.15
   return(TvC)
}
el <- function(Tv) {
   TvC <- convertTv(Tv)
   el <- 0.611 * exp((17.502 * TvC)/(TvC + 240.97))
   return(el)
}

ea <- function (relHum, Tv) {
   relHum /(100 * el(Tv))  # changed from multiplication to divide; Campbell and Norman page 42
}

input.df <- dummy
farquhar_solver <- function (input.df, stomata = c('open', 'closed')) { 
   ## input.df needs to have columns for Tv, PAR, Ca and relhum. the rest can be calculated... 
   
   ## defaults that we'd like to avoid including in the giant function call
   gamma <- 0.015
   Vmo <- 23 #92 #35
   b <- 2
   alpha <- 0.04
   Do <- 1.3 
   M <- 9  #stomatal slope; unitless; Raczka et al 2016, Biogeosciences; (also Heroult 2013, Plant Cell & Environment)
   Tvlo <- 277.85
   
   ## Functions used to calculate other inputs from input.df
   input.df$Vcmax <- Vcmax(Vmo = Vmo, Tvlo = Tvlo, Tv = input.df$Tv )
   input.df$comp <- comp(Tv = input.df$Tv)
   input.df$el <- el(Tv = input.df$Tv)
   input.df$ea <- ea(Tv = input.df$Tv, relHum = input.df$relHum)
   
   # calculated inputs
   input.df$Rd <- gamma * input.df$Vcmax
   input.df$k1.dat <- K1(input.df$Tv)
   input.df$k2.dat <- K2(input.df$Tv)
   
   # substitutions for simplification in quadratics
   input.df$X = input.df$k1.dat * (1 + input.df$k2.dat) # X = K1*(1+K2)
   input.df$Y <- alpha * input.df$PAR  # Y = alpha*PAR
   # F = (Ca - comp) * (1 + ((el - ea)/Do))
   input.df$FF <- ((input.df$Ca - input.df$comp) * (1 + ((input.df$el - input.df$ea)/Do)))  
   
   
   ### relevant functions
   # A1 = [(Vcmax * (Ci - comp)) / (Ci + K1 * (1+K2)) ] - Rd
   # A1 = [(Vcmax * (Ci - comp)) / (Ci + X) ] - Rd
   A1<- function(Ci, input.df) {
      A <- ((input.df$Vcmax * (Ci - input.df$comp)) / (Ci + input.df$X)) - input.df$Rd
      return(A)
   }

   # A2 = [ alpha * PAR * (Ci - comp) / (Ci + 2comp) ] - Rd
   A2 <- function(Ci, input.df) {
      A <- (alpha * input.df$PAR * ((Ci - input.df$comp) / (Ci + 2 * input.df$comp))) - input.df$Rd
      return(A)
   }
   # 
   
   # #Actually want to be solving for and minimizing J (post-quadratic) and minimizing, THEN subtracting R
   # 
   # J1 <- function(Ci, input.df) {
   #    J <- ((input.df$Vcmax * (Ci - input.df$comp)) / (Ci + input.df$X))
   #    return(J)
   # }
   # 
   # J2 <- function(Ci, input.df) {
   #    J <- (alpha * input.df$PAR * ((Ci - input.df$comp) / (Ci + 2 * input.df$comp))) 
   #    return(J)
   # }
   # 
   ### stomata closed ###
   if(stomata == 'closed') {
      ### CO2 limited case, coefficients for polyroot function
      aa <- (b * input.df$X * input.df$Ca/1.6) + (input.df$Vcmax * input.df$comp) + (input.df$Rd * input.df$X)
      bb <- (input.df$Rd) + (b * input.df$Ca / 1.6) - (b * input.df$X / 1.6) - (input.df$Vcmax)
      cc <- (-b / 1.6)
      
      z <- data.frame(aa = aa, bb = bb, cc = cc)  # where aa + bb*c1 + cc*c1^2
      #solve the polynomial
      roots <- apply(z, 1, polyroot)
      
      if(round(Im(roots[1]), 10) != 0) { 
         stop("quadratic roots are imaginary") 
      }
      
      #coerce into non-imaginary components
      roots.num <- Re(roots)
      #extract the non-negative value
      Ci.extract.A1 <- apply(roots.num, 2, max)
      
      # calculate A
      AA1 <- A1(Ci = Ci.extract.A1, input.df = input.df)
      
      ### Light limited case, coefficients for polyroot function
      aa <- ((b * 2 * input.df$comp * input.df$Ca / 1.6) + (input.df$Rd * 2 * input.df$comp) + (input.df$Y * input.df$comp))
      bb <- (input.df$Rd + (b * input.df$Ca / 1.6) - (b * 2 * input.df$comp / 1.6) - input.df$Y) 
      cc <- (-b / 1.6) 
      
      # define polynomial roots for each data point
      z <- data.frame(aa = aa, bb = bb, cc = cc)  # where aa + bb*c1 + cc*c1^2
      #solve the polynomial
      roots <- apply(z, 1, polyroot)
      
      if(round(Im(roots[1]), 10) != 0) { 
         stop("quadratic roots are imaginary") 
      }
      
      #coerce into non-imaginary components
      roots.num <- Re(roots)
      # extract the non-negative value
      Ci.extract.A2 <- apply(roots.num, 2, max)
      
      # calculate A2
      AA2 <- A2(Ci = Ci.extract.A2, input.df = input.df)  # only works if PAR has values 6 orders of magnitude higher
      
      
      ### Build output data frame
      # pick minimum for each time point: 
      A.df <- data.frame (AA1 = AA1, AA2 = AA2, Ci.A1 = Ci.extract.A1, Ci.A2 = Ci.extract.A2)
      A.df$A.min <- apply(A.df[,1:2], 1, min)
      A.df$min.eq <- apply(A.df[,1:2], 1, which.min)
      
      ## Solve for gsw ## 
      # stomata closed, gsw = b
      A.df$gsw <- rep(b, dim(A.df)[1])
   }
   
   
   
   #### Stomata Open ###
   if(stomata == 'open') {
      ############ vvvvv THIS IS WRONG!!!!!!! vvvvv ************
      ### CO2 limited coefficients
      aa <- ((b * input.df$X * input.df$Ca * input.df$FF) + (1.6 * input.df$FF * input.df$Vcmax * input.df$comp) - (1.6 * input.df$FF * input.df$Rd * input.df$X) - 
            (M * input.df$Vcmax * input.df$comp * input.df$Ca) + (M * input.df$Rd * input.df$X * input.df$Ca))
      bb <- ((-input.df$FF * b * input.df$Ca) - (b * input.df$X) - (1.6 * input.df$FF * input.df$Vcmax) + (1.6 * input.df$FF * input.df$Rd) + (M * input.df$Vcmax * input.df$Ca) +
            (M * input.df$Vcmax * input.df$comp) - (input.df$Rd * M * input.df$Ca) + (M * input.df$Rd * input.df$X))
      cc <- ((-b * input.df$FF) + (M * input.df$Vcmax) - (M * input.df$Rd))
      ############ ^^^^^ THIS IS WRONG!!!!!!! ^^^^^ ************
      
      # define polynomial roots for each data point
      z <- data.frame(aa = aa, bb = bb, cc = cc)  # where aa + bb*c1 + cc*c1^2
      #solve the polynomial
      roots <- apply(z, 1, polyroot)
      
      if(round(Im(roots[1]), 10) != 0) { 
         stop("quadratic roots are imaginary") 
      }
      
      #coerce into non-imaginary components
      roots.num <- Re(roots)
      #extract the non-negative value
      Ci.extract.A1 <- apply(roots.num, 2, max)
      
      #calculate A1
      AA1 <- A1(Ci = Ci.extract.A1, input.df = input.df)
      
      
      ### Light-limited coefficients
      aa <- ((b * input.df$FF * input.df$Ca * 2 * input.df$comp) - (1.6 * input.df$FF * input.df$Y * input.df$comp) + (1.6 * input.df$FF * input.df$Rd * 2 * input.df$comp) - 
            (M * input.df$Ca * input.df$Y * input.df$comp) + (M * input.df$Rd * 2 * input.df$comp * input.df$Ca))
      bb <- ((b * input.df$FF * input.df$Ca) - (b * input.df$FF * 2 * input.df$comp) - (1.6 * input.df$FF * input.df$Y) + (1.6 * input.df$FF * input.df$Rd) + (M * input.df$Y * input.df$Ca) +
            (M * input.df$Y * input.df$comp) - (M * input.df$Rd * input.df$Ca) - (M * input.df$Rd * 2 * input.df$comp))
      cc <- ((-b * input.df$FF) - (M * input.df$Y) + (M * input.df$Rd))
      
      # define polynomial roots for each data point
      z <- data.frame(aa = aa, bb = bb, cc = cc)  # where aa + bb*c1 + cc*c1^2
      #solve the polynomial
      roots <- apply(z, 1, polyroot)
      
      if(round(Im(roots[1]), 10) != 0) { 
         stop("quadratic roots are imaginary") 
      }
      
      # coerce into non-imaginary components
      roots.num <- Re(roots) 
      # extract the non-negative value
      Ci.extract.A2 <- apply(roots.num, 2, max)
      
      # calculate A2
      AA2 <- A2(Ci = Ci.extract.A2, input.df = input.df)  # only works if PAR has values 6 orders of magnitude higher
      
      
      ### build output data frame
      # pick minimum for each time point: 
      A.df <- data.frame (AA1 = AA1, AA2 = AA2, Ci.A1 = Ci.extract.A1, Ci.A2 = Ci.extract.A2)
      A.df$A.min <- apply(A.df[,1:2], 1, min)
      A.df$min.eq <- apply(A.df[,1:2], 1, which.min)
      
      ### solve for gsw ###
      gsw.solve <- ((M * A.df$A.min)/((input.df$Ca - input.df$comp) * (1 + ((input.df$el - input.df$ea)/Do)))) + b
      
      A.df$gsw <- gsw.solve
   }
   
   
   return(A.df)
}


##### Applying the farquhar_solver function #####
farquhar_solver(input.df = dummy, stomata = 'closed')
farquhar_solver(input.df = dummy, stomata = 'open') 

#### Check against plantecophs package... ####
dummy$el <- el(Tv = dummy$Tv)
dummy$ea <- ea(relHum = dummy$relHum, Tv = dummy$Tv)
dummy$VPD<- dummy$el-dummy$ea
dummy$Tair<-dummy$Tv-273.15

FARAO(Ca=dummy$Ca, VPD=dummy$VPD, Tair=dummy$Tair)


#### Now look at real data to compare ####
dat <- read.csv('~/Documents/git/scaling_farquhar/Aggregated_Climate_Data.csv', header=TRUE)
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

library(dplyr)
dat.check<-dat%>%filter(year==2013)%>%filter(month==8)%>%filter(day==13)%>% ### choose this date due to LAI information
   dplyr::select(hour, Tv, relHum, Ca, comp, Vcmax, el, ea, PAR, time)
dat.check<-dat.check[!duplicated(dat.check),]

dat.check<-dat.check[(dat.check$time<24),]

library(plantecophys)
dat.check$VPD<- dat.check$el-dat.check$ea
#dat.check$VPDx<-RHtoVPD(dat.check$relHum, TdegC = dat.check$Tv-273.15)
dat.check$Tair<-dat.check$Tv-273.15

FARAO(Ca=dat.check$Ca, VPD=dat.check$VPD, Tair=dat.check$Tair)

##### Applying the farquhar_solver function to real data #####
farquhar_solver(input.df = dat.check, stomata = 'closed')
farquhar_solver(input.df = dat.check, stomata = 'open') 


### And then here is the LAI info for scaling! ###
lai <- read.csv('~/Documents/git/scaling_farquhar/hf150-01-hem-lai.csv', header=TRUE)
lai<-lai%>%filter(date=="2013-08-13")
lai$lai.sum<-ave(lai$lai.masked, FUN=sum) ## 54.28

###make a data frame for FARAO out put
eco.phys.sol<-as.data.frame(FARAO(Ca=dat.check$Ca, VPD=dat.check$VPD, Tair=dat.check$Tair))
##add a column for the lai sums from harvard forest
eco.phys.sol$HF.lai<-lai$lai.sum
###multiply them together to get the footprint level A
eco.phys.sol$A.footprint<-eco.phys.sol$HF.la*eco.phys.sol$ALEAF
write.csv(eco.phys.sol, file="~/Documents/git/scaling_farquhar/plantecophys_output.csv", row.names = FALSE)
