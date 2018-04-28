####################################################
#                 Farquhar_solver                  #
#        OEB203 project code, spring 2018          #
####################################################

##### The system of equations #####
### Equation 1: Farquhar Model - take minimum of light-limited vs. temperature limited
## CO2 limited
# A1 = (Vcmax * (Ci-comp)) / (Ci + Kc * (1 + O/K0)) - Rd   # units: umol CO2 / m2 / sec
#     (Vcmax * (Ci - comp)) / (Ci + K1*(1+K2)) - Rd

# Rd = gamma * vcmax; gamma = 0.015
# vcmax in this context is at 25 C 
# CO2-limited assimilation

A1<- function(Ci, Tv, Vmo = 35, gamma = 0.015) {
   VCmax_fixedT <- Vcmax(Vmo = Vmo, Tv = Tv)
   Rd <- gamma * VCmax_fixedT
   A <- (VCmax_fixedT * (Ci - comp(Tv))) / (Ci + K1(Tv)*(1+K2(Tv))) - Rd
   return(A)
}


## Light-limited 
# A2 = alpha * PAR * ((Ci - comp) / (Ci + 2 * comp)) - Rd
# alpha = quantum yield = 0.04

A2 <- function(Ci, Tv, PAR = 1000, Vmo = 35, alpha = 0.04, gamma = 0.015) {
   VCmax_fixedT <- Vcmax(Vmo = Vmo, Tv = Tv)
   Rd <- gamma * VCmax_fixedT
   A <- alpha * PAR * ((Ci - comp(Tv)) / (Ci + 2 * comp(Tv))) - Rd
   return(A)
}


### Equation 2: Stomatal conductance Model (Leunning); Medvigy
# gsw = ((M * A) / ((Cs - comp) * (1 + (el - es)/Do))) + b    # open stomata
#       b                                                     # closed stomata
# b <- function (Tv) {}

gsw <- function(Amin, M = 9, Do = 1.3, b = 2, Ca, Tv, open = TRUE) {
   if(open == TRUE) {
      gsw <- ((M * Amin) / ((Ca - comp(Tv)) * (1 + (el - ea)/Do))) + b
   }
   else {
      gsw <- b
   }
}


## Equation 3: 
# A = (-gsw / 1.6) * (Ci  - Ca)

## Equation 4: Water flux (not necessary?) 
# Phi = gsw * (el - ea) 


##### The calculated variables #####
#### Parameters we can calculate or get trait values for: 
## Vcmax is from Medvigy 
# Vcmax = (Vmo * exp(3000*(1/288.15 -1/Tv))) / ((1+exp(.4(Tvlo - Tv))) * (1 + exp(0.4 * (Tv - 318.15))))
Vcmax <- function (Vmo, Tvlo, Tv) {  #Vmo = 23, Tvlo = 277.15 are constants
   (Vmo * (exp(3000 * (1/288.15 -1/Tv)))) /
   ((1 + exp(.4 * (Tvlo - Tv))) * (1 + exp(0.4 * (Tv - 318.15))))
}

##### VCmax equation with Paul
Tv <- 270:320
C1 <- 1 # 2.12e-5  #missing in original coding of this equation
C2 <- 3000
Vmo <- 23
num<- Vmo * (C1 *exp(C2 * (1/288.15 -1/Tv)))
denom<-1/((1 + exp(.4 * (277.15 - Tv))) * (1 + exp(0.4 * (Tv - 318.15))))
plot(Tv, num*denom)

# Vcmax <- function(Tv) {
#    C1 <- 1
#    C2 <- 3000
#    Vmo <- 23
#    num<- Vmo * (C1 *exp(C2 * (1/288.15 -1/Tv)))
#    denom<-1/((1 + exp(.4 * (277.15 - Tv))) * (1 + exp(0.4 * (Tv - 318.15))))
#    Vc <- num * denom
#    return(Vc) 
# }

# 
# Vcmax <- function (Vmo, Tvlo, Tv) {  #Vmo, Tvlo are constants
#    (Vmo * (exp(3000 * (1/288.15 -1/Tv)))) / 
#       ((1 + exp(.4 * (Tv - Tvlo))) * (1 + exp(0.4 * (Tv - 318.15))))
# }


# Vmo = 6.3 umol/m2/sec (Medvigy 2009, table 3)
# Vmo = Vc at 25C or 15C  #Vcmax should be ~40 at 25 degrees
# Tvlo = Low temperature threshold (Kelvin) = 4.7 C = 277.85 Kelvin
# Tv = leaf temperature (Kelvin)

## new Vcmax equation (for Vcmax at 25c):
# Kattge 2007, Plant, Cell and Environment
# Vcmax = Vmo * exp((36380 * (Tv - 298))/(298*8.314*Tv))
# Vmo is between 58 and 92, for different connifers (eqn 9)
# Vcmax <- function (Vmo, Tv) {
#    Vmo * exp((36380 * (Tv - 298)) / (298 * 8.314 * Tv))
# }


# Compensation point equation from Medgivy (b.7)
# comp = 21.2 * exp(5000 * ((1/288.15) - (1/Tv)))  #CO2 compensation point where photosynthesis = respiration 
comp <- function(Tv) { 21.2 * exp(5000 * ((1/288.15) - (1/Tv))) } #units: umol/mol

# K1 = 150 * exp(6000 * ((1/288.15) - (1/Tv)))
K1 <- function(Tv) { 150 * exp(6000 * ((1/288.15) - (1/Tv))) } 


# K2 = 0.836 * exp(-1400 * ((1/288.15) - (1/Tv)))
K2 <- function(Tv) {0.836 * exp(-1400 * ((1/288.15) - (1/Tv))) }

# alpha = 0.06 # mol CO2 / mol photons   ## quantum efficiency

# Rd = gamma * Vcmax # ?? gamma = 0.015?

# el # saturated leaf concentration - 100
# ea = atmospheric water concentration - meterology?
# Saturation vapor pressure equation: el = a * exp(b*T/(T+c)) 
# T is temperature (degrees C) of the *leaf*
# a = 0.611 kPa
# b = 17.502
# c = 240.97 degrees C
# function to convert kelvin to celcius for input into the saturated leaf function (el)
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


# PAR - photon flux density, 500-2000
# J = 0 - 250 (put alpha*PAR instead of J)
# Ca = [atm CO2]

##### The Farquhar Solver Function #####
farquhar_solver <- function (input.df, stomata = c('open', 'closed')) { 
   ## input.df needs to have columns for Tv and Ca and relhum. the rest can be calculated... 
   
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
   
   
   # calculated in function
   Rd <- gamma * input.df$Vcmax
   k1.dat <- K1(input.df$Tv)
   k2.dat <- K2(input.df$Tv)
   X = k1.dat * (1 + k2.dat)
   Y <- alpha * input.df$PAR
   FF <- ((input.df$Ca - input.df$comp) * (1 + ((input.df$el - input.df$ea)/Do)))
   
   
   # relevant functions
   A1<- function(Ci, Tv, Vmo = 35, gamma = 0.015) {
      Rd <- gamma * input.df$Vcmax
      A <- (input.df$Vcmax * (Ci - comp(Tv))) / (Ci + K1(Tv)*(1+K2(Tv))) - Rd
      return(A)
   }
   
   A2 <- function(Ci, Tv, PAR, Vmo, alpha = 0.04, gamma = 0.015) {
      Rd <- gamma * input.df$Vcmax
      A <- (alpha * PAR * ((Ci - comp(Tv)) / (Ci + 2 * comp(Tv)))) - Rd
      return(A)
   }
   
   ### stomata closed ###
   if(stomata == 'closed') {
      ### CO2 limited case, coefficients for polyroot function
      aa <- (b * X * input.df$Ca/1.6) + (input.df$Vcmax * input.df$comp) + (Rd * X)
      bb <- (Rd) + (b * input.df$Ca / 1.6) - (b * X / 1.6) - (input.df$Vcmax)
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
      AA1 <- A1(Ci.extract.A1, input.df$Tv, Vmo = Vmo, gamma = gamma)
      
      ### Light limited case, coefficients for polyroot function
      aa <- ((b * 2 * input.df$comp * input.df$Ca / 1.6) + (Rd * 2 * input.df$comp) + (Y * input.df$comp))
      bb <- (Rd + (b * input.df$Ca / 1.6) - (b * 2 * input.df$comp / 1.6) - Y) 
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
      AA2 <- A2(Ci.extract.A2, Tv = input.df$Tv, PAR = input.df$PAR, Vmo = Vmo, alpha = alpha)  # only works if PAR has values 6 orders of magnitude higher
      
      
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
      aa <- ((b * X * input.df$Ca) + (1.6 * FF * input.df$Vcmax * input.df$comp) + (1.6 * FF * Rd * X) + 
            (M * input.df$Vcmax * input.df$comp * input.df$Ca) + (M * Rd * X * input.df$Ca))
      bb <- ((b * input.df$Ca) - (b * X) - (1.6 * FF * input.df$Vcmax) + (1.6 * FF * Rd) + (M * input.df$Vcmax * input.df$Ca) +
            (M * input.df$Vcmax * input.df$comp) + (Rd * M * input.df$Ca) + (M * Rd * X))
      cc <- ((-b * FF) + (M * input.df$Vcmax) - (M * Rd))
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
      AA1 <- A1(Ci.extract.A1, input.df$Tv, Vmo = Vmo, gamma = gamma)
      
      
      ### Light-limited coefficients
      aa <- ((b * FF * input.df$Ca * 2 * input.df$comp) + (1.6 * FF * Y * input.df$comp) + (1.6 * FF * Rd * 2 * input.df$comp) - 
            (M * input.df$Ca * Y * input.df$comp) - (M * Rd * 2 * input.df$comp * input.df$Ca))
      bb <- ((b * FF * input.df$Ca) - (b * FF * 2 * input.df$comp) - (1.6 * FF * Y) + (1.6 * FF * Rd) + (M * Y * input.df$Ca) +
            (M * Y * input.df$comp) - (M * Rd * input.df$Ca) + (M * Rd * 2 * input.df$comp))
      cc <- ((-b * FF) - (M * Y) + (M * Rd))
      
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
      AA2 <- A2(Ci.extract.A2, Tv = input.df$Tv, PAR = input.df$PAR, Vmo = Vmo, alpha = alpha)  # only works if PAR has values 6 orders of magnitude higher

      
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

### GSW troubleshooting 
M<-9
input.df <- dummy
A.min <- closed.sol$A.min
Do <- 1.3
((M * A.min)/((input.df$Ca - input.df$comp) * (1 + ((input.df$el - input.df$ea)/Do))))

##### The Data #####
### Dummy data ###
dummy <- data.frame(time = c(0:23))
dummy$Tv <- c(seq(280,300, length.out = 12), seq(300, 280, length.out = 12))
dummy$relHum <- c(rep(70, 3), 80, 80, 90, 100, 100, 90, 80, 70, 60, 50, 40, 50, 60, rep(70, 8))
dummy$Ca <- rep(400, 24)
dummy$comp <- comp(Tv = dummy$Tv)
dummy$Vcmax <- Vcmax(Vmo = 6.3, Tv = dummy$Tv, Tvlo = 277.85)
dummy$el <- el(Tv = dummy$Tv)
dummy$ea <- ea(relHum = dummy$relHum, Tv = dummy$Tv)
dummy$PAR <- c(rep(0, 5), seq(0, 0.002, length.out = 7), seq(0.002, 0, length.out = 7), rep(0,5))
# more realistic dummy PAR: 
dummy$PAR <- c(rep(0, 5), seq(0, 2000, length.out = 7), seq(2000, 0, length.out = 7), rep(0, 5))



### Real data ###
dat <- read.csv('Aggregated_Climate_Data.csv', header=TRUE)
dat$Tv <- dat$Air_Temp_K 
dat$relhum <- dat$Relative_Humidity_Percent
dat$Ca <- dat$Atmospheric_CO2
dat$comp <- comp(Tv = dat$Tv)
dat$Vcmax <- Vcmax(Vmo = 35, Tv = dat$Tv)
dat$el <- el(Tv = dat$Tv)
dat$ea <- ea(relHum = dat$relhum, Tv = dat$Tv)
dat$PAR <- dat$Par_moles_m2_s * 1000000 ### only keep this until push updated to git 

##### Applying the farquhar_solver function #####
farquhar_solver(input.df = dummy, stomata = 'closed')
farquhar_solver(input.df = dummy, stomata = 'open') 
