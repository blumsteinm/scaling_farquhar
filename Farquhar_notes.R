##### Farquhar Modeling equations

### Equation 1: Farquhar Model - take minimum of light-limited vs. temperature limited
# A1 = (Vcmax * (Ci-comp)) / (Ci + Kc * (1 + O/K0)) - Rday   # units: umol CO2 / m2 / sec
#     (Vcmax * (Ci - comp)) / (Ci + K1*(1+K2)) - Rday



A1<- function(Ci, Tv, Vmo = 35, Tvlo = 300, gamma = 0.015) {
   VCmax_fixedT <- Vcmax(Vmo = Vmo, Tvlo = Tvlo, Tv = Tv)
   Rday <- gamma * VCmax_fixedT
   A <- (VCmax_fixedT * (Ci - comp(Tv))) / (Ci + K1(Tv)*(1+K2(Tv))) - Rday
   return(A)
}


# A2 = alpha * PAR * ((Ci - comp) / (Ci + 2 * comp)) - Rday
# alpha = quantum yield = 0.04

A2 <- function(Ci, Tv, PAR = 1000, Vmo = 35, Tvlo = 300, alpha = 0.04, gamma = 0.015) {
   VCmax_fixedT <- Vcmax(Vmo = Vmo, Tvlo = Tvlo, Tv = Tv)
   Rday <- gamma * VCmax_fixedT
   A <- alpha * PAR * ((Ci - comp(Tv)) / (Ci + 2 * comp(Tv))) - Rday
   return(A)
}

Amin <- function (Ci, Tv, PAR = 100, Vmo = 35, Tvlo = 300, gamma = 0.015, alpha = 0.04) {
   A1val <- A1(Ci = Ci, Vmo = Vmo, Tvlo = Tvlo, Tv = Tv, gamma = gamma)
   A2val <- A2(Ci = Ci, PAR = PAR, Vmo = Vmo, Tvlo = Tvlo, Tv = Tv, alpha = alpha, gamma = gamma)
   minval <- min(A1val, A2val)
   return(minval)
}

Amin(Ci, Tv = 300)


### Equation 2: Stomatal conductance Model
# gsw = ((M * A) / ((Cs - comp) * (1 + (el - es)/Do))) + b    # open stomata
#       b                                                     # closed stomata
# b <- function (Tv) {}

gsw <- function(M, Do, b, Ca, Tv, open = TRUE) {
   if(open == TRUE) {
      gsw <- ((M * Amin) / ((Ca - comp(Tv)) * (1 + (el - es)/Do))) + b
   }
   else {
      gsw <- b
   }
}

# M = stomatal slope
# b = cuticular conductance
# Ca = atmospheric composition of CO2


## Equation 3: 
# A = (-gsw / 1.6) * (Ci  - Ca)

## Equation 4: Water flux
# Phi = gsw * (el - ea) 

# Solve for: A, Ci, gsw, Phi





#### Parameters we can calculate or get trait values for: 
# Vcmax = (Vmo * exp(3000*(1/288.15 -1/Tv))) / ((1+exp(.4(Tvlo - Tv))) * (1 + exp(0.4 * (Tv - 318.15))))
Vcmax <- function (Vmo, Tvlo, Tv) {  #Vmo, Tvlo are constants 
   (Vmo * exp(3000 * (1/288.15 -1/Tv))) / ((1 + exp(.4 * (Tvlo - Tv))) * (1 + exp(0.4 * (Tv - 318.15))))
}
# Vmo = Vc at 25C or 15C  #Vcmax should be ~40 at 25 degrees
# Tvlo = Low temperature threshold (Kelvin)
# Tv = leaf temperature (Kelvin)

# comp = 21.2 * exp(5000 * ((1/288.15) - (1/Tv)))  #CO2 compensation point where photosynthesis = respiration 
comp <- function(Tv) { 21.2 * exp(5000 * ((1/288.15) - (1/Tv))) } #units: umol/mol

# K1 = 150 * exp(6000 * ((1/288.15) - (1/Tv)))
K1 <- function(Tv) { 150 * exp(6000 * ((1/288.15) - (1/Tv))) } 


# K2 = 0.836 * exp(-1400 * ((1/288.15) - (1/Tv)))
K2 <- function(Tv) {0.836 * exp(-1400 * ((1/288.15) - (1/Tv))) }

# alpha = 0.06 # mol CO2 / mol photons   ## quantum efficiency

# Rday = gamma * Vcmax # ?? gamma = 0.015?

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
   relHum / el(Tv)
}


# PAR - photon flux density, 500-2000
# J = 0 - 250 (put alpha*PAR instead of J)
# Ca = [atm CO2]

# Do = optimized 

##### Testing Farquhar (Equation 1 above) 
Tv  = seq(274, 315, 1)
Ci = seq(100, 400, 1)


plot(Ci, A1(Ci, Tv = 298))


plot(Ci, A2(Ci, Tv = 298))


points(Ci, A1(Ci))
plot(Ci, A2(Ci), col='red', ylim=c(0, max(A2(Ci))))

plot(Tv, Vcmax(Vmo = 30, Tvlo = 300, Tv = Tv), col='blue')



##### Dummy data #####
dummy <- data.frame(time = c(0:23))
dummy$Tv <- c(seq(280,300, length.out = 12), seq(300, 280, length.out = 12))
dummy$relHum <- c(70, 70, 70, 80, 90, 100, 100, 100, 90, 80, 70, 60, 50, 40, 50, 60, rep(70, 8))
dummy$comp <- comp(Tv = dummy$Tv)
dummy$Vcmax <- Vcmax(Vmo = 35, Tvlo = 300, Tv = dummy$Tv)
dummy$el <- el(Tv = dummy$Tv)
dummy$ea <- ea(relHum = dummy$relHum, tV = dummy$Tv)
dummy$Ca <- rep(400, 24)

