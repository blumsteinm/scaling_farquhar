##### Farquhar Modeling equations

## Equation 1: Farquhar Model - take minimum
# A1 = (Vcmax * (Ci-comp)) / (Ci + Kc * (1 + O/K0)) - Rday   # units: umol CO2 / m2 / sec
#     (Vcmax * (Ci - comp)) / (Ci + K1*(1+K2)) - Rday

# A2 = alpha * PAR * ((Ci - comp) / (Ci + 2 * comp)) - Rday

## Equation 2: Stomatal conductance Model
# gsw = ((M * A) / ((Cx - comp) * (1 + (el - es)/Do))) + b    # open stomata
#       b                                                     # closed stomata


## Equation 3: 
# A = -gsw / 1.6 * (Ci  - Ca)

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

plot(Tv, Vcmax(Vmo = 30, Tvlo = 300, Tv = Tv), col='blue')

# comp = 21.2 * exp(5000 * ((1/288.15) - (1/Tv)))  #CO2 compensation point where photosynthesis = respiration 
comp <- function(Tv) { 21.2 * exp(5000 * ((1/288.15) - (1/Tv))) } #units: umol/mol

# K1 = 150 * exp(6000 * ((1/288.15) - (1/Tv)))
K1 <- function(Tv) { 150 * exp(6000 * ((1/288.15) - (1/Tv))) } 


# K2 = 0.836 * exp(-1400 * ((1/288.15) - (1/Tv)))
K2 <- function(Tv) {0.836 * exp(-1400 * ((1/288.15) - (1/Tv))) }

# alpha = 0.06 # mol CO2 / mol photons   ## quantum efficiency

# Rday = gamma * Vcmax # ?? gamma = 0.015?

# el # saturated leaf concentration - 100
# ea # ea = atmospheric water concentration - meterology?
# PAR - photon flux density, 500-2000
# J = 0 - 250 (put alpha*PAR instead of J)
# Ca = [atm CO2]

# Do = optimized 
Tv  = seq(274, 315, 1)


K1 <- function(Tv) { 150 * exp(6000 * ((1/288.15) - (1/Tv)))} 
K1 (Tv)


##### Testing Farquhar (Equation 1 above) 
Tv  = seq(274, 315, 1)
Ci = seq(100, 400, 1)

A1<- function(Ci = Ci, Vmo = 35, Tvlo = 300, Tv = Tv, gamma = 0.015) {
   VCmax_fixedT <- Vcmax(Vmo = Vmo, Tvlo = Tvlo, Tv = Tv)
   Rday <- gamma * VCmax_fixedT
   
   A <- (VCmax_fixedT * (Ci - comp(Tv))) / (Ci + K1(Tv)*(1+K2(Tv))) - Rday
   return(A)
}

plot(Ci, A1(Ci))


A2 <- function(Ci = Ci, PAR = 1000, Vmo = 35, Tvlo = 300, Tv = Tv, alpha = 0.06, gamma = 0.015) {
   VCmax_fixedT <- Vcmax(Vmo = Vmo, Tvlo = Tvlo, Tv = Tv)
   Rday <- gamma * VCmax_fixedT
   
   A <- alpha * PAR * ((Ci - comp(Tv)) / (Ci + 2 * comp(Tv))) - Rday
   return(A)
}

plot(Ci, A2(Ci))


points(Ci, A1(Ci))
plot(Ci, A2(Ci), col='red', ylim=c(0, max(A2(Ci))))