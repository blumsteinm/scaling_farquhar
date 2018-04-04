# How well does farquhur model scale? (key uncertainties, errors, etc.)


farquhar <- function(V.max=50, J.max=100, APAR=500, c.i=30) {
  
  # Model inputs:
  # V.max = maximum rubisco-limited rate in micromoles per (m^2 sec)
  # J.max = maximum light-limited rate in micromoles per (m^2 sec)
  # APAR = absorbed photosynthetically-active radiation in micromoles per (m^2 sec)
  # c.i = intercellular CO2 partial pressure in Pascals (roughly ppm/10)
  
  # Some local parameters we need
  p.sfc <- 101325  # surface air pressure (Pascals)
  gamma <- 3. # CO2 compensation point (Pascals)
  O.i <- 0.209 * p.sfc  # oxygen partial pressure in chloroplast
  K.c <- 30 # Michaelis-Menten constant for carboxylation (Pascals)
  K.o <- 30000 # Michaelis-Menten constant for oxidation (Pascals)
  
  # Solution of quadratic (Bonan 17.8)
  a <- 0.7
  b <- -(J.max + 0.385 * APAR)
  c <- 0.385 * J.max * APAR
  J.1 <- (-b + sqrt(b^2 - 4 * a * c) ) / (2 * a)
  J.2 <- (-b - sqrt(b^2 - 4 * a * c) ) / (2 * a)
  J <- min(J.1, J.2)
  
  # Rubisco-limited rate of photosynthesis
  w.c <- V.max * (c.i - gamma) / (c.i + K.c * (1 + O.i/K.o))  # Bonan 17.6
  
  # Light-limited rate of photosynthesis
  w.j <- J * (c.i - gamma) / (4 * (c.i + 2 * gamma))            # Bonan 17.7
  
  # Sink-limited rate of photosynthesis
  w.s <- V.max / 2
  
  # Dark respiration 
  R.d <- 0.015 * V.max
  
  # Net assimilation 
  A.n <- min(w.c, w.j, w.s)-R.d
  
  return(A.n)
}


plot.light.response <- function(V.max=50, J.max=100, c.i=30) {
  A.n <- replicate(100, NA)
  light <- replicate(100, NA)
  for (i in 1:100) {
    light[i] <- 10 * i
    A.n[i] <- farquhar(V.max=V.max, J.max=J.max, c.i=30,  APAR=light[i])
  }
  orig.par <- par(no.readonly=TRUE)
  par(mar=c(5,5,4,2)+0.1)
  plot(light, A.n, main='Light Response of Photosynthesis', 
       xlab=expression(paste('PPFD (',mu, 'mol photons ',m^{-2},s^{-1},')')),
       ylab=expression(paste(mu, 'mol ', CO[2], ' ',m^{-2}, s^{-1},)),
       cex.main=2, cex.axis=1.5, cex.lab=1.5)
  par(orig.par)
}

plot.CO2.response <- function(V.max=50, J.max=100, APAR=500) {
  A.n <- replicate(100, NA)
  CO2 <- replicate(100, NA)
  for (i in 1:100) {
    CO2[i] <-  10 * i
    A.n[i] <- farquhar(V.max=V.max, J.max=J.max, c.i=CO2[i]/10,  APAR=500)
  }
  orig.par <- par(no.readonly=TRUE)
  par(mar=c(5,5,4,2)+0.1)
  plot(CO2, A.n, main='CO2 Response of Photosynthesis', xlab='Intercellular CO2 (ppm)',
       ylab=expression(paste(mu, 'mol ', CO[2], ' ',m^{-2}, s^{-1},)),
       cex.main=2, cex.axis=1.5, cex.lab=1.5)
  par(orig.par)
}






