####### Deep-time scaling #####

# call in modeled CO2 and Tv data
petm12 <- read.csv('petm12.csv', header = TRUE)
petm22 <- read.csv('petm22.csv', header = TRUE)
petm60 <- read.csv('petm60.csv', header = TRUE)

# add PAR values: 
petm12$PAR <- 600
petm22$PAR <- 600
petm60$PAR <- 600


open.12 <- farquhar_solver(input.df = petm12, stomata = 'open')
closed.12 <- farquhar_solver(input.df = petm12, stomata = 'closed')
open.22 <- farquhar_solver(input.df = petm22, stomata = 'open')
closed.22 <- farquhar_solver(input.df = petm22, stomata = 'closed')
open.60 <- farquhar_solver(input.df = petm60, stomata = 'open')
closed.60 <- farquhar_solver(input.df = petm60, stomata = 'closed')



plot(petm12$time_yr, closed.12$A.min, ylim=c(0, max(closed.12$A.min)))  
points(petm12$time_yr, open.12$A.min, col='blue')   #woah, we actually become light-limited in the middle of this!
points(petm22$time_yr, closed.22$A.min, col = 'blue')
points(petm60$time_yr, closed.60$A.min, col='red')


petm12$PAR <- 600
petm22$PAR <- 600
petm60$PAR <- 600

closed.12 <- farquhar_solver(input.df = petm12, stomata = 'closed')
closed.22 <- farquhar_solver(input.df = petm22, stomata = 'closed')
closed.60 <- farquhar_solver(input.df = petm60, stomata = 'closed')



### Make plots ###

# environment
jpeg('PETM_environment.jpg', height = 4, width = 10, units = 'in', res = 200)
par(mfrow = c(1,2))

plot(petm12$time_yr, petm12$Ca, type = 'l', col = 'black', xlab = 'Time post-carbon injection (years)', ylab = 'atmospheric CO2 (ppm)', ylim = c(500,5000), main = 'Atmospheric CO2 across the PETM')
points(petm22$time_yr, petm22$Ca, type = 'l', col = 'blue')
points(petm60$time_yr, petm60$Ca, type = 'l', col = 'red')

legend('topright', legend = c('Volcano', 'Fossil fuels', 'Methane'), lty = 1, col = c('black', 'blue', 'red'), title = 'Carbon source', cex = 0.8)

plot(petm12$time_yr, petm12$temp_C, type = 'l', col = 'black', xlab = 'Time post-carbon injection (years)', ylab = 'Temperature (C)', ylim = c(20, max(petm12$temp_C)), main = 'Temperature change across the PETM')
points(petm22$time_yr, petm22$temp_C, type = 'l', col = 'blue')
points(petm60$time_yr, petm60$temp_C, type = 'l', col = 'red')

legend('topright', legend = c('Volcano', 'Fossil fuels', 'Methane'), lty = 1, col = c('black', 'blue', 'red'), title = 'Carbon source', cex = 0.8)
dev.off()

#assimilation 
jpeg('PETM_assimilation.jpg', height = 4, width = 10, units = 'in', res = 200)

par(mfrow = c(1,2))

petm12$PAR <- 600
petm22$PAR <- 600
petm60$PAR <- 600

closed.12 <- farquhar_solver(input.df = petm12, stomata = 'closed')
closed.22 <- farquhar_solver(input.df = petm22, stomata = 'closed')
closed.60 <- farquhar_solver(input.df = petm60, stomata = 'closed')

plot(petm12$time_yr, closed.12$A.min, type = 'l', ylim=c(16,32), xlab = 'Time post-carbon injection (years)', ylab = 'Carbon Assimilation', main = 'Leaf-based Carbon assimilation, PAR=600')  #woah, we actually become light-limited in the middle of this!
#points(petm12$time_yr, open.12$A.min, col='blue')   
points(petm22$time_yr, closed.22$A.min, col = 'blue', type = 'l')
points(petm60$time_yr, closed.60$A.min, col='red', type = 'l')
legend('topright', legend = c('Volcano', 'Fossil fuels', 'Methane'), lty = 1, col = c('black', 'blue', 'red'), title = 'Carbon source', cex = 0.8)

petm12$PAR <- 800
petm22$PAR <- 800
petm60$PAR <- 800

closed.12 <- farquhar_solver(input.df = petm12, stomata = 'closed')
closed.22 <- farquhar_solver(input.df = petm22, stomata = 'closed')
closed.60 <- farquhar_solver(input.df = petm60, stomata = 'closed')

plot(petm12$time_yr, closed.12$A.min, type = 'l', ylim=c(16,32), xlab = 'Time post-carbon injection (years)', ylab = 'Carbon Assimilation', main = 'Leaf-based Carbon assimilation, PAR=800')  #woah, we actually become light-limited in the middle of this!
#points(petm12$time_yr, open.12$A.min, col='blue')   
points(petm22$time_yr, closed.22$A.min, col = 'blue', type = 'l')
points(petm60$time_yr, closed.60$A.min, col='red', type = 'l')
legend('topright', legend = c('Volcano', 'Fossil fuels', 'Methane'), lty = 1, col = c('black', 'blue', 'red'), title = 'Carbon source', cex = 0.8)


dev.off()



