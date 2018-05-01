####### Deep-time scaling #####

# call in modeled CO2 and Tv data
petm12 <- read.csv('petm12.csv', header = TRUE)
petm22 <- read.csv('petm22.csv', header = TRUE)
petm60 <- read.csv('petm60.csv', header = TRUE)

# add PAR values: 
petm12$PAR <- 800
petm22$PAR <- 800
petm60$PAR <- 800


# CO2 plot
plot(petm12$time_yr, petm12$Ca)
points(petm22$time_yr, petm22$Ca, col= 'red')
points(petm60$time_yr, petm60$Ca, col='blue')


open.12 <- farquhar_solver(input.df = petm12, stomata = 'open')
closed.12 <- farquhar_solver(input.df = petm12, stomata = 'closed')

plot(petm12$time_yr, closed.12$A.min, ylim=c(0, max(closed.12$A.min)))  
points(petm12$time_yr, open.12$A.min, col='blue')   #woah, we actually become light-limited in the middle of this!


