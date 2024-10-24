
#####simulate averages of: generation at which minimal population size is reached, minimal population size, number of simulated generations, maximum population size ##########
data <- read.csv("analysis2.0.csv")

# get the average generation at which the minimum population size is reached for each number of migrants
i = 1
ava_gen_minpop <- c()
for ( i in 0:10){
  ava_gen_minpop <- c(ava_gen_minpop,mean(data$generation_with_min_pop_size[data$migration == i]))
}

i = 1
ava_minpop <- c()
for ( i in 0:10){
  ava_minpop <- c(ava_minpop,mean(data$min_pop_size[data$migration == i]))
}

i = 1
ava_sim_gen <- c()
for ( i in 0:10){
  ava_sim_gen <- c(ava_sim_gen,mean(data$simulated_gen[data$migration == i]))
}

i = 1
ava_maxpop <- c()
for ( i in 0:10){
  ava_maxpop <- c(ava_maxpop,mean(data$maximal_pop_size[data$migration == i]))
}


par(mfrow = c(2,2))
plot(0:10, ava_gen_minpop, type = "l", xlab = "migrants", ylab = "gens to minpop", main = "average gens to min population size")
plot(0:10, ava_sim_gen, type = "l", xlab = "migrants", ylab = "simulated gens", main = " average number of simulated generations")
plot(0:10, ava_maxpop, type = "l", xlab = "migrants", ylab = "maxpop", main = "average max population size")
plot(0:10, ava_minpop, type = "l", xlab = "migrants", ylab = "minpop", main = "average min population size")

#####



