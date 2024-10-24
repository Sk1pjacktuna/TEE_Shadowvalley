
data <- read.csv("analysis.csv")

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

p1 <- plot(0:10, ava_gen_minpop, type = "l", xlab = "No. migrants", ylab = "av. gens to minpop")
p2<- plot(0:10, ava_maxpop, type = "l", xlab = "No. migrants", ylab = "av. maxpop")
p3 <- plot(0:10, ava_minpop, type = "l", xlab = "No. migrants", ylab = "av.minpop")
p4 <- plot(0:10, ava_sim_gen, type = "l", xlab = "No. migrants", ylab = "av. simulated gens")
library(patchwork)


gridExtra::grid.arrange(p1,p2,ncol=2)



