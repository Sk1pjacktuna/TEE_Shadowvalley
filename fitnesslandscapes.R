simulate_one_gen <- function(N_a, N_A, decay_rate, sel_coeff, mut_rate) {
  # draw offspring according to Poisson distribution
  offsp_a <- rpois(1, N_a * (1-decay_rate))
  offsp_A <- rpois(1, N_A * (1-decay_rate+sel_coeff))
  # draw new mutants according to Poisson distribution
  mut_a_to_A <- rpois(1, offsp_a * mut_rate)
  
  # determine new population sizes of wild type and mutant
  N_a_new <- max(offsp_a - mut_a_to_A, 0)
  N_A_new <-  offsp_A + mut_a_to_A
  
  return(c(N_a_new, N_A_new))
}
simulate_pop <- function(N_init_a, N_init_A, decay_rate, sel_coeff, mut_rate, t_max) {
  # Create the vector in which to save the results
  pop_vector <- c(N_init_a, N_init_A)
  # initiate the variables
  pop_new <- c(N_init_a, N_init_A)
  
  # run the simulation until generation t_max
  for (i in 1:t_max+1) {
    # redefine the current population one generation later
    pop_new <- simulate_one_gen(pop_new[1],pop_new[2], decay_rate, sel_coeff, mut_rate)
    # add the new population sizes to the output vector
    pop_vector <- rbind(pop_vector,pop_new)
    # condition to stop the simulation before t_max: either the population exceeds 1.5 times the original population size, or it goes extinct
    if (pop_new[1]+pop_new[2]>=1.5*(N_init_a+ N_init_A) | pop_new[1]+pop_new[2]==0) break
  }
  
  # define the row and column names of the output vector
  rownames(pop_vector) <- (0:t_max)[1:length(pop_vector[,1])] # note that the vector has to be cut if the simulation stopped early
  colnames(pop_vector) <- c("a","A")
  # return the result
  return(pop_vector)	
}

simulation_length <- 20 #needs to be divisible by 2*fit_flip
fitflip <- 5
season_changes <- simulation_length/(2*fitflip)
fit_shadowA <- 1.3
fit_shadowB <- 0.9
fit_sunA <- 0.9
fit_sunB <- 1.3
x <- seq(0,30,by = 0.1)
fitnessA <- rep(c(rep(fit_shadowA,fitflip),rep(fit_sunA,fitflip)), season_changes)
fitnessB <- rep(c(rep(fit_shadowB,fitflip),rep(fit_sunB,fitflip)), season_changes)
plot(fitnessA, col = "blue", xlab = "generation",ylab = "fitness", type = "l", main = "four (2) season valley")
lines(fitnessB, col = "red")
legend(0,1.1, legend=rep(c("mut","wt"),5), col = c("blue","red"), pch = 1:10)
fitnessA <- fit_shadowA
fitnessB <- fit_shadowB

plot(c(fitnessA,fitnessB), ylab = "fitness", axes = F, xlab = "Genotype", main="shadow valley" )
axis(1,at = c(1,2), labels = c("Mutated genotype", "wt"))
axis(2)


