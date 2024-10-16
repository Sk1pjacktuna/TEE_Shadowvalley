simulate_one_gen <- function(N_a, N_A, fitnessa, fitnessA, mut_rate) {
  # draw offspring according to Poisson distribution
  offsp_a <- rpois(1, N_a * fitnessa)
  offsp_A <- rpois(1, N_A * fitnessA)
  # draw new mutants according to Poisson distribution
  mut_a_to_A <- rpois(1, offsp_a * mut_rate)
  mut_A_to_a <- rpois(1, offsp_A * mut_rate)
  
  # determine new population sizes of wild type and mutant
  N_a_new <- max(offsp_a - mut_a_to_A, 0)+mut_A_to_a
  N_A_new <- max(offsp_A - mut_A_to_a, 0) + mut_a_to_A
  
  return(c(N_a_new, N_A_new))
}

simulate_pop <- function(N_init_a, N_init_A, maxfit, minfit, fitflip, mut_rate, t_max) {
  # Create the vector in which to save the results
  pop_vector <- c(N_init_a, N_init_A)
  # initiate the variables
  pop_new <- c(N_init_a, N_init_A)
  v <- 0
  # run the simulation until generation t_max
  for (i in 1:t_max+1) {
    #define fitness of a and A for each 
    # redefine the current population one generation later
    pop_new <- simulate_one_gen(pop_new[1],pop_new[2], fitnessa, fitnessA, mut_rate)
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