simulate_one_gen_migrants <- function(N_a, N_A, decay_rate, sel_coeff, mut_rate, lambda_a,lambda_A) {
  # draw offspring according to Poisson distribution
  offsp_a <- rpois(1, N_a * (1-decay_rate))
  offsp_A <- rpois(1, N_A * (1-decay_rate+sel_coeff))
  # draw new mutants according to Poisson distribution
  mut_a_to_A <- rpois(1, offsp_a * mut_rate)
  mut_A_to_a <- rpois(1, offsp_A * mut_rate)
  # draw migrants for each genotype
  migrants_a <- rpois(1,lambda_a)
  migrants_A <- rpois(1,lambda_A)
  # determine new population sizes of wild type and mutant
  N_a_new <- max(offsp_a - mut_a_to_A, 0) + mut_A_to_a + migrants_a
  N_A_new <- max(offsp_A - mut_A_to_a, 0) + mut_a_to_A + migrants_A
  
  return(c(N_a_new, N_A_new))
}

simulate_pop_switching_migrants <- function(N0_a, N0_A, decay_rate, sel_coeff, flipswitch, migration_min, migration_max, mut_rate, t_max) {
  # Create the vector in which to save the results
  pop_vector <- c(N0_a, N0_A)
  # initiate the variables
  pop_new <- c(N0_a, N0_A)
  
  migration <- sample(c(migration_min,migration_max))
  migration_a <- migration[1]
  migration_A <- migration[2]
  # run the simulation until generation t_max
  v <- 0
  for (i in 1:t_max+1) {
    #define fitness of a and A for each iteration
    if (v<flipswitch){
      v <- v+1
    } else {
      migration <- migration_a
      migration_a <- migration_A
      migration_A <- migration
      v <- 0
    }
    
    # redefine the current population one generation later
    pop_new <- simulate_one_gen_migrants(pop_new[1],pop_new[2], decay_rate, sel_coeff, mut_rate,migration_a,migration_A)
    # add the new population sizes to the output vector
    pop_vector <- rbind(pop_vector,pop_new)
    # condition to stop the simulation before t_max: either the population exceeds 1.5 times the original population size, or it goes extinct
    
    if (pop_new[1]+pop_new[2]>=1.5*(N0_a+ N0_A) | pop_new[1]+pop_new[2]==0) break
  }
  
  # define the row and column names of the output vector
  rownames(pop_vector) <- (0:t_max)[1:length(pop_vector[,1])] # note that the vector has to be cut if the simulation stopped early
  colnames(pop_vector) <- c("a","A")
  # return the result
  return(pop_vector)
}

# set some parameters to fixed values
init_a <- 1000
init_A <- 0
m_rate <- 0.001
max_gen <- 1000
decay_rate <- 0.1
sel_coef <- 0.2
migration_min <- 0.1
migration_max <- 5
mut_rate <- 0

# determine how often to run the simulation for each set of parameters
no_replicates <- 100

# set parameters to vary
switch_gens <- c(1,3,5,10,15)

# initialize data table - where to collect the results
data_table <- c()
# run the simulation across all chosen parameters
# loop over switch generations
for(flipswitch in switch_gens){
    # reset counter
    i<-1
    repeat {
      # increase counter by one
      i<-i+1
      # run the simulation once
      one_run <- simulate_pop_switching_migrants(init_a,init_A,decay_rate,sel_coef,flipswitch, migration_min, migration_max, mut_rate, max_gen)
      # determine total population sizes
      total_size <- one_run[,1]+one_run[,2]
      #number of generation at which 1.5*n0 is achieved
      number_of_generations <- length(one_run[,1])
      # determine minimum population size
      min_size <- min(total_size)
      # determine (first) generation at which this population size occurred
      min_gen <- as.numeric(which(total_size==min_size)[1])
      # enter the data into the table
      data_table <- rbind(data_table,c(flipswitch,min_gen,min_size,number_of_generations)) # note that we add the varying parameters (decay rate and selection coefficient) to the table too
      # stop the repeated computation after no_replicates times
      if(i>no_replicates) break
    }
}
colnames(data_table) <- c("flipswitch","min_gen","min_size","number_of_generations")

m_data_table <- c("flipswitch_av" = mean(data_table[,1]),"min_gen_av" = mean(data_table[,2]),"min_size_av" = mean(data_table[,3]),"av_no_gen" = mean(data_table[,4]))
print(m_data_table)

