##########################averages##################

simulate_one_gen_hardy_weinberg <- function(N_aa,N_Aa, N_AA, fitnessaa, fitnessAa, fitnessAA, avgmigrants, mut_rate) {
  # calculate total pop size
  tot_pop <- N_aa+N_Aa+N_AA
  # calculate average fitness
  avg_fit <- weighted.mean(c(fitnessaa, fitnessAa, fitnessAA), c(N_aa,N_Aa, N_AA))
  # calculate deterministic next gen pop size
  next_gen_tot_pop <- tot_pop*avg_fit
  # calculate the allele frequnecies
  a_tot <- 2*N_aa+N_Aa
  A_tot <- 2*N_AA+N_Aa
  a_to_A_mut <- max(rpois(1,a_tot*mut_rate),a_tot)
  A_to_a_mut <- max(rpois(1,A_tot*mut_rate),A_tot)
  p <- (a_tot+A_to_a_mut-a_to_A_mut)/(2*tot_pop)
  q <- 1-p
  # Next generation genotype frequencies; selection is applied here
  det_aa_next <- (p^2)*fitnessaa/avg_fit
  det_Aa_next <- (2*p*q)*fitnessAa/avg_fit
  det_AA_next <- (q^2)*fitnessAA/avg_fit
  
  #migrationrate
  aa_migrants <- rpois(1, avgmigrants)
  # draw offspring according to Poisson distribution
  offsp_aa <- rpois(1, next_gen_tot_pop*det_aa_next) + aa_migrants
  offsp_Aa <- rpois(1, next_gen_tot_pop*det_Aa_next)
  offsp_AA <- rpois(1, next_gen_tot_pop*det_AA_next)
  
  
  return(c(offsp_aa, offsp_Aa,offsp_AA))
}
simulate_pop_HW <- function(N_init_aa, N_init_Aa, N_init_AA, fitnessaa, fitnessAa, fitnessAA, avgmigrants, mut_rate, t_max) {
  # Create the vector in which to save the results
  pop_vector <- c(N_init_aa, N_init_Aa, N_init_AA)
  # initiate the variables
  pop_new <- c(N_init_aa, N_init_Aa, N_init_AA)
  v <- 0
  # run the simulation until generation t_max
  for (i in 1:t_max+1) {
    #define fitness of a and A for each iteration
    # redefine the current population one generation later
    pop_new <- simulate_one_gen_hardy_weinberg(pop_new[1],pop_new[2], pop_new[3], fitnessaa, fitnessAa, fitnessAA,avgmigrants, mut_rate)
    # add the new population sizes to the output vector
    pop_vector <- rbind(pop_vector,pop_new)
    # condition to stop the simulation before t_max: either the population exceeds 1.5 times the original population size, or it goes extinct
    if (pop_new[1]+pop_new[2]+pop_new[3]>=5*(N_init_aa+ N_init_Aa+ N_init_AA) | pop_new[1]+pop_new[2]+pop_new[3]==0) break
  }
  
  # define the row and column names of the output vector
  rownames(pop_vector) <- (0:t_max)[1:length(pop_vector[,1])] # note that the vector has to be cut if the simulation stopped early
  colnames(pop_vector) <- c("aa","Aa","AA")
  # return the result
  return(pop_vector)	
}

# Test the function

# set some parameters to fixed values
#init_a <- 1000
#init_A <- 0
#m_rate <- 0.001
#max_gen <- 1000
#decay_rate <- 0.1
#sel_coef <- 0.2
#igration_min <- 0.1
#migration_max <- 5
#mut_rate <- 0
N_aa <- 1000
N_Aa <- 0
N_AA <- 100
fitnessaa <- 0.3
fitnessAa <- 0.3
fitnessAA <- 1
avgmigrants <- 4
mut_rate <- 0.01
t_max <- 100

# determine how often to run the simulation for each set of parameters
no_replicates <- 100

# initialize data table - where to collect the results
data_table <- c()
# run the simulation across all chosen parameters
# loop over switch generations
migrants <- 0:10

for(avgmigrants in migrants){
    # reset counter
    i<-1
    repeat {
      # increase counter by one
      i<-i+1
      # run the simulation once
      one_run <- simulate_pop_HW(N_aa,N_Aa, N_AA, fitnessaa, fitnessAa, fitnessAA, avgmigrants, mut_rate, t_max)
      # determine total population sizes
      total_size <- one_run[,1] + one_run[,2] + one_run[,3]
      #number of generation at which 1.5*n0 is achieved
      number_of_generations <- length(one_run[,1])
      # determine minimum population size
      min_size <- min(total_size)
      # determine (first) generation at which this population size occurred
      min_gen <- as.numeric(which(total_size==min_size)[1])
      # enter the data into the table
      data_table <- rbind(data_table,c(avgmigrants,min_gen,min_size,number_of_generations)) # note that we add the varying parameters (decay rate and selection coefficient) to the table too
      # stop the repeated computation after no_replicates times
      if(i>no_replicates) break
    }
}


# colnames(data_table) <- c("flipswitch","min_gen","min_size","number_of_generations")

print(one_run)
print(one_run[1])
print(data_table)

View(data_table)
# gives average of the values recieved in one run (for now, it's going to do it for the whole population soon)
data_table_av <- c("No_migrants" = mean(data_table[,1]), "min_gen_av" = mean(data_table[,2]), "min_size_av" = mean(data_table[,3]), "No_gen" = mean(data_table[,4]))
print(data_table_av)
