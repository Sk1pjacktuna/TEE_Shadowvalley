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
  a_to_A_mut <- min(rpois(1,a_tot*mut_rate),a_tot)
  A_to_a_mut <- min(rpois(1,A_tot*mut_rate),A_tot)
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
  success <- 2
  # run the simulation until generation t_max
  for (i in 1:t_max+1) {
    #define fitness of a and A for each iteration
    # redefine the current population one generation later
    pop_new <- simulate_one_gen_hardy_weinberg(pop_new[1],pop_new[2], pop_new[3], fitnessaa, fitnessAa, fitnessAA,avgmigrants, mut_rate)
   # add the new population sizes to the output vector
    pop_vector <- rbind(pop_vector,pop_new) 
    # condition to stop the simulation before t_max: either the population exceeds 1.5 times the original population size, or it goes extinct
    if (pop_new[1]+pop_new[2]+pop_new[3]>=5*(N_init_aa+ N_init_Aa+ N_init_AA)){
      success <- 3
      break
    }else if (pop_new[1]+pop_new[2]+pop_new[3]<=2*avgmigrants){
      success <- 1
      break
    }
  }
  pop_vector <- rbind(pop_vector,success)
  # define the row and column names of the output vector
  rownames(pop_vector) <- (0:(t_max+1))[1:length(pop_vector[,1])] # note that the vector has to be cut if the simulation stopped early
  colnames(pop_vector) <- c("aa","Aa","AA")
  # return the result
  return(pop_vector)	
}
output <- simulate_pop_HW(1000,0,0,0.9,0.9,1.1, 5,0.01,t_max)
tail(output[,1],1)
as.numeric(output[,1][1])+as.numeric(output[,2][1])
head(output[,1],-1)
View(output)
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
N_AA <- 0
fitnessaa <- 0.9
fitnessAa <- 0.9
fitnessAA <- 1.1
mut_rate <- 0.01
t_max <- 1000

# determine how often to run the simulation for each set of parameters
no_replicates <- 100

# initialize data table - where to collect the results
# run the simulation across all chosen parameters
# loop over switch generations
migrants <- 0:10
data_table <- c()
for(avgmigrants in migrants){
    # reset counter
    i<-1
    repeat {
      # increase counter by one
      i<-i+1
      # run the simulation once
      one_run <- simulate_pop_HW(N_aa,N_Aa, N_AA, fitnessaa, fitnessAa, fitnessAA, avgmigrants, mut_rate, t_max)
      # determine total population sizes
      total_size <- head(one_run[,1],-1) + head(one_run[,2],-1) + head(one_run[,3],-1)
      #number of generation at which 1.5*n0 is achieved
      number_of_generations <- length(one_run[,1])-1
      # determine minimum population size
      min_size <- min(total_size)
      # determine (first) generation at which this population size occurred
      min_gen <- as.numeric(which(total_size==min_size)[1])
      # determine maximal pop size
      max_pop <- max(total_size)
      # determine success
      success <- tail(one_run[,1],1)
      # enter the data into the table
      data_table <- rbind(data_table,c(avgmigrants,min_gen,min_size,number_of_generations,max_pop,success)) # note that we add the varying parameters (decay rate and selection coefficient) to the table too
      # stop the repeated computation after no_replicates times
      if(i>no_replicates) break
    }
}

analysis <- data.frame(migration = data_table[,1], generation_with_min_pop_size = data_table[,2], min_pop_size = data_table[,3],simulated_gen = data_table[,4],maximal_pop_size = data_table[,5],success = data_table[,6] )

analysis$success <- factor(analysis$success, 1:3, c("extinction","medium","reached_max"))

mean(analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration==5])
mean(analysis$simulated_gen[analysis$success == "extinction"])
mean(analysis$simulated_gen[analysis$success == "medium"])

boxplot(analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 0],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 1],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 2],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 3],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 4],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 5],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 6],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 7],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 8],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 9],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 10])

boxplot(analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 0],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 1],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 2],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 3],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 4],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 5],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 6],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 7],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 8],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 9],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 10])

boxplot(data_table[,4][data_table[,1]==0], data_table[,4][data_table[,1]==1],data_table[,4][data_table[,1]==2],data_table[,4][data_table[,1]==3],data_table[,4][data_table[,1]==4],data_table[,4][data_table[,1]==5],data_table[,4][data_table[,1]==6],data_table[,4][data_table[,1]==7],data_table[,4][data_table[,1]==8],data_table[,4][data_table[,1]==9],data_table[,4][data_table[,1]==10], main = "numbers of generations run")

boxplot(data_table[,3][data_table[,1]==0],data_table[,3][data_table[,1]==1],data_table[,3][data_table[,1]==2],data_table[,3][data_table[,1]==3],data_table[,3][data_table[,1]==4],data_table[,3][data_table[,1]==5],data_table[,3][data_table[,1]==6],data_table[,3][data_table[,1]==7],data_table[,3][data_table[,1]==8],data_table[,3][data_table[,1]==9],data_table[,3][data_table[,1]==10],main ="minimal pop size for each avgmigrants")

boxplot(data_table[,5][data_table[,1]==0], data_table[,5][data_table[,1]==1],data_table[,5][data_table[,1]==2],data_table[,5][data_table[,1]==3],data_table[,5][data_table[,1]==4],data_table[,5][data_table[,1]==5],data_table[,5][data_table[,1]==6],data_table[,5][data_table[,1]==7],data_table[,5][data_table[,1]==8],data_table[,5][data_table[,1]==9],data_table[,5][data_table[,1]==10],main = "maximal pop size")

install.packages("ggplot2")
library(ggplot2)
datenggplot <- read.csv("analysis.csv")
