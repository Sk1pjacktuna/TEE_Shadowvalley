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
t_max = 1000
output <- simulate_pop_HW(1000,0,0,0.9,0.9,1.1, 5,0.01,t_max)
tail(output,2)
output[length(output[,1])-1]
output[1001,2]
as.numeric(output[,1][1])+as.numeric(output[,2][1])
output[,1]
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
t_max <- 10000

# determine how often to run the simulation for each set of parameters
no_replicates <- 100
# initialize data table - where to collect the results
# run the simulation across all chosen parameters
# loop over switch generations
migrants <- 0:10
heterozygote <- c(0.9,1,1.1)
iteration_values <- sort(rep(migrants,no_replicates))
data_table <- data.frame(matrix(NA,nrow = no_replicates*length(migrants)*length(fitnessAa),ncol = 3))
zahl <- 0
for (fitnessAa in heterozygote){
  for(i in 1:length(iteration_values)){
    one_run <- simulate_pop_HW(N_aa,N_Aa, N_AA, fitnessaa, fitnessAa, fitnessAA, iteration_values[i], mut_rate, t_max)
    # determine total population sizes
    #total_size <- one_run[,1] + one_run[,2] + one_run[,3]
    #number of generation at which 1.5*n0 is achieved
    #number_of_generations <- length(one_run[,1])-1
    # determine minimum population size
    #min_size <- min(total_size)
    # determine (first) generation at which this population size occurred
    #min_gen <- as.numeric(which(total_size==min_size)[1])
    # determine maximal pop size
    #max_pop <- max(total_size)
    # determine final allele frequencies
    #p_final <- (0.5*one_run[number_of_generations-1,1]+one_run[number_of_generations,2])/sum(one_run[number_of_generations,1:3])
    # determine success
    success <- tail(one_run[,1],1)
    # enter the data into the table
    data_table[i+zahl,] <- c(iteration_values[i],fitnessAa,success) # note that we add the varying parameters (decay rate and selection coefficient) to the table too
  }
  zahl <- zahl+length(iteration_values)
}
data_matrix <- matrix(NA,nrow = 3,ncol = 11)
for (i in 1:ncol(data_matrix)){
  data_matrix[,i] <- c(sum(data_table$X3[data_table$X1 == i-1&data_table$X2 ==0.9] == 1),sum(data_table$X3[data_table$X1 == i-1&data_table$X2 ==0.9] == 2),sum(data_table$X3[data_table$X1 == i-1&data_table$X2 ==0.9] == 3))
}
install.packages("viridis")
library("viridis")
colors <- magma(3,begin = 0.2,end = 0.8)
barplot(data_matrix,names.arg = 0:10, xlab = "Average Migrants per generation", col = colors, legend = c("extinction","happy ever after","reached max"), ylab = "proportion in %")


analysis <- data.frame(migration = data_table[,1], generation_with_min_pop_size = data_table[,2], min_pop_size = data_table[,3],simulated_gen = data_table[,4],maximal_pop_size = data_table[,5],success = data_table[,6] )

analysis$success <- factor(analysis$success, 1:3, c("extinction","medium","reached_max"))

mean(analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration==5])
mean(analysis$simulated_gen[analysis$success == "extinction"])
mean(analysis$simulated_gen[analysis$success == "medium"])

boxplot(analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 0],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 1],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 2],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 3],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 4],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 5],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 6],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 7],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 8],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 9],analysis$simulated_gen[analysis$success == "reached_max"&analysis$migration== 10])

boxplot(analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 0],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 1],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 2],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 3],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 4],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 5],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 6],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 7],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 8],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 9],analysis$simulated_gen[analysis$success == "extinction"&analysis$migration== 10])


filename ="analysis.csv"
write.csv(analysis,file = filename)
boxplot(data_table[,4][data_table[,1]==0], data_table[,4][data_table[,1]==1],data_table[,4][data_table[,1]==2],data_table[,4][data_table[,1]==3],data_table[,4][data_table[,1]==4],data_table[,4][data_table[,1]==5],data_table[,4][data_table[,1]==6],data_table[,4][data_table[,1]==7],data_table[,4][data_table[,1]==8],data_table[,4][data_table[,1]==9],data_table[,4][data_table[,1]==10], main = "numbers of generations run")

boxplot(data_table[,3][data_table[,1]==0],data_table[,3][data_table[,1]==1],data_table[,3][data_table[,1]==2],data_table[,3][data_table[,1]==3],data_table[,3][data_table[,1]==4],data_table[,3][data_table[,1]==5],data_table[,3][data_table[,1]==6],data_table[,3][data_table[,1]==7],data_table[,3][data_table[,1]==8],data_table[,3][data_table[,1]==9],data_table[,3][data_table[,1]==10],main ="minimal pop size for each avgmigrants")

boxplot(data_table[,5][data_table[,1]==0], data_table[,5][data_table[,1]==1],data_table[,5][data_table[,1]==2],data_table[,5][data_table[,1]==3],data_table[,5][data_table[,1]==4],data_table[,5][data_table[,1]==5],data_table[,5][data_table[,1]==6],data_table[,5][data_table[,1]==7],data_table[,5][data_table[,1]==8],data_table[,5][data_table[,1]==9],data_table[,5][data_table[,1]==10],main = "maximal pop size")

