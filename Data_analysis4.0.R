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
    if (pop_new[1]+pop_new[2]+pop_new[3]>=1.5*(N_init_aa+ N_init_Aa+ N_init_AA)){
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
tail(output,5)
number_of_generations <- length(output[,1])-1
p_final <- (0.5*output[number_of_generations,2]+output[number_of_generations,1])/sum(tail(output,2)[1,])
output[1,][1]
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
mut_rate <- 0.0005
t_max <- 10000

# determine how often to run the simulation for each set of parameters
no_replicates <- 1000
# initialize data table - where to collect the results
# run the simulation across all chosen parameters
# loop over switch generations
migrants <- 0:10
iteration_values <- sort(rep(migrants,no_replicates))
data_table <- data.frame(matrix(NA,nrow = no_replicates*length(migrants)*length(fitnessAa),ncol = 8))
for(i in 1:length(iteration_values)){
    one_run <- simulate_pop_HW(N_aa,N_Aa, N_AA, fitnessaa, fitnessAa, fitnessAA, iteration_values[i], mut_rate, t_max)
    # determine total population sizes
    total_size <- head(one_run[,1],-1) + head(one_run[,2],-1) + head(one_run[,3],-1)
    # number of generation at which 1.5*n0 or extinction is achieved
    number_of_generations <- length(one_run[,1])-1
    # determine final pop size
    final_pop_size <- sum(tail(one_run,2)[1,])    
    # determine minimum population size
    min_size <- min(total_size)
    # determine (first) generation at which this population size occurred
    min_gen <- as.numeric(which(total_size==min_size)[1])
    # determine maximal pop size
    max_pop <- max(total_size)
    # determine final allele frequencies
    p_final <- (0.5*one_run[number_of_generations,2]+one_run[number_of_generations,1])/final_pop_size

    # determine success
    success <- tail(one_run[,1],1)
    # enter the data into the table
    data_table[i,] <- c(iteration_values[i],fitnessAa,number_of_generations,min_size,min_gen,final_pop_size,p_final,success) # note that we add the varying parameters (decay rate and selection coefficient) to the table too
}
colnames(data_table) <-c("migration_per_gen","fitnessAa","number_of_generations","min_size","min_gen","final_pop_size","p_final","success")
data_table$success <- factor(data_table$success, level = 1:3, labels = c("extinction","medium","reached_max"))

setwd("C:/Users/flori/OneDrive/Dokumente/Studium/THEE researchpractical/TEE_Shadowvalley")
write.csv(data_table,"Data_29_10_2024.csv")


hist(data_table$final_pop_size)
hist(data_table$final_pop_size[data_table$migration_per_gen == 0])
hist(data_table$final_pop_size[data_table$migration_per_gen == 1])
hist(data_table$final_pop_size[data_table$migration_per_gen == 2])
hist(data_table$final_pop_size[data_table$migration_per_gen == 3])
hist(data_table$final_pop_size[data_table$migration_per_gen == 4])
hist(data_table$final_pop_size[data_table$migration_per_gen == 5])
hist(data_table$final_pop_size[data_table$migration_per_gen == 6])
hist(data_table$final_pop_size[data_table$migration_per_gen == 7])
hist(data_table$final_pop_size[data_table$migration_per_gen == 8])
hist(data_table$final_pop_size[data_table$migration_per_gen == 9])
hist(data_table$final_pop_size[data_table$migration_per_gen == 10])

hist(data_table$p_final)
hist(data_table$p_final[data_table$migration_per_gen == 0])
hist(data_table$p_final[data_table$migration_per_gen == 1])
hist(data_table$p_final[data_table$migration_per_gen == 2])
hist(data_table$p_final[data_table$migration_per_gen == 3])
hist(data_table$p_final[data_table$migration_per_gen == 4])
hist(data_table$p_final[data_table$migration_per_gen == 5])
hist(data_table$p_final[data_table$migration_per_gen == 6])
hist(data_table$p_final[data_table$migration_per_gen == 7])
hist(data_table$p_final[data_table$migration_per_gen == 8])
hist(data_table$p_final[data_table$migration_per_gen == 9])
hist(data_table$p_final[data_table$migration_per_gen == 10])


mig <- matrix(NA,nrow = 3, ncol = 11)
for (i in 0:(ncol(mig)-1)){
  mig[,i+1] <- c(sum(data_table$success[data_table$migration_per_gen == i] == "extinction"),
               sum(data_table$success[data_table$migration_per_gen == i] == "medium"),
               sum(data_table$success[data_table$migration_per_gen == i] == "reached_max"))
}
par(mfrow = c(2,1))
barplot(mig, col = c("red","orange","green"), names.arg = 0:10)
barplot(mig[3,], names.arg = 0:10)

sum(data_table$success[data_table$migration_per_gen == 10] =="extinction")
gen_for_success3 <- data_table$number_of_generations[data_table$success== "reached_max"&data_table$migration_per_gen==3]
gen_for_success4 <- data_table$number_of_generations[data_table$success== "reached_max"&data_table$migration_per_gen==4]
gen_for_success5 <- data_table$number_of_generations[data_table$success== "reached_max"&data_table$migration_per_gen==5]
gen_for_success6 <- data_table$number_of_generations[data_table$success== "reached_max"&data_table$migration_per_gen==6]
gen_for_success7 <- data_table$number_of_generations[data_table$success== "reached_max"&data_table$migration_per_gen==7]
boxplot(gen_for_success3,gen_for_success4,gen_for_success5,gen_for_success6,gen_for_success7, names = 3:7)
