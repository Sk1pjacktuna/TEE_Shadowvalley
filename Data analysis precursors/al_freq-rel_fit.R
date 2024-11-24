

al_freq <- seq(0,1,length = 20)

rel_fit <- 1 - (1 - al_freq ^ 2) * (0.1) - al_freq ^ 2 * (-0.1) 

plot(al_freq,rel_fit, type = "l")
abline(h = 1)
abline(v = 0.707)


