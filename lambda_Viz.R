# Data
wd = "C:/Users/dac23/Downloads"
setwd(wd)

Prev <- read.csv("Output_Prev_lambda_Angambiae.csv")
Pos <- read.csv("Output_Pos_lambda_Angambiae.csv")

timesteps <- seq(0, 1000, by = 1)
lambda <- seq(0.01, 90, by = 0.1) # Assume biting rates can increase tenfolds outdoor
length(lambda)

plot(timesteps, Prev[,900], # 90 bites (89.91)
     xlab = "Proportion of LF in human population (with microfilaremia)",
     ylab = "Proportion of infective mosquitoes",
     main = "Proportion of infective mosquitoes",
     # xlim = c(0,1), ylim = c(0,.45),
     type = "l", col = "black")
lines(timesteps, Prev[,600], col = "black") # 60 bites
lines(timesteps, Prev[,300], col = "black") # 30 bites
lines(timesteps, Prev[,200], col = "black") # 20 bites
lines(timesteps, Prev[,100], col = "black") # 10 bites
lines(timesteps, Prev[,10], col = "red") # 1 bite






plot(timesteps, Pos[,900], # 90 bites (89.91)
     xlab = "Proportion of LF in human population (with microfilaremia)",
     ylab = "Proportion of infective mosquitoes",
     main = "Proportion of infective mosquitoes",
     # xlim = c(0,1), ylim = c(0,.45),
     type = "l", col = "black")
lines(timesteps, Pos[,600], col = "black") # 60 bites
lines(timesteps, Pos[,300], col = "black") # 30 bites
lines(timesteps, Pos[,200], col = "black") # 20 bites
lines(timesteps, Pos[,100], col = "black") # 10 bites
lines(timesteps, Pos[,10], col = "black") # 1 bite
