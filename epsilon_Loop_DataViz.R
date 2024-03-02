# Epsilon ######################################################################
dat_V <- read.csv("Output_V_epsilon_Angambiae.csv")
dat_Pos <- read.csv("Output_Pos_epsilon_Angambiae.csv")
dat_Prev <- read.csv("Output_Prev_epsilon_Angambiae.csv")

# Recall timesteps & epsilon
timesteps <- seq(0, 1000, by = 1)
epsilon <- seq(0.1875, 0.5, by = 0.01) # 1/(14/3) = 0.2142857;
# we wanna test 1/(7/3) = 0.4285714 to 1/(16/3) = 0.1875
length(epsilon)
epsilon[1] # epsilon[1] = 0.1875 = 1/(16/3)

# Too lazy to calculate the numbers ############################################
calcEpsilon <- function(day){
  calcDecimal = 1/(day/3)
  epsilon = seq(0.1875, 0.5, by = 0.01)
  
  for (i in seq_along(epsilon)){
    if (abs(calcDecimal - epsilon[i]) <= 0.01) {
      epsilonVal <- i
      break
    }
  }
  
  # return(calcDecimal)
  return(list(calcDecimal, epsilonVal))
}

calcEpsilon(7) # if day = 7, 1/(day/3) = 0.4285714, then epsilon[25]
epsilon[25] # trial epsilon[25]

calcEpsilon(12)
epsilon[7]
tail(dat_Prev[,7])
tail(dat_Pos[,7])

calcEpsilon(13)
epsilon[5]
tail(dat_Prev[,5])
tail(dat_Pos[,5])

calcEpsilon(14)
epsilon[3]
tail(dat_Prev[,3])
tail(dat_Pos[,3])

calcEpsilon(15)
epsilon[2]
tail(dat_Prev[,2])
tail(dat_Pos[,2])

calcEpsilon(16)
epsilon[1]


ress_vec <- list()
for (i in 7:14) {
  trial <- i
  ress <- calcEpsilon(i)
  
  ress_vec[[i-6]] <- list(res = ress, day = i) # save the value of ress_vec to the i-th epsilon
  
  print(paste("Dayy", trial))
  print(ress)
}

# PLOT! ########################################################################

par(mfrow = c(1,2))
# Infectives
plot(timesteps, dat_Prev[,2],
     xlab = "timesteps",
     ylab = "Proportion of infective mosquitoes",
     main = "The proportion of infective mosquitoes on different EIPs",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Prev[,2]), paste(15, " days"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

for (j in 1:length(ress_vec)) {
  lines(timesteps, dat_Prev[, ress_vec[[j]]$res[[2]]], col = "black")
  text(200, tail(dat_Prev[, ress_vec[[j]]$res[[2]]], 1), paste(ress_vec[[j]]$day, " days"),
       adj = c(-0.1, -0.2), col = "black", cex = .8)
}

# Positives
plot(timesteps, dat_Pos[,2],
     xlab = "timesteps",
     ylab = "Proportion of positive mosquitoes",
     main = "The proportion of positive mosquitoes on different EIPs",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Pos[,2]), paste(15, " days"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

for (j in 1:length(ress_vec)) {
  lines(timesteps, dat_Pos[, ress_vec[[j]]$res[[2]]], col = "black")
  text(200, tail(dat_Pos[, ress_vec[[j]]$res[[2]]], 1), paste(ress_vec[[j]]$day, " days"),
       adj = c(-0.1, -0.2), col = "black", cex = .8)
}

par(mfrow = c(1,1))

# Parous Adults
plot(timesteps, dat_V[,1],
     xlab = "timesteps",
     ylab = "Proportion of parous adult mosquitoes",
     main = "The proportion of parous adult mosquitoes on different EIPs",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
# text(200, tail(dat_V[,1]), paste(16, " days"), adj = c(-0.1, -0.2), col = "black")

for (j in 1:length(ress_vec)) {
  lines(timesteps, dat_V[, ress_vec[[j]]$res[[2]]], col = "black")
  # text(200, tail(dat_V[, ress_vec[[j]]$res[[2]]], 1), paste(ress_vec[[j]]$day, " days"), adj = c(-0.1, -0.2), col = "black")
}
