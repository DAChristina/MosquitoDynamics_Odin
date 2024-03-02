# mu0 #######################################################################
dat_V <- read.csv("Output_V_mu0_Angambiae.csv")
dat_Pos <- read.csv("Output_Pos_mu0_Angambiae.csv")
dat_Prev <- read.csv("Output_Prev_mu0_Angambiae.csv")

# Recall timesteps & mu0
timesteps <- seq(0, 1000, by = 1)
mu0 <- seq(0.025, 0.045, by = 0.001) # from White et al. (2011)
# mu0 value used = 1 (10 per-months in TRANSFIL, 1 per-cycle)
length(mu0)
mu0[1] # mu0[11] = 1

par(mfrow = c(1,2))
# Infectives
plot(timesteps, dat_Prev[,1],
     xlab = "timesteps",
     ylab = "Proportion of infective mosquitoes",
     main = "The proportion of infective mosquitoes on different mortality rate of E",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Prev[,1]), paste("Mortality rate of E = vary"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

for (i in 1:length(mu0)){
  lines(timesteps, dat_Prev[, i], col = "black")
}

# Positives
plot(timesteps, dat_Pos[,1],
     xlab = "timesteps",
     ylab = "Proportion of positive mosquitoes",
     main = "The proportion of positive mosquitoes on different mortality rate of E",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Pos[,2]), paste("Mortality rate of E = vary"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

for (i in 1:length(mu0)){
  lines(timesteps, dat_Pos[, i], col = "black")
}

par(mfrow = c(1,1))

# Parous Adults
plot(timesteps, dat_V[,1],
     xlab = "timesteps",
     ylab = "Proportion of parous adult mosquitoes",
     main = "The proportion of parous adult mosquitoes on different mortality rate of E",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
# text(200, tail(dat_V[,1]), paste(16, " days"), adj = c(-0.1, -0.2), col = "black")

for (j in 1:length(mu0)) {
  lines(timesteps, dat_V[, mu0[j]], col = "black")
  # text(200, tail(dat_V[, ress_vec[[j]]$res[[2]]], 1), paste(ress_vec[[j]]$day, " days"), adj = c(-0.1, -0.2), col = "black")
}


for (j in 1:length(mu0)) {
  print(tail(dat_V[, mu0[j]],1))
}

