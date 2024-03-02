# Gamma ########################################################################
dat_V <- read.csv("Output_V_gamma_Angambiae.csv")
dat_Pos <- read.csv("Output_Pos_gamma_Angambiae.csv")
dat_Prev <- read.csv("Output_Prev_gamma_Angambiae.csv")

# Recall timesteps & gamma
timesteps <- seq(0, 1000, by = 1)
sg <- seq(9.82, 18, by = 0.01) # Density-dependent used = 13.25
# What if I wanna test 10 values of sg
length(sg)
sg[1] # sg[1] = 9.82

par(mfrow = c(1,2))
# Infectives
plot(timesteps, dat_Prev[,1],
     xlab = "timesteps",
     ylab = "Proportion of infective mosquitoes",
     main = "The proportion of infective mosquitoes on different sigma",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Prev[,2]), paste("sigma = vary"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

for (i in 1:length(sg)){
  lines(timesteps, dat_Prev[, i], col = "black")
}

# Positives
plot(timesteps, dat_Pos[,1],
     xlab = "timesteps",
     ylab = "Proportion of positive mosquitoes",
     main = "The proportion of positive mosquitoes on different sigma",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Pos[,2]), paste("sigma = vary"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

for (i in 1:length(sg)){
  lines(timesteps, dat_Pos[, i], col = "black")
}

par(mfrow = c(1,1))

# Parous Adults
plot(timesteps, dat_V[,1],
     xlab = "timesteps",
     ylab = "Proportion of parous adult mosquitoes",
     main = "The proportion of parous adult mosquitoes on different sigma",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
# text(200, tail(dat_V[,1]), paste(16, " days"), adj = c(-0.1, -0.2), col = "black")

for (j in 1:length(sg)) {
  lines(timesteps, dat_V[, sg[j]], col = "black")
  # text(200, tail(dat_V[, ress_vec[[j]]$res[[2]]], 1), paste(ress_vec[[j]]$day, " days"), adj = c(-0.1, -0.2), col = "black")
}
