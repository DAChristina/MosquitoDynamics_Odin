# Lambda #######################################################################
dat_V <- read.csv("Output_V_lambda_Angambiae.csv")
dat_Pos <- read.csv("Output_Pos_lambda_Angambiae.csv")
dat_Prev <- read.csv("Output_Prev_lambda_Angambiae.csv")

# Recall timesteps & lambda
timesteps <- seq(0, 1000, by = 1)
lambda <- seq(0, 100, by = 0.1) # Assume biting rates can increase tenfolds outdoor
# lambda value used = 1 (10 per-months in TRANSFIL, 1 per-cycle)
length(lambda)
lambda[1] # lambda[11] = 1

# I wanna use:
# lambda[1] = 0 bite
# lambda[11] = 1 bite
# lambda[51] = 5 bites
# lambda[101] = 10 bites
# lambda[1001] = 100 bites

par(mfrow = c(1,2))
# Infectives
plot(timesteps, dat_Prev[,101],
     xlab = "timesteps",
     ylab = "Proportion of infective mosquitoes",
     main = "The proportion of infective mosquitoes on different biting rates",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Prev[,101]), paste("bites = 10"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)
lines(timesteps, dat_Prev[,51], col = "black")
text(200, tail(dat_Prev[, 51], 1), paste("bites = 5"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)
lines(timesteps, dat_Prev[,11], col = "black")
text(200, tail(dat_Prev[, 11], 1), paste("bites = 1"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

# Positives
plot(timesteps, dat_Pos[,101],
     xlab = "timesteps",
     ylab = "Proportion of positive mosquitoes",
     main = "The proportion of positive mosquitoes on different biting rates",
     xlim = c(0,400),
     # ylim = c(.001,.007),
     type = "l", col = "black")
text(200, tail(dat_Pos[,101]), paste("bites = 10"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)
lines(timesteps, dat_Pos[,51], col = "black")
text(200, tail(dat_Pos[, 51], 1), paste("bites = 5"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)
lines(timesteps, dat_Pos[,11], col = "black")
text(200, tail(dat_Pos[, 11], 1), paste("bites = 1"),
     adj = c(-0.1, -0.2), col = "black", cex = .8)

par(mfrow = c(1,1))
