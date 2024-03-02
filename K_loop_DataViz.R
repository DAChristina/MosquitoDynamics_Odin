# K ############################################################################
dat_V <- read.csv("Output_V_K_Angambiae.csv")
dat_Pos <- read.csv("Output_Pos_K_Angambiae.csv")
dat_Prev <- read.csv("Output_Prev_K_Angambiae.csv")

K_values <- seq(0, 268000, by = 100)


head(dat_V)
tail(dat_V[,2679])
K_values[2679]

# LOG Plot for total (parous) mosquitoes
plot(timesteps,log10(dat_V[,2679]), type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "timesteps",
     ylab = "log10(total parous mosquitoes)",
     main = "The log proportion of total parous mosquitoes on various K")
lines(timesteps,log10(dat_V[,269]))  # V/H = 10
lines(timesteps,log10(dat_V[,28]))   # V/H = 1
lines(timesteps,log10(dat_V[,4]))    # V/H = 0.1

legend("topleft",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 0.8, bty = "n")

par(mfrow = c(1,2))
# LOG Plot for total positive (E_v + I_v) mosquitoes
plot(timesteps,dat_Prev[,2679], type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "timesteps",
     ylab = "Proportion of infective mosquitoes",
     main = "The proportion of infective mosquitoes on various K")
lines(timesteps,dat_Prev[,269])  # V/H = 10
lines(timesteps,dat_Prev[,28])   # V/H = 1
lines(timesteps,dat_Prev[,4])    # V/H = 0.1

legend("bottomright",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 0.8, bty = "n")

# Plot for total positive (E_v + I_v) mosquitoes
plot(timesteps,dat_Pos[,2679], type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "timesteps",
     ylab = "Proportion of positive mosquitoes)",
     main = "The proportion of positive mosquitoes on various K")
lines(timesteps,dat_Pos[,269])  # V/H = 10
lines(timesteps,dat_Pos[,28])   # V/H = 1
lines(timesteps,dat_Pos[,4])    # V/H = 0.1

legend("bottomright",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 0.8, bty = "n")

par(mfrow = c(1,1))

# Plot for total infective (I_v) mosquitoes
plot(timesteps,Predat_V[,2679], type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "timesteps",
     ylab = "Proportion of infective mosquitoes)")
lines(timesteps,Predat_V[,269])  # V/H = 10
lines(timesteps,Predat_V[,28])   # V/H = 1
lines(timesteps,Predat_V[,4])    # V/H = 0.1

legend("bottomright",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 0.8, bty = "n")
