# K_loop data visualisation

head(V_mtx)
tail(V_mtx[,2679])
K_values[2679]

# LOG Plot for total (parous) mosquitoes
plot(timesteps,log10(V_mtx[,2679]), type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "Time (day)",
     ylab = "log10(total parous mosquitoes)")
lines(timesteps,log10(V_mtx[,269]))  # V/H = 10
lines(timesteps,log10(V_mtx[,28]))   # V/H = 1
lines(timesteps,log10(V_mtx[,4]))    # V/H = 0.1

legend("bottomright",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 1)

# LOG Plot for total positive (E_v + I_v) mosquitoes
plot(timesteps,log10(Pos_mtx[,2679]), type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "Time (day)",
     ylab = "log10(total parous mosquitoes)")
lines(timesteps,log10(Pos_mtx[,269]))  # V/H = 10
lines(timesteps,log10(Pos_mtx[,28]))   # V/H = 1
lines(timesteps,log10(Pos_mtx[,4]))    # V/H = 0.1

legend("bottomright",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 1)

# Plot for total positive (E_v + I_v) mosquitoes
plot(timesteps,Pos_mtx[,2679], type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "Time (day)",
     ylab = "Proportion of positive mosquitoes)")
lines(timesteps,Pos_mtx[,269])  # V/H = 10
lines(timesteps,Pos_mtx[,28])   # V/H = 1
lines(timesteps,Pos_mtx[,4])    # V/H = 0.1

legend("bottomright",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 1)

# Plot for total infective (I_v) mosquitoes
plot(timesteps,Prev_mtx[,2679], type = "l", # V/H = 100
     xlim = c(0,50),
     xlab = "Time (day)",
     ylab = "Proportion of infective mosquitoes)")
lines(timesteps,Prev_mtx[,269])  # V/H = 10
lines(timesteps,Prev_mtx[,28])   # V/H = 1
lines(timesteps,Prev_mtx[,4])    # V/H = 0.1

legend("bottomright",
       legend = c(paste0("K = ", K_values[2679]), paste0("K = ", K_values[269]), paste0("K = ", K_values[28]), paste0("K = ", K_values[4])),
       cex = 1)