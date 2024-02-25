library(tidyverse)

# wd =  # Change the working directory to the location of the *.csv files if downloaded to PC
# setwd(wd)

An_g <- read.csv("Output_InfHuman_Angambiae.csv")
An_a <- read.csv("Output_InfHuman_Anarabiensis.csv")

dat <- left_join(An_g, An_a, by = "InfHuman_values", relationship = "many-to-many")
# x = An. gambiae (red)
# y = An. arabiensis (blue)

# Given all values of InfHuman
par(mfrow = c(1,2))
plot(dat$InfHuman_values, dat$Prev_loop.x, # An. gambiae
     xlab = "Proportion of LF in human population (with microfilaremia)",
     ylab = "Proportion of infective mosquitoes",
     main = "Proportion of infective mosquitoes",
     xlim = c(0,1), ylim = c(0,.45), type = "l", col = "red")
lines(dat$InfHuman_values, dat$Prev_loop.y, col = "blue") # An. arabiensis

plot(dat$InfHuman_values, dat$Pos_loop.x, # An. gambiae
     xlab = "Proportion of LF in human population (with microfilaremia)",
     ylab = "Proportion of positive mosquitoes",
     main = "Proportion of positive mosquitoes",
     xlim = c(0,1), ylim = c(0,.45), type = "l", col = "red")
lines(dat$InfHuman_values, dat$Pos_loop.y, col = "blue") # An. arabiensis
par(mfrow = c(1,1))


# ONLY InfHuman below WHO threshold (either it is .01 or .02)
dat <- dat %>% 
  filter(InfHuman_values <= .020)

par(mfrow = c(1,2))
plot(dat$InfHuman_values, dat$Prev_loop.x, # An. gambiae
     xlab = "Proportion of LF in human population below WHO thresholds",
     ylab = "Proportion of infective mosquitoes",
     main = "Proportion of infective mosquitoes",
     xlim = c(0,.021), ylim = c(0,.014), type = "l", col = "red")
lines(dat$InfHuman_values, dat$Prev_loop.y, col = "blue") # An. arabiensis

plot(dat$InfHuman_values, dat$Pos_loop.x, # An. gambiae
     xlab = "Proportion of LF in human population below WHO thresholds",
     ylab = "Proportion of positive mosquitoes",
     main = "Proportion of positive mosquitoes",
     xlim = c(0,.021), ylim = c(0,.014), type = "l", col = "red")

# Burkina Faso data points:
points(0.00471, 0.033898) # (Human data = Dano (Nablegane 2016 TAS1), Mosquito data = Ouessa)
text(0.00471, 0.033898, "(Dano, Ouessa)", adj = c(-0.1, -0.5), col = "black")

points(0.061, 0.003367) # (Human data = Diebougou (Danko-Tanzou 2017 spot-check FTS 19/316 participants), Mosquito data = Bapla)
text(0.061, 0.003367, "(Diebougou, Bapla)", adj = c(-0.1, -0.5), col = "black")

lines(dat$InfHuman_values, dat$Pos_loop.y, col = "blue") # An. arabiensis
par(mfrow = c(1,1))
