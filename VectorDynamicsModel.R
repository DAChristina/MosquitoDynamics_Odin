transition <- odin::odin({
  N_cycle <- length(cycle_width)
  # Nulliparous means day 0, mosquitoes emerged from aquatic stages
  cycle_width[] <- user() # First gonotrophic cycle, 3 days
  
  cycle_rate[1:(N_cycle - 1)] <- 1 / cycle_width[i]
  cycle_rate[N_cycle] <- 0
  
  # den = density of the population in each age group
  den[1] <- 1 / (1 + cycle_rate[1] / (1/beta))
  ## to work out the % of the population in each cycle group
  den[2:N_cycle] <- cycle_rate[i - 1] * den[i - 1] / (cycle_rate[i] + (1/beta))
  
  # 1. PARAMETERS ################################################################
  # Gompertz mortality rate have already in cycle (Clements & Paterson, 1981)
  g1 <- user()
  g2 <- user()
  # With Survival(i) = (exp(-g1/g2*(exp(i*g2)-1))
  
  # Population Dynamics, Eggs -> Larvae -> Mature (White et al., 2011)
  beta <- 21.19*3 # Egg deposition per-capita, per-day * 3 days for 1 gonotrophic cycles
  mu0 <- user() # Per-capita daily mortality rate of Eggs & early instar larvae ^3 days for 1 gonotrophic cycles
  mu1 <- .035 # Per-capita daily mortality rate of Late instar larvae ^3 days for 1 gonotrophic cycles
  mu2 <- .25 # Per-capita daily mortality rate of pupae
  gE <- 1/(6.64/3) # Time required for growth of Eggs to early instar larvae (in 1 cycle)
  gL <- 1/((3.72+0.64)/3)  # Time required for growth of early instar larvae to adult mosquitoes (in 1 cycle)
  
  K <- user() # Saturation coefficient
  sg <- 13.25 # Effects of density-dependence on late instars (L) relative to early instars (E)
  
  ## S,E,I are arrays
  lambda <- user() # biting rate of mosquitoes per cycle (source: TRANSFIL, 1 month of TRANSFIL has 10 cycles)
  InfecMosq <- 0.37 # Vector competence, the proportion of mosquitoes which pick up the infection when biting an infective host (source: TRANSFIL)
  epsilon <- user() # incubation rate of LF in mosquitoes (per-cycle)
  InfHuman <- user()
  
  # 2. INITIAL VALUES ############################################################
  initial(E) <- 115
  initial(L) <- 60
  initial(N) <- 52
  
  E_v0 <- 0
  I_v0 <- 0
  initial(S_v[1:N_cycle]) <- den[i]*(V_tot -(E_v_tot+I_v_tot))
  initial(E_v[1:N_cycle]) <- den[i]*E_v0
  initial(I_v[1:N_cycle]) <- den[i]*I_v0
  
  # Dimension & values of arrays
  dim(S_v) <- N_cycle
  dim(E_v) <- N_cycle
  dim(I_v) <- N_cycle
  dim(den) <- N_cycle
  dim(cycle_width) <- user()
  dim(cycle_rate) <- N_cycle
  
  # Define mortality rates in 1 cycle (3 days)
  muE <- mu0*3*(1+(E+L)/K)
  muL <- mu1*3*(1+sg*(E+L)/K)
  
  # 3. DERIVATIVES #############################################################
  temp_deriv_E <- beta*V_tot -E*(gE+muE)
  temp_deriv_L <- E*gE -L*(gL+muL)
  
  deriv(E) <- if (temp_deriv_E > 0) temp_deriv_E else 0 # Checking of the result < 0, throw 0
  deriv(L) <- if (temp_deriv_L > 0) temp_deriv_L else 0 # Checking of the result < 0, throw 0
  deriv(N) <- L*gL*(1-mu2)*.5 -N*(g1*exp((0)*3*g2)) # given mortality as age = 1(*3 for 1 cycle), assume g1*exp((0)*3*g2) = baseline mortality rate at day 0
  
  # Insert deriv(N) as the first cycle of susceptible mosquitoes
  # Disease dynamics
  deriv(S_v[1])         <- -S_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*((1+lambda*InfecMosq*InfHuman)) +(N                                                        - cycle_rate[i]*S_v[i]*(exp(-g1/g2*(exp(i*g2)-1))))
  deriv(S_v[2:N_cycle]) <- -S_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*((1+lambda*InfecMosq*InfHuman)) +(cycle_rate[i-1]*S_v[i-1]*(exp(-g1/g2*(exp((i-1)*g2)-1))) - cycle_rate[i]*S_v[i]*(exp(-g1/g2*(exp(i*g2)-1))))
  # dS = [previous state] + [recent state]*(1-exposed-death)
  
  deriv(E_v[1])         <- S_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*(lambda*InfecMosq*InfHuman) -E_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*(1+epsilon) +(                                                         - cycle_rate[i]*E_v[i]*(exp(-g1/g2*(exp(i*g2)-1))))
  deriv(E_v[2:N_cycle]) <- S_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*(lambda*InfecMosq*InfHuman) -E_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*(1+epsilon) +(cycle_rate[i-1]*E_v[i-1]*(exp(-g1/g2*(exp((i-1)*g2)-1))) - cycle_rate[i]*E_v[i]*(exp(-g1/g2*(exp(i*g2)-1))))
  # dE = [previous state] + [recent state]*(1-infected-death)
  
  deriv(I_v[1])         <- E_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*(1-epsilon) -I_v[i]*(exp(-g1/g2*(exp(i*g2)-1))) +(                                                         - cycle_rate[i]*I_v[i]*(exp(-g1/g2*(exp(i*g2)-1))))
  deriv(I_v[2:N_cycle]) <- E_v[i]*(exp(-g1/g2*(exp(i*g2)-1)))*(1-epsilon) -I_v[i]*(exp(-g1/g2*(exp(i*g2)-1))) +(cycle_rate[i-1]*I_v[i-1]*(exp(-g1/g2*(exp((i-1)*g2)-1))) - cycle_rate[i]*I_v[i]*(exp(-g1/g2*(exp(i*g2)-1))))
  # dI = [previous state] + [recent state]*(1-death)
  
  # Cumulative S-E-I
  S_v_tot <- sum(S_v)
  E_v_tot <- sum(E_v)
  I_v_tot <- sum(I_v)
  V_tot <- S_v_tot + E_v_tot + I_v_tot
  
  output(S_v_tot) <- S_v_tot
  output(E_v_tot) <- E_v_tot
  output(I_v_tot) <- I_v_tot
  output(V_tot) <- V_tot
  
  output(prev) <- I_v_tot / V_tot
  output(pos) <- (E_v_tot + I_v_tot) / V_tot
  
  config(base) <- "transition"
})

cycle_width_values <- seq(3, 30, by = 3)
# InfHuman_values <- c(seq(0, 0.02, by = 0.001), seq(0.02, 1, by = 0.01))  # Separate InfHuman into 2 prevalence ranges; below WHO threshold and above
Gompz_pars1 <- c(.356, .339) # 1 = An. gambiae, 2 = An. arabiensis
Gompz_pars2 <- c(.097, .225) # 1 = An. gambiae, 2 = An. arabiensis

# K_values <- seq(0, 268000, by = 100) # Assume V/H = 0.01 to 100,000 (given human pop H = 1,000 in TRANSFIL)
lambda <- seq(0.01, 90, by = 0.1) # Assume biting rates can increase tenfolds outdoor
epsilon <- seq(0.1875, 0.5, by = 0.01) # 1/(14/3) = 0.2142857; we wanna test 1/(7/3) to 1/(16/3)
mu0 <- seq(0.025, 0.045, by = 0.001)
timesteps <- seq(0, 1000, by = 1)

# MATRIX dimension storage (tried as hard as I can to avoid matrix but failed)
E_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))
L_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))
N_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))

S_v_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))
E_v_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))
I_v_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))

V_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))

Prev_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))
Pos_mtx <- matrix(NA, nrow=length(timesteps), ncol=length(mu0))

for (i in seq_along(mu0)) {
  mu <- mu0[i]
  
  pars <- list(cycle_width =cycle_width_values,
               InfHuman =0.01, # Assume 1% human prevalence
               K =267800,
               lambda = 1,
               epsilon = 1/(14/3),
               mu0 = mu, # stored coz' I wanna print the progress
               g1 =Gompz_pars1[1], # choose[1] for An. gambiae, [2] for An. arabiensis
               g2 =Gompz_pars2[1]) # choose[1] for An. gambiae, [2] for An. arabiensis)
  
  mod <- transition$new(user = pars)
  y <- mod$run(timesteps)
  
  # Store the results in the matrix
  E_mtx[, i] <- y[, "E"]
  L_mtx[, i] <- y[, "L"]
  N_mtx[, i] <- y[, "N"]
  
  S_v_mtx[, i] <- y[, "S_v_tot"]
  E_v_mtx[, i] <- y[, "E_v_tot"]
  I_v_mtx[, i] <- y[, "I_v_tot"]
  
  V_mtx[, i] <- y[, "V_tot"]
  
  Prev_mtx[, i] <- y[,"prev"]
  Pos_mtx[, i] <- y[,"pos"]
  
  # SHOW the progresss to avoid the computer's sleeping/logout -_-)
  print(paste("Processing mu0 =", mu))
}

# How to access the result?
# How about save it to *csvs???
write.csv(E_mtx, file="Output_E_mu0_Angambiae.csv", row.names =F)
write.csv(L_mtx, file="Output_L_mu0_Angambiae.csv", row.names =F)
write.csv(N_mtx, file="Output_N_mu0_Angambiae.csv", row.names =F)
write.csv(S_v_mtx, file="Output_S_v_mu0_Angambiae.csv", row.names =F)
write.csv(E_v_mtx, file="Output_E_v_mu0_Angambiae.csv", row.names =F)
write.csv(I_v_mtx, file="Output_I_v_mu0_Angambiae.csv", row.names =F)
write.csv(V_mtx, file="Output_V_mu0_Angambiae.csv", row.names =F)
write.csv(Prev_mtx, file="Output_Prev_mu0_Angambiae.csv", row.names =F)
write.csv(Pos_mtx, file="Output_Pos_mu0_Angambiae.csv", row.names =F)
