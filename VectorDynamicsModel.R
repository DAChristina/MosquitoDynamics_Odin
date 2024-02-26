if (!require(odin, quietly=T)){
  install.packages("odin")
  library(odin)
}

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
  
  # 1. PARAMETERS ##############################################################
  # Gompertz mortality rate have already in cycle (Clements & Paterson, 1981)
  g1 <- user()
  g2 <- user()
  # With Survival(i) = (exp(-g1/g2*(exp(i*g2)-1))
  
  # Population Dynamics, Eggs -> Larvae -> Mature (White et al., 2011)
  beta <- 21.19*3 # Egg deposition per-capita, per-day * 3 days for 1 gonotrophic cycles
  mu0 <- .034 # Per-capita daily mortality rate of Eggs & early instar larvae ^3 days for 1 gonotrophic cycles
  mu1 <- .035 # Per-capita daily mortality rate of Late instar larvae ^3 days for 1 gonotrophic cycles
  mu2 <- .25 # Per-capita daily mortality rate of pupae
  gE <- 1/(6.64/3) # Time required for growth of Eggs to early instar larvae (in 1 cycle)
  gL <- 1/((3.72+0.64)/3)  # Time required for growth of early instar larvae to adult mosquitoes (in 1 cycle)
  
  K <- 200 # Saturation coefficient
  sg <- 13.25 # Effects of density-dependence on late instars (L) relative to early instars (E)
  
  ## S,E,I are arrays
  lambda <- 10/10 # biting rate of mosquitoes per cycle (source: TRANSFIL, 1 month of TRANSFIL has 10 cycles)
  InfecMosq <- 0.37 # Vector competence, the proportion of mosquitoes which pick up the infection when biting an infective host (source: TRANSFIL)
  epsilon <- 1/(14/3) # incubation rate of LF in mosquitoes (per-cycle)
  InfHuman <- user() # Proportion of infected humans in the population with detectable microfilariae
  
  # 2. INITIAL VALUES ##########################################################
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
  
  output(prev) <- I_v_tot/V_tot
  output(pos) <- (E_v_tot+I_v_tot)/V_tot
  
  config(base) <- "transition"
})

# Trial2 Various cycle_width AND InfHuman values: ##############################
cycle_width_values <- seq(3, 30, by = 3)
InfHuman_values <- c(seq(0, 0.02, by = 0.001), seq(0.02, 1, by = 0.01))  # Separate InfHuman into 2 prevalence ranges; below WHO threshold and above
Gompz_pars1 <- c(.356, .339) # 1 = An. gambiae, 2 = An. arabiensis
Gompz_pars2 <- c(.097, .225) # 1 = An. gambiae, 2 = An. arabiensis

# Vector dimension storage
S_v_loop <- numeric(length(InfHuman_values))
E_v_loop <- numeric(length(InfHuman_values))
I_v_loop <- numeric(length(InfHuman_values))

# Additional coz' I'm curious about the result
V_loop <- numeric(length(InfHuman_values))
Prev_loop <- numeric(length(InfHuman_values)) # Iv/V OR proportion of all infective mosquitoes
Pos_loop <- numeric(length(InfHuman_values)) # (Ev+Iv)/V OR proportion of all positive mosquitoes

# Trial loop function
for (i in seq_along(InfHuman_values)) {
  pars <- list(cycle_width = cycle_width_values,
               InfHuman = InfHuman_values[i],
               g1 = Gompz_pars1[1], # choose[1] for An. gambiae, [2] for An. arabiensis
               g2 = Gompz_pars2[1]) # choose[1] for An. gambiae, [2] for An. arabiensis)
  
  mod <- transition$new(user = pars)
  timesteps <- seq(0, 2000, by = 1)
  y <- mod$run(timesteps)
  
  S_v_loop[i] <- tail(y[,"S_v_tot"], 1) # All susceptibles from 10 gonotrophic cycle, given InfHuman
  E_v_loop[i] <- tail(y[,"E_v_tot"], 1) # All exposed from 10 gonotrophic cycle, given InfHuman
  I_v_loop[i] <- tail(y[,"I_v_tot"], 1) # All infectives from 10 gonotrophic cycle, given InfHuman
  
  V_loop[i] <- tail(y[,"V_tot"], 1) # total parous mosquitoes, given InfHuman
  
  Prev_loop[i] <- tail(y[,"prev"], 1) # Infective only
  Pos_loop[i] <- tail(y[,"pos"], 1) # All positives
}


# dataframe storage
Output_InfHuman_Angambiae <- data.frame(
  # timesteps not required because it gives repetitions in each column
  InfHuman_values = InfHuman_values,
  S_v_loop = S_v_loop,
  E_v_loop = E_v_loop,
  I_v_loop = I_v_loop,
  
  # Additional coz' I'm curious about the result
  V_loop = V_loop,
  Prev_loop = Prev_loop, # Iv/V OR proportion of all infective mosquitoes
  Pos_loop = Pos_loop # All positives
)

# Output_InfHuman_Angambiae

write.csv(Output_InfHuman_Angambiae, file = "Output_InfHuman_Angambiae.csv", row.names = FALSE)
