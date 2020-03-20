# Laura's re-working of Kristen's Leaky or Loyal model
# March 20, 2020

require(tidyverse)

# State Variables 
# 
# P = individual plant biomass (tree)
# A = individual fungus biomass
# N_P = Nutrient pool in the plant
# N_S = Nutrient pool in the soil
# 
# Parameters
# 
# g = growth rate of tree
# s = senescence. density independent mortality of tree (loss of leaves, branches, etc)
# d_P = density dependent mortality of tree (as a tree gets bigger, more likely to get hit by lightning etc.)
# r_A = rewards (in C) of tree to fungus. cost for tree, gain for fungus. 
# e_A = efficiency of conversion for tree C to fungus C. 
# m_A = density independent mortality of fungus
# d_A = density independent mortality of fungus
# u_P = nutrient uptake rate of tree 
# u_A = nutrient uptake rate of fungus
# l = leakiness of tree. nutrients that leave the tree and reenter the soil pool (decomposed leaves, root leakiness, etc.)
# 
# Parameter Set

g = 5 #growth rate of tree is higher than fungus
s = 1 #senescence is low compared to growth rate
e_A = .1 #standard energy flow between two trophic levels
r_A = .8 #close to optimal investment (see IV.)
m_A = 0 #set to zero so ZNGI goes through origin (see III.)
u_A = 1 #uptake of fungus is greater than plant to provide benefit ... ensure mutualism
u_P = 0.01 #uptake of plant lower than fungus 
l = 1 #small compared to total pool
d_A = 1 #easy value
d_P = 0 #set to zero to prevent equations from being insane
N_T = 5 #initial soil nutrient pool level

tset <- seq(from = 0, to = 1000, length.out = 5000) #create a time vector going from 0 to 1000. these numbers are fairly arbitrary, but good to show multiple "seasons" of length 200 later.
P.simu <- NaN*tset; P.simu[1] <- .1
A.simu <- NaN*tset; A.simu[1] <- .1 #both tree and plant start small, but not at zero. model if for growth, not generation.
Np.simu <- NaN*tset; Np.simu[1] <- .1 #relatively arbitrary start number. 
Ns.simu <- NaN*tset; Ns.simu[1] <- N_T - Np.simu[1] #Nutrients in plant (Np) and nutrients in soil (Ns) add to total nutrients in system (Nt)

for(i in 2:length(tset)){ #for loop with the equations for N, P, Np, and Ns within the model. 
  dt <- tset[i] - tset[i-1]
  P <- P.simu[i-1]
  A <- A.simu[i-1]
  N_P <- Np.simu[i-1]
  N_S <- Ns.simu[i-1]
  dP <- (P* (g*N_P/P - r_A*A - s*(1+d_P*P)))*dt
  dA <- (A* (e_A*r_A*P-m_A-d_A*A))*dt
  dNp <- (u_P*P*N_S + u_A*A*N_S - l*N_P)*dt
  dNs <- (-u_P*P*N_S - u_A*A*N_S + l*N_P)*dt
  P.simu[i] <- P + dP
  A.simu[i] <- A + dA
  Np.simu[i] <- N_P + dNp
  Ns.simu[i] <- N_S + dNs
  
}

Pcol <- 'forestgreen' #color of P in graphs
Acol <- 'darksalmon' #color of A in graphs

simulation_results = data.frame(time = tset, Plant = P.simu, Fungus = A.simu)

ggplot(simulation_results, aes(x = time, y = Plant)) +
  geom_line()

# plot(x = tset, y = P.simu, type = 'l', col=Pcol, las = 1, lwd=2, xlab = 'Time', ylab = 'Biomass in C', ylim = c(0, max(c(P.simu,A.simu))))
# abline(h = 0) # I like to have this line as a reference
# lines(x = tset, y = A.simu, col = Acol, lwd = 2)
# legend(x = max(tset)*0.5, y = .35*max(c(P.simu,A.simu)), legend=c('P','A'),lwd= 2, col=c(Pcol,Acol), horiz=TRUE)
# 
# tail(P.simu)
# tail(A.simu) #tail functions allow you to see what the plant / fungus grows to in the end.