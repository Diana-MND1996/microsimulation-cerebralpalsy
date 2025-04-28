############################################################################################
################# Microsimulation modeling using R: a tutorial #### 2025 ###################
############################################################################################
# This code forms the basis for the microsimulation model of the article: 
#
# Diana Marcela Nova Díaz, Sergio Aguilera Albesa, Eduardo Sánchez Iriso. 
# Cost-Effectiveness of Complementary and Alternative Therapies for Children with Cerebral Palsy: 
# Evidence from Real-World Data and Microsimulation Modelling

# Please cite the article when using this code

############################################################################################
################# Code of Appendix A.1 #######################################################
############################################################################################
#rm(list = ls())  # remove any variables in R's memory 

##################################### Model input #########################################
# Model input
n.i   <- 100000                     # Number of individuals, here is the number of patients in the study
n.t   <- 30                       # 30-year time horizon
v.n   <- c("GMFCS I-II", "GMFCS III", "GMFCS IV-V", "Dead")  # Health statuses: GMFCS I-II (Near-healthy: H), GMFCS III (Sick: S1), GMFCS IV-V (Sicker: S2), Dead (D)
n.s   <- length(v.n)               # The number of health states
v.M_1 <- rep("GMFCS I-II", n.i)       # All start in the healthy state (H: GMFCS I-II)
d.c   <- d.e <- 0.03               # Descuento de costos y QALYs al 3%
v.Trt <- c("Standard treatment", "Hippo", "Homo", "Peto", "Ther")

# Transition probabilities (per cycle ajusted to GMFCS level)
p.HD    <- 0.002                   # probability to die when GMFCS I-II
p.HS1   <- 0.08                    # probability to become GMFCS III when GMFCS I-II
p.S1H   <- 0.015           	       # probability to become GMFCS I-II when GMFCS III
p.S1S2  <- 0.05         	         # probability to become GMFCS IV-V when GMFCS III
rr.S1   <- 3             	         # rate ratio of death when GMFCS III vs GMFCS I-II
rr.S2   <- 10            	         # rate ratio of death when GMFCS IV-V vs GMFCS I-II 
r.HD    <- -log(1 - p.HD) 	       # rate of death when GMFCS I-II 
r.S1D   <- rr.S1 * r.HD  	         # rate of death when GMFCS III
r.S2D   <- rr.S2 * r.HD  	         # rate of death when GMFCS IV-V
p.S1D   <- 1 - exp(- r.S1D)        # probability to die when GMFCS III
p.S2D   <- 1 - exp(- r.S2D)        # probability to die when GMFCS IV-V

# Cost and utility inputs (per cycle adjusted to GMFCS level)
c.H     <- 4992.11                 # cost of remaining one cycle GMFCS I-II
c.S1    <- 5611.4                  # cost of remaining one cycle GMFCS III
c.S2    <- 5756.06                 # cost of remaining one cycle GMFCS IV-V
c.HiI   <- 5121.66                 # cost of Hippotherapy treatment  (per cycle)
c.HiIII   <- 10040                 # cost of Hippo treatment (per cycle)
c.HiV <- 12073.1                   # cost of Hippo treatment (per cycle)
c.HoI   <- 8310.85                 # cost of Homo treatment (per cycle)
c.HoIII    <- 9770                 # cost of Homo  therapies (per cycle)
c.HoV  <- 11228.38                 # cost of Homo therapies (per cycle)
c.PI   <- 8350.66                  # cost of Peto treatment (per cycle)
c.PIII    <- 7920.75               # cost of Peto therapies (per cycle)
c.PV  <- 10843.16                  # cost of Peto therapies (per cycle)
c.TI   <- 7887                     # cost of Ther treatment (per cycle)
c.TIII    <- 9348                  # cost of Ther therapies (per cycle)
c.TV  <- 11573                     # cost of Ther therapies (per cycle)



u.H     <- 0.7042                   # utility when GMFCS I-II 
u.S1    <- 0.5197                   # utility when GMFCS III
u.S2    <- 0.0527                    # utility when GMFCS IV-V
u.HiI   <- 0.3250                    # utility when GMFCS I and being treated with Hippo
u.HiIII   <- 0.8783                    # utility when GMFCS III and being treated with Hippo
u.HiV <- 0.2970                    # utility when GMFCS IV-V and being treated with Hippo
u.HoI   <- 0.8792                    # utility when GMFCS I and being treated with Homo
u.HoIII    <- 0.7751                    # utility when GMFCS III and being treated with Homo
u.HoV   <- 0.3246                    # utility when GMFCS IV-V and being treated with Homo
u.PI    <- 0.8702                    # utility when GMFCS I and being treated with Peto 
u.PIII   <- 0.7963                    # utility when GMFCS III and being treated with Peto 
u.PV   <- 0.3390                    # utility when GMFCS IV-V and being treated with Peto
u.TI    <- 1                    # utility when GMFCS I and being treated with Ther
u.TIII   <- 0.7549                    # utility when GMFCS III and being treated with Ther 
u.TV   <- 0.3993                    # utility when GMFCS IV-V and being treated with Ther 



##################################### Functions ###########################################

# The MicroSim function for the simple microsimulation of the 'GMFCS III-GMFCS IV-V' model keeps track of what happens to each individual during each cycle. 


MicroSim <- function(v.M_1, n.i, n.t, v.n, d.c, d.e, TR.out = TRUE, TS.out = TRUE, Trt = 1, seed = 1) {
  # Arguments:  
  # v.M_1:   vector of initial states for individuals 
  # n.i:     number of individuals
  # n.t:     total number of cycles to run the model
  # v.n:     vector of health state names
  # d.c:     discount rate for costs
  # d.e:     discount rate for health outcome (QALYs)
  # TR.out:  should the output include a microsimulation trace? (default is TRUE)
  # TS.out:  should the output include a matrix of transitions between states? (default is TRUE)
  # Trt:     are the n.i individuals receiving treatment? (scalar with a Boolean value, default is FALSE)
  # seed:    starting seed number for random number generator (default is 1)
  # Makes use of:
  # Probs:   function for the estimation of transition probabilities
  # Costs:   function for the estimation of cost state values
  # Effs:    function for the estimation of state specific health outcomes (QALYs)
  
  v.dwc <- 1 / (1 + d.c) ^ (0:n.t)   # # calculate the cost discount weight based on the discount rate d.c
  v.dwe <- 1 / (1 + d.e) ^ (0:n.t)   # calculate the QALY discount weight based on the discount rate d.e
  
  # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <-  matrix(nrow = n.i, ncol = n.t + 1, 
                               dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                               paste("cycle", 0:n.t, sep = " ")))  
  
  m.M[, 1] <- v.M_1                     # indicate the initial health state  
  
  for (i in 1:n.i) {
    set.seed(seed + i)                  # set the seed for every individual for the random number generator
    m.C[i, 1] <- Costs(m.M[i, 1], Trt)  # estimate costs per individual for the initial health state conditional on treatment
    m.E[i, 1] <- Effs(m.M[i, 1], Trt)   # estimate QALYs per individual for the initial health state conditional on treatment
    
    for (t in 1:n.t) {
      v.p <- Probs(m.M[i, t])            # calculate the transition probabilities at cycle t 
      
      m.M[i, t + 1] <- sample(v.n, prob = v.p, size = 1)  # sample the next health state and store that state in matrix m.M 
      m.C[i, t + 1] <- Costs(m.M[i, t + 1], Trt)   # estimate costs per individual during cycle t + 1 conditional on treatment
      m.E[i, t + 1] <- Effs(m.M[i, t + 1], Trt)   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    }
    if(i/100 == round(i/100,0)) {      
      cat('\r', paste(i/n.i * 100, "% done", sep = " "))
    }
  } # close the loop for the individuals 
  
  tc <- m.C %*% v.dwc        # total (discounted) cost per individual
  te <- m.E %*% v.dwe        # total (discounted) QALYs per individual 
  
  tc_hat <- mean(tc)         # average (discounted) cost
  te_hat <- mean(te)         # average (discounted) QALYs
  
  if (TS.out == TRUE) {  # create a  matrix of transitions across states
    TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
    TS <- matrix(TS, nrow = n.i)
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ") # name the rows 
    colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")  # name the columns   
  } else {
    TS <- NULL
  }
  
  if (TR.out == TRUE) {  # create a trace from the individual trajectories
    TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
    TR <- TR / n.i                                      # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")    # name the rows
    colnames(TR) <- v.n                                 # name the columns 
  } else {
    TR <- NULL
  }
  
  results <- list(m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR) # store the results from the simulation in a list
  return(results) # return the results
} # end of the MicroSim function 

#### Probability function
# The Probs function that updates the transition probabilities of every cycle is shown below

Probs <- function(M_it) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  v.p.it <- rep(NA, n.s)     # create vector of state transition probabilities
  names(v.p.it) <- v.n       # name the vector
  
  v.p.it[M_it == "GMFCS I-II"]  <- c(1 - p.HS1 - p.HD, p.HS1, 0, p.HD)                   # transition probabilities when GMFCS I-II               
  v.p.it[M_it == "GMFCS III"]     <- c(p.S1H, 1- p.S1H - p.S1S2 - p.S1D, p.S1S2, p.S1D)  # transition probabilities when GMFCS III
  v.p.it[M_it == "GMFCS IV-V"]   <- c(0, 0, 1 - p.S2D, p.S2D)                            # transition probabilities when GMFCS IV-V
  v.p.it[M_it == "Dead"]     <- c(0, 0, 0, 1)                                            # transition probabilities when dead 
  
  ifelse(sum(v.p.it) == 1, return(v.p.it), print("Probabilities do not sum to 1"))       # return the transition probabilities or produce an error
}       

### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function(M_it, Trt = 1) {
  c.it <- 0
  if (M_it == "GMFCS I-II") {
    c.it <- c.H
    if (Trt == 1) c.it <- c.HiI
    if (Trt == 2) c.it <- c.HoI
    if (Trt == 3) c.it <- c.PI
    if (Trt == 4) c.it <- c.TI

  } else if (M_it == "GMFCS III") {
    c.it <- c.S1
    if (Trt == 1) c.it <- c.HiIII
    if (Trt == 2) c.it <- c.HoIII
    if (Trt == 3) c.it <- c.PIII
    if (Trt == 4) c.it <- c.TIII
  } else if (M_it == "GMFCS IV-V") {
    c.it <- c.S2
    if (Trt == 1) c.it <- c.HiV
    if (Trt == 2) c.it <- c.HoV
    if (Trt == 3) c.it <- c.PV
    if (Trt == 4) c.it <- c.TV
  }
  return(c.it)
}


### Health outcome function 
# The Effs function to update the utilities at every cycle.

Effs <- function(M_it, Trt = 1) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual treated? (default is FALSE) 
  # cl:   cycle length (default is 1)
  
  u.it <- 0
  if (M_it == "GMFCS I-II") {
    u.it <- u.H
   if (Trt == 1) u.it <- u.HiI
    if (Trt == 2) u.it <- u.HoI
    if (Trt == 3) u.it <- u.PI
    if (Trt == 4) u.it <- u.TI

  } else if (M_it == "GMFCS III") {
   if (Trt == 1) u.it <- u.HiIII
    if (Trt == 2) u.it <- u.HoIII
    if (Trt == 3) u.it <- u.PIII
    if (Trt == 4) u.it <- u.TIII
    if (Trt == 0) u.it <- u.S1
  } else if (M_it == "GMFCS IV-V") {
   if (Trt == 1) u.it <- u.HiV
    if (Trt == 2) u.it <- u.HoV
    if (Trt == 3) u.it <- u.PV
    if (Trt == 4) u.it <- u.TV
    if (Trt == 0) u.it <- u.S2
  }
  return(u.it)
}


##################################### Run the simulation ##################################

sim_S_trt  <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 0)   # run for Standard treatment
sim_Hi     <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 1)   # run for Hi
sim_Ho     <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 2)   # run for Ho therapies
sim_P    <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 3)   # run for P therapies
sim_T    <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 4)   # run for T therapies



cat("Average cost with standard treatment:", sim_S_trt$tc_hat, "\n")
cat("Average cost with Hi therapy:", sim_Hi$tc_hat, "\n")
cat("Average cost with Ho therapy:", sim_Ho$tc_hat, "\n")
cat("Average cost with P therapy:", sim_P$tc_hat, "\n")
cat("Average cost with T therapy:", sim_T$tc_hat, "\n")
cat("Average QALYs with standard treatment:", sim_S_trt$te_hat, "\n")
cat("Average QALYs with Hi therapy:", sim_Hi$te_hat, "\n")
cat("Average QALYs with Ho therapy:", sim_Ho$te_hat, "\n")
cat("Average QALYs with P  therapy:", sim_P$te_hat, "\n")
cat("Average QALYs with T  therapy:", sim_T$te_hat, "\n")

################################# Cost-effectiveness analysis #############################

# store the mean costs (and the MCSE) of each strategy in a new variable v.C (vector costs)
v.C  <- c(sim_S_trt$tc_hat, sim_Hi$tc_hat, sim_Ho$tc_hat, sim_P$tc_hat, sim_T$tc_hat)  # The costs for the four strategies
se.C <- c(sd(sim_S_trt$tc), sd(sim_Hi$tc), sd(sim_Ho$tc), sd(sim_P$tc), sd(sim_T$tc)) / sqrt(n.i)  # MCSE of costs

# store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)
v.E  <- c(sim_S_trt$te_hat, sim_Hi$te_hat, sim_Ho$te_hat, sim_P$te_hat, sim_T$te_hat)  # The QALYs for the four strategies
se.E <- c(sd(sim_S_trt$te), sd(sim_Hi$te), sd(sim_Ho$te), sd(sim_P$te), sd(sim_T$te)) / sqrt(n.i)  # # MCSE of QALYs

# Calcular los costos y QALYs incrementales
delta.C <- v.C[2] - v.C[1]    # calculate incremental costs
delta.E <- v.E[2] - v.E[1]    # calculate incremental QALYs
se.delta.C <- sd(sim_Hi$tc - sim_S_trt$tc) / sqrt(n.i)  # # Monte Carlo squared error (MCSE) of incremental costs
se.delta.E <- sd(sim_Hi$te - sim_S_trt$te) / sqrt(n.i)  # Monte Carlo squared error (MCSE) of incremental QALYs
ICER <- delta.C / delta.E    # calculate the ICER

# Calculate earned QALYs
QALYs_gained_Hi_vs_S_trt <- sim_Hi$te_hat - sim_S_trt$te_hat   # QALYs gained between hi and st Treatment
QALYs_gained_Ho_vs_S_trt <- sim_Ho$te_hat - sim_S_trt$te_hat      # QALYs gained between Ho and St Treatment
QALYs_gained_P_vs_S_trt <- sim_P$te_hat - sim_S_trt$te_hat    # QALYs gained between P Therapy and St Treatment
QALYs_gained_T_vs_S_trt <- sim_T$te_hat - sim_S_trt$te_hat      # QALYs gained between T Therapy and St Treatment

# Show earned QALYs
cat("QALYs gained with Hippo vs. St treatment:", QALYs_gained_Hi_vs_S_trt, "\n")
cat("QALYs gained with Homo therapy vs. St treatment:", QALYs_gained_Ho_vs_S_trt, "\n")
cat("QALYs gained with Peto therapy vs. St treatment:", QALYs_gained_P_vs_S_trt, "\n")
cat("QALYs gained with Thera therapy vs. St treatment:", QALYs_gained_T_vs_S_trt, "\n")


# Print the table
table_micro <- data.frame(
  Treatment = c("Standard treatment", "Hippo", "Homo", "Peto", "Ther"),
  Costs = c(round(v.C[1], 0), round(v.C[2], 0), round(v.C[3], 0), round(v.C[4], 0), round(v.C[5], 0)),
  "MCSE Costs" = c(round(se.C[1], 0), round(se.C[2], 0), round(se.C[3], 0), round(se.C[4], 0), round(se.C[5], 0)),
  QALYs = c(round(v.E[1], 3), round(v.E[2], 3), round(v.E[3], 3), round(v.E[4], 3), round(v.E[5], 3)),
  "MCSE QALYs" = c(round(se.E[1], 3), round(se.E[2], 3), round(se.E[3], 3), round(se.E[4], 3), round(se.E[5], 3)),
  "Incremental Costs" = c("", round(delta.C, 0), round(v.C[3] - v.C[1], 0), round(v.C[4] - v.C[1], 0), round(v.C[5] - v.C[1], 0)),
  "Incremental QALYs" = c("", round(delta.E, 3), round(v.E[3] - v.E[1], 3), round(v.E[4] - v.E[1], 3), round(v.E[5] - v.E[1], 3)),
  "QALYs Gained (vs St Treatment)" = c("", QALYs_gained_Hi_vs_S_trt, QALYs_gained_Ho_vs_S_trt, QALYs_gained_P_vs_S_trt, QALYs_gained_T_vs_S_trt),
  ICER = c("", round(ICER, 0), round((v.C[3] - v.C[1]) / (v.E[3] - v.E[1]), 0), round((v.C[4] - v.C[1]) / (v.E[4] - v.E[1]), 0), round((v.C[5] - v.C[1]) / (v.E[5] - v.E[1]), 0))
)

# Showing results
print(table_micro)


