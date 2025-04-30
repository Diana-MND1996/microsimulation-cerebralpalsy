############################################################################################
################# Microsimulation modeling using R: a tutorial #### 2025 ###################
############################################################################################
# This code forms the basis for the microsimulation model of the article: 
#
# Diana Marcela Nova Díaz, Sergio Aguilera Albesa, Eduardo Sánchez-Iriso. 
# Cost-Effectiveness Analysis of Complementary and Alternative Therapies for Children with Cerebral Palsy: 
# Evidence from Real-World Data and Microsimulation Modelling

# Please cite the article when using this code

############################################################################################
################# Code of Appendix A.2 #######################################################
############################################################################################
#rm(list = ls())  # remove any variables in R's memory 

##################################### Model input #########################################
# Model input
n.i   <- 100000                    # Number of individuals, here is the number of patients in the study
n.t   <- 30                        # 30-year time horizon
v.n   <- c("GMFCS I-II", "GMFCS III", "GMFCS IV-V", "Dead")  # Health statuses: GMFCS I-II (Near-healthy: H), GMFCS III (Sick: S1), GMFCS IV-V (Sicker: S2), Dead (D)
n.s   <- length(v.n)               # The number of health states
v.M_1 <- rep("GMFCS I-II", n.i)    # All start in the healthy state (H: GMFCS I-II)
d.c   <- d.e <- 0.03               # Discount rate for costs and QALYs 3%
v.Trt <- c("Standard treatment", "SpeechTH", "Rehab","Occupth")

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
c.H     <- 4838.75                 # cost of remaining one cycle GMFCS I-II
c.S1    <- 5693.22                 # cost of remaining one cycle GMFCS III
c.S2    <- 5757.18                 # cost of remaining one cycle GMFCS IV-V
c.SpI   <- 7640                    # cost of Speech therapy (per cycle)  being GMFCS I-II
c.SpIII   <- 12913.33              # cost of Speech therapy (per cycle)  being GMFCS III
c.SpV <- 13579.90                  # cost of Speech therapy (per cycle)  being GMFCS IV-V
c.RI   <- 10664.58                 # cost of Rehabilitation (per cycle)  being GMFCS I-II
c.RIII    <- 12153.64              # cost of Rehabilitation (per cycle)  being GMFCS III
c.RV  <- 14715.89                  # cost Rehabilitation (per cycle)     being GMFCS IV-V
c.OI   <- 11847.5                  # cost of Occupational therapy (per cycle) being GMFCS I-II
c.OIII    <- 12533.48              # cost of Occupational therapy (per cycle) mean value of SpIII and RIII (missing) being GMFCS III
c.OV  <- 14985.23                  # cost of Occupational therapy (per cycle) being GMFCS IV-V



u.H     <- 0.6879                  # utility when GMFCS I-II 
u.S1    <- 0.4943                  # utility when GMFCS III
u.S2    <- 0                       # utility when GMFCS IV-V
u.SpI   <- 0.8901                  # utility when GMFCS I-II and being treated with Speech therapy
u.SpIII   <- 0.7907                # utility when GMFCS III and being treated with Speech therapy
u.SpV <- 0.3034                    # utility when GMFCS IV-V and being treated with Speech therapy
u.RI   <- 0.8932                   # utility when GMFCS I-II and being treated with Rehabilitation
u.RIII    <- 0.7791                # utility when GMFCS III and being treated with Rehabilitation
u.RV   <- 0.34                     # utility when GMFCS IV-V and being treated with Rehabilitation
u.OI    <- 0.8584                  # utility when GMFCS I-II and being treated with Occupational therapy 
u.OIII   <- 0.7849                 # utility when GMFCS III and being treated with Occupational therapy mean value of SpIII and RIII (missing) 
u.OV   <- 0.3593                   # utility when GMFCS IV-V and being treated with Occupational therapy 

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
      m.C[i, t + 1] <- Costs(m.M[i, t + 1], Trt)          # estimate costs per individual during cycle t + 1 conditional on treatment
      m.E[i, t + 1] <- Effs(m.M[i, t + 1], Trt)           # estimate QALYs per individual during cycle t + 1 conditional on treatment
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
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
    colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")   # name the columns   
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
    if (Trt == 1) c.it <- c.SpI
    if (Trt == 2) c.it <- c.RI
    if (Trt == 3) c.it <- c.OI
  } else if (M_it == "GMFCS III") {
    c.it <- c.S1
    if (Trt == 1) c.it <- c.SpIII
    if (Trt == 2) c.it <- c.RIII
    if (Trt == 3) c.it <- c.OIII
  } else if (M_it == "GMFCS IV-V") {
    c.it <- c.S2
    if (Trt == 1) c.it <- c.SpV
    if (Trt == 2) c.it <- c.RV
    if (Trt == 3) c.it <- c.OV
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
    if (Trt == 1) u.it <- u.SpI
    if (Trt == 2) u.it <- u.RI
    if (Trt == 3) u.it <- u.OI
  } else if (M_it == "GMFCS III") {
    if (Trt == 1) u.it <- u.SpIII
    if (Trt == 2) u.it <- u.RIII
    if (Trt == 3) u.it <- u.OIII
    if (Trt == 0) u.it <- u.S1
  } else if (M_it == "GMFCS IV-V") {
    if (Trt == 1) u.it <- u.SpV
    if (Trt == 2) u.it <- u.RV
    if (Trt == 3) u.it <- u.OV
    if (Trt == 0) u.it <- u.S2
  }
  return(u.it)
}


##################################### Run the simulation ##################################

sim_S_trt  <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 0)   # run for Standard treatment
sim_Sp     <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 1)   # run for Speech therapy
sim_Reh     <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 2)   # run for Rehabilitation
sim_Ocup    <- MicroSim(v.M_1, n.i, n.t, v.n, d.c, d.e, Trt = 3)   # run for Occupational


cat("Average cost with standard treatment:", sim_S_trt$tc_hat, "\n")
cat("Average cost with Speech therapy:", sim_Sp$tc_hat, "\n")
cat("Average cost with Rehabilitation therapy:", sim_Reh$tc_hat, "\n")
cat("Average cost with Occupational therapy:", sim_Ocup$tc_hat, "\n")
cat("Average QALYs with standard treatment:", sim_S_trt$te_hat, "\n")
cat("Average QALYs with Speech therapy:", sim_Sp$te_hat, "\n")
cat("Average QALYs with Rehabilitation therapy:", sim_Reh$te_hat, "\n")
cat("Average QALYs with Occupational  therapy:", sim_Ocup$te_hat, "\n")

################################# Cost-effectiveness analysis #############################

# store the mean costs (and the MCSE) of each strategy in a new variable v.C (vector costs)
v.C  <- c(sim_S_trt$tc_hat, sim_Sp$tc_hat, sim_Reh$tc_hat, sim_Ocup$tc_hat)  # The costs for the four strategies
se.C <- c(sd(sim_S_trt$tc), sd(sim_Sp$tc), sd(sim_Reh$tc), sd(sim_Ocup$tc)) / sqrt(n.i)  # MCSE of costs

# store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)
v.E  <- c(sim_S_trt$te_hat, sim_Sp$te_hat, sim_Reh$te_hat, sim_Ocup$te_hat)  # The QALYs for the four strategies
se.E <- c(sd(sim_S_trt$te), sd(sim_Sp$te), sd(sim_Reh$te), sd(sim_Ocup$te)) / sqrt(n.i)  # # MCSE of QALYs

# Calcular los costos y QALYs incrementales
delta.C <- v.C[2] - v.C[1]    # calculate incremental costs
delta.E <- v.E[2] - v.E[1]    # calculate incremental QALYs
se.delta.C <- sd(sim_Sp$tc - sim_S_trt$tc) / sqrt(n.i)  # # Monte Carlo squared error (MCSE) of incremental costs
se.delta.E <- sd(sim_Sp$te - sim_S_trt$te) / sqrt(n.i)  # Monte Carlo squared error (MCSE) of incremental QALYs
ICER <- delta.C / delta.E    # calculate the ICER

# Calculate earned QALYs
QALYs_gained_Sp_vs_S_trt <- sim_Sp$te_hat - sim_S_trt$te_hat   # QALYs gained between Speech therapy and Standard Treatment
QALYs_gained_Reh_vs_S_trt <- sim_Reh$te_hat - sim_S_trt$te_hat      # QALYs gained between Rehabilitation and Standard Treatment
QALYs_gained_Ocup_vs_S_trt <- sim_Ocup$te_hat - sim_S_trt$te_hat    # QALYs gained between Occupational therapy and Standard Treatment

# Show earned QALYs
cat("QALYs gained with Speech therapy vs. St treatment:", QALYs_gained_Sp_vs_S_trt, "\n")
cat("QALYs gained with Rehabilitation vs. St treatment:", QALYs_gained_Reh_vs_S_trt, "\n")
cat("QALYs gained with Occupational therapy vs. St treatment:", QALYs_gained_Ocup_vs_S_trt, "\n")

# Print the table
table_micro <- data.frame(
  Treatment = c("Standard treatment", "SpeechTH", "Rehab","Occupth"),
  Costs = c(round(v.C[1], 0), round(v.C[2], 0), round(v.C[3], 0), round(v.C[4], 0)),
  "MCSE Costs" = c(round(se.C[1], 0), round(se.C[2], 0), round(se.C[3], 0), round(se.C[4], 0)),
  QALYs = c(round(v.E[1], 3), round(v.E[2], 3), round(v.E[3], 3), round(v.E[4], 3)),
  "MCSE QALYs" = c(round(se.E[1], 3), round(se.E[2], 3), round(se.E[3], 3), round(se.E[4], 3)),
  "Incremental Costs" = c("", round(delta.C, 0), round(v.C[3] - v.C[1], 0), round(v.C[4] - v.C[1], 0)),
  "Incremental QALYs" = c("", round(delta.E, 3), round(v.E[3] - v.E[1], 3), round(v.E[4] - v.E[1], 3)),
  "QALYs Gained (vs St Treatment)" = c("", QALYs_gained_Sp_vs_S_trt, QALYs_gained_Reh_vs_S_trt, QALYs_gained_Ocup_vs_S_trt),
  ICER = c("", round(ICER, 0), round((v.C[3] - v.C[1]) / (v.E[3] - v.E[1]), 0), round((v.C[4] - v.C[1]) / (v.E[4] - v.E[1]), 0))
)

# Showing results
print(table_micro)


