library(DEoptim)
library(parallelly)
library(dplyr)
library(tidyr)


# True parameters
p1.true <- 0.55
c1.true <- 0.005
m1.true <- 4


p2.true <- 0.15
c2.true <- 0.012
m2.true <- 12

p4.true <- 0.2
c4.true <- 0.008
m4.true <- 20

s0.true <- 200
r.true <- 0.2

# Simulation set-up
min_p_t <- 0.0001
max_it <- 50000
seed0 <- 922125
t_max <- 100
reps <- 100

f_t1 <- function(t, c1, m1) { 1 / (1 + c1 * (t - m1)^2)}

f_t2 <- function(t, c2, m2) { 1 / (1 + c2 * (t - m2)^2)}

f_t4 <- function(t, c4, m4) { 1 / (1 + c4 * (t - m4)^2)}

probs_t <- function(t, p1, p2, p4, c1, c2, c4, m1, m2, m4) {
  f1t <- f_t1(t, c1, m1)
  f2t <- f_t2(t, c2, m2)
  f4t <- f_t4(t, c4, m4)
  
  p1t <- p1 * f1t
  p2t <- p2 * f2t
  p4t <- p4 * f4t
  p3t <- 1 - p1t - p2t - p4t
  
  if (any(p3t < 0)) warning("Some p3(t) values are negative — check parameters!")
  cbind(p1t, p2t, p3t, p4t)
}

# simulate data
set.seed(seed0)

sim.markov <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric())

for (i in 1:reps) {
  sc.current <- s0.true
  ec.current <- dc.current <- time.current <- 0
  sc.steps <- c(s0.true)
  ec.steps <- c(ec.current)
  dc.steps <- c(dc.current)
  time.steps <- c(time.current)
  itnum <- 0
  while (sc.current > 0) {
    itnum <- itnum + 1
    time.current <- time.current + rexp(1, rate = r.true*sc.current)
    probs <- probs_t(t = time.current, p1 = p1.true, p2 = p2.true,
                     p4 = p4.true, c1 = c1.true, c2 = c2.true, c4 = c4.true, m1 = m1.true, m2 = m2.true, m4 = m4.true)
    
    if ((sum(probs < min_p_t) > 0 | (itnum > max_it))) {
      cat(q_t, "\n")
      break}
    if (sum(probs) > 1) {
      probs <- probs/sum(probs)
    }
    event_type <- sample(1:4, 1, prob=probs)
    switch(event_type,
           "1" = {
             sc.current <- sc.current + 1
             ec.current <- ec.current
             dc.current <- dc.current
           },
           "2" = {
             sc.current <- sc.current
             ec.current <- ec.current + 1
             dc.current <- dc.current
           },
           "3" = {
             sc.current <- sc.current - 1
             ec.current <- ec.current + 2
             dc.current <- dc.current
           },
           "4" = {
             sc.current <- sc.current
             ec.current <- ec.current
             dc.current <- dc.current + 1
           })
    sc.steps <- c(sc.steps, sc.current)
    ec.steps <- c(ec.steps, ec.current)
    dc.steps <- c(dc.steps, dc.current)
    time.steps <- c(time.steps, time.current)
  }
  sim.markov <- rbind(sim.markov,cbind(reps = rep(i, length(sc.steps)), time.steps, sc.steps, ec.steps, dc.steps))
}
head(sim.markov,10)


loglik <- function(params, data) {
  
  # safety wrapper to prevent NaN while sampling
  safe_fail <- function() return(1e12)
  
  params <- as.numeric(params)
  if (length(params) != 10) return(safe_fail())
  
  names(params) <- c("p1", "p2", "p4", "c1", "c2", "c4", "m1", "m2", "m4", "r")
  
  # extract data
  M <- as.numeric(data$M)
  Y <- as.numeric(data$Y)
  t <- as.numeric(data$t)

  n_steps <- length(t)

  # parameters
  p1 <- params["p1"]; p2 <- params["p2"]; p4 <- params["p4"]
  c1 <- params["c1"]; c2 <- params["c2"]; c4 <- params["c4"]
  m1 <- params["m1"]; m2 <- params["m2"]; m4 <- params["m4"]
  r  <- params["r"]
  
  # Hidden state set
  hidden_states <- 0
  alpha_prev <- 1
  
  scale_factor <- rep(1, n_steps)
  
  for (step in 2:n_steps) {
    deltaT <- t[step] - t[step-1]
    deltaY <- Y[step] - Y[step-1]
    deltaM <- M[step] - M[step-1]
    
   
    # event probabilities (clamped to [1e-12,1-1e-12]) ---
    p1_t <- p1 / (1 + c1*(t[step] - m1)^2)
    p2_t <- p2 / (1 + c2*(t[step] - m2)^2)
    p4_t <- p4 / (1 + c4*(t[step] - m4)^2)
    
    if (!is.finite(p1_t) || !is.finite(p2_t) || !is.finite(p4_t)) return(safe_fail())
    
    clamp <- function(x) max(min(x, 1 - 1e-12), 1e-12)
    p1_t <- clamp(p1_t)
    p2_t <- clamp(p2_t)
    p4_t <- clamp(p4_t)
    
    p_other <- 1 - p1_t - p2_t - p4_t
    p_other <- clamp(p_other)
    
    rates <- r * (M[step-1] - hidden_states)
    
    rates[!is.finite(rates)] <- 1e-8
    rates[rates < 1e-8] <- 1e-8
    
    # exponential PDF
    exp_term <- dexp(deltaT, rates)
    if (any(!is.finite(exp_term))) return(safe_fail())
    
    # p2
    if (deltaY == 1 && deltaM == 0) {
      alpha_new <- alpha_prev * exp_term * p2_t
    }
    # p3
    else if (deltaY == 2 && deltaM == -1) {
      alpha_new <- alpha_prev * exp_term * p_other
      
      # Remove illegal hidden states
      valid <- (hidden_states <= M[step])
      if (!any(valid)) return(safe_fail())
      
      alpha_new <- alpha_new[valid]
      hidden_states <- hidden_states[valid]
    }
    # p1 or p4
    else if (deltaY == 0 && deltaM == 1) {
      
      diag_part  <- alpha_prev * exp_term * p1_t
      shift_part <- alpha_prev * exp_term * p4_t
      
      if (any(!is.finite(diag_part)) || any(!is.finite(shift_part)))
        return(safe_fail())
      
      # If diag_part length=0 -> no valid state
      L <- length(diag_part)
      if (L == 0) return(safe_fail())
      
      # Create combined vector
      alpha_new <- numeric(L + 1)
      alpha_new[1:L]     <- diag_part
      alpha_new[2:(L+1)] <- alpha_new[2:(L+1)] + shift_part
      
      # Add new hidden state (increment)
      hidden_states <- c(hidden_states, max(hidden_states)+1)
      
      # Remove illegal states
      valid <- (hidden_states <= M[step])
      if (!any(valid)) return(safe_fail())
      
      alpha_new     <- alpha_new[valid]
      hidden_states <- hidden_states[valid]
    }
    
    else {
      return(safe_fail())
    }
    
    s <- sum(alpha_new)
    if (!is.finite(s) || s <= 0) return(safe_fail())
    
    scale_factor[step] <- s / scale_factor[step-1]
    alpha_prev <- alpha_new / s
  }
  
  total <- sum(log(scale_factor))
  if (!is.finite(total)) return(safe_fail())
  
  return(-total)
}

res.prm <- data.frame()
start_time <- Sys.time()
for (j in 1:reps) {
  sim.data <- data.frame(t = sim.markov$time.steps[sim.markov$reps == j],
                         Y = sim.markov$ec.steps[sim.markov$reps == j],
                         M = sim.markov$sc.steps[sim.markov$reps == j] +
                           sim.markov$dc.steps[sim.markov$reps == j])
  t_min <- min(sim.data$t)
  t_max <- max(sim.data$t)
  lower <- c(p1 = 1e-5, p2 = 1e-5, p4 = 1e-5,
             c1 = 1e-5, c2 = 1e-5, c4 = 1e-5,
             m1 = t_min, m2 = t_min, m4 = t_min,
             r = 1e-5)
  upper <- c(p1 = 0.99999, p2 = 0.99999, p4 = 0.99999,
             c1 = 100, c2 = 100, c4 = 100,
             m1 = t_max, m2 = t_max, m4 = t_max,
             r = 0.99999)
  set.seed(1234 + j)
  DE_res <- DEoptim(
    fn = loglik,
    lower = lower,
    upper = upper,
    data = sim.data,
    # DEoptim.control(NP = 150, itermax = 500, F = 0.9, CR = 0.8, strategy = 6,
    #                 trace = FALSE, paralleType = 1)
    DEoptim.control(NP = 150, itermax = 600, F = 0.9, CR = 0.8, strategy = 6,
                    trace = FALSE, parallelType = 1)
  )
  par.est <- DE_res$optim$bestmem
  res.prm <- rbind(res.prm, par.est)
  print(j)
}
end_time <- Sys.time()
end_time-start_time

colnames(res.prm) <- c("p1", "p2", "p4", "c1", "c2", "c4",
                       "m1", "m2", "m4", "r")

res.prm %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(mean = mean(Value),
            median = median(Value),
            percentile025 = quantile(Value, probs = 0.025),
            percentile975 = quantile(Value, probs = 0.975))

write.csv(res.prm, "start200_055_008_04_015_012_12_020_008_20_02.csv", row.names = FALSE)


###
### original function
forward_log_likelihood <- function(params, data){
  params <- as.numeric(params)
  names(params) <- c("p1", "p2", "p4", "c1", "c2", "c4", "m1", "m2", "m4", "r")
  
  M <- as.numeric(data$M)
  Y <- as.numeric(data$Y)
  t <- as.numeric(data$t)
  
  p1 <- params["p1"]; p2 <- params["p2"]; p4 <- params["p4"]
  c1 <- params["c1"]; c2 <- params["c2"]; c4 <- params["c4"]
  m1 <- params["m1"]; m2 <- params["m2"]; m4 <- params["m4"]
  r <- params["r"]
  
  hidden_states <- c(0) 
  alpha_prev <- matrix(1, ncol = length(hidden_states)) 
  
  n_steps <- length(t)
  scale_factor <- numeric(n_steps)
  scale_factor[1] <- 1 
  
  for (step in 2:n_steps) {
    
    deltaT <- t[step] - t[step-1]
    deltaY <- Y[step] - Y[step-1]
    deltaM <- M[step] - M[step-1]
    
    # Event probabilities
    p1_t <- p1 / (1 + c1*(t[step] - m1)^2)
    p2_t <- p2 / (1 + c2*(t[step] - m2)^2)
    p4_t <- p4 / (1 + c4*(t[step] - m4)^2)
    
    if (deltaY == 1 && deltaM == 0) {
      trans_mat <- diag(length(hidden_states))
      
      for (i in seq_along(hidden_states)) {
        rate <- max(r * (M[step-1] - hidden_states[i]), 1e-8)
        trans_mat[i,i] <- dexp(deltaT, rate) * p2_t
      }
      
      alpha_new <- alpha_prev %*% trans_mat
    }
    
    else if (deltaY == 2 && deltaM == -1) {
      trans_mat <- diag(length(hidden_states))
      p3_t <- max(1 - p1_t - p2_t - p4_t, 1e-12)
      
      for (i in seq_along(hidden_states)) {
        rate <- max(r * (M[step-1] - hidden_states[i]), 1e-8)
        trans_mat[i,i] <- dexp(deltaT, rate) * p3_t
      }
      
      alpha_new <- alpha_prev %*% trans_mat
      
      # Remove hidden states > M(step)
      hidden_states <- hidden_states[hidden_states <= M[step]]
      alpha_new <- alpha_new[1:length(hidden_states)]
    }
    
    else if (deltaY == 0 && deltaM == 1) {
      
      # Build matrix with one potential new state
      trans_mat <- matrix(0, nrow = length(hidden_states), ncol = length(hidden_states) + 1)
      
      for (i in seq_along(hidden_states)) {
        rate <- max(r * (M[step-1] - hidden_states[i]), 1e-8)
        trans_mat[i,i]   <- dexp(deltaT, rate) * p1_t
        trans_mat[i,i+1] <- dexp(deltaT, rate) * p4_t
      }
      
      alpha_new <- alpha_prev %*% trans_mat
      
      # Add new hidden state but only if <= M
      new_state <- max(hidden_states) + 1
      hidden_states <- c(hidden_states, new_state)
      hidden_states <- hidden_states[hidden_states <= M[step]]
      
      alpha_new <- alpha_new[1:length(hidden_states)]
    }
    alpha_sum <- sum(alpha_new)
    if (alpha_sum == 0 || is.na(alpha_sum)) alpha_sum <- 1e-12
    
    scale_factor[step] <- alpha_sum/scale_factor[step-1]
    alpha_prev <- alpha_new / alpha_sum
  }
  
  log_likelihood <- sum(log(scale_factor))
  return(-log_likelihood)
}

loglik(params = c(p1.true, p2.true, p4.true, c1.true, c2.true, c4.true,
                  m1.true, m2.true, m4.true, r.true), data = sim.data)
forward_log_likelihood(params = c(p1.true, p2.true, p4.true, c1.true, c2.true, c4.true,
                                  m1.true, m2.true, m4.true, r.true), data = sim.data)

