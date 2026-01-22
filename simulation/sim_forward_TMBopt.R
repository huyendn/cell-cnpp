library(TMB)
library(tidyr)
library(dplyr)
library(ggplot2)
library(DEoptim)
library(parallelly)

compile("loglik.cpp")       # compile C++ code for automatic differentiation
dyn.load(dynlib("loglik"))  # load compiled code for automatic differentiation


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

reps <- 3
simulate_markov <- function(
    reps = 100,
    p1, p2, p4,
    c1, c2, c4,
    m1, m2, m4,
    s0 = 200,
    r = 0.2,
    min_p_t = 1e-4,
    max_it = 50000,
    t_max = 100,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  f_t <- function(t, c, m) 1 / (1 + c * (t - m)^2)
  
  probs_t <- function(t) {
    f1 <- f_t(t, c1, m1)
    f2 <- f_t(t, c2, m2)
    f4 <- f_t(t, c4, m4)
    
    p1t <- p1 * f1
    p2t <- p2 * f2
    p4t <- p4 * f4
    p3t <- 1 - p1t - p2t - p4t
    
    p3t <- pmax(p3t, 0)
    
    probs <- c(p1t, p2t, p3t, p4t)
    
    s <- sum(probs)
    if (s <= 0) stop("All transition probabilities zero — invalid parameters.")
    if (s != 1) probs <- probs / s
    
    # avoid extremely small values
    probs[probs < min_p_t] <- 0
    probs <- probs / sum(probs)
    
    probs
  }
  
  all_out <- list()
  idx <- 1
  
  for (rep_i in 1:reps) {
    
    sc <- s0
    ec <- 0
    dc <- 0
    t <- 0
    itnum <- 0
    
    sc_steps <- sc
    ec_steps <- ec
    dc_steps <- dc
    time_steps <- t
    
    while (sc > 0) {
      
      itnum <- itnum + 1
      if (itnum > max_it || t > t_max) {
        warning("Stopped early: max_it or time bound exceeded")
        break
      }
      
      t <- t + rexp(1, rate = r * sc)
      
      probs <- probs_t(t)
      
      event <- sample(1:4, size = 1, prob = probs)
      
      if (event == 1) { sc <- sc + 1 }
      if (event == 2) { ec <- ec + 1 }
      if (event == 3) { sc <- sc - 1; ec <- ec + 2 }
      if (event == 4) { dc <- dc + 1 }
      
      sc_steps <- c(sc_steps, sc)
      ec_steps <- c(ec_steps, ec)
      dc_steps <- c(dc_steps, dc)
      time_steps <- c(time_steps, t)
    }
    
    all_out[[idx]] <- data.frame(
      reps = rep(rep_i, length(sc_steps)),
      time.steps = time_steps,
      sc.steps = sc_steps,
      ec.steps = ec_steps,
      dc.steps = dc_steps
    )
    idx <- idx + 1
  }
  
  do.call(rbind, all_out)
}

sim.markov <- simulate_markov(
  reps = reps,
  p1 = p1.true, c1 = c1.true, m1 = m1.true,
  p2 = p2.true, c2 = c2.true, m2 = m2.true,
  p4 = p4.true, c4 = c4.true, m4 = m4.true,  
  s0 = s0.true,
  r = r.true,
  seed = 922125
)

head(sim.markov)


loglik_R <- function(params, data) {
  
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
    
    # scale_factor[step] <- s / scale_factor[step-1]
    scale_factor[step] <- s 
    alpha_prev <- alpha_new / s
  }
  
  total <- sum(log(scale_factor))
  
  if (!is.finite(total)) return(safe_fail())
  
  return(-total)
}

## de optimization

res.prm.de <- data.frame()

start_time <- Sys.time()
for (j in 1:reps) {
  sim.data <- data.frame(t = sim.markov$time.steps[sim.markov$reps == j],
                         Y = sim.markov$ec.steps[sim.markov$reps == j],
                         M = sim.markov$sc.steps[sim.markov$reps == j] +
                           sim.markov$dc.steps[sim.markov$reps == j])

  lower <- c(p1 = 1e-5, p2 = 1e-5, p4 = 1e-5,
             c1 = 1e-5, c2 = 1e-5, c4 = 1e-5,
             m1 = 1e-5, m2 = 1e-5, m4 = 1e-5,
             r = 1e-5)
  upper <- c(p1 = 0.99999, p2 = 0.99999, p4 = 0.99999,
             c1 = 100, c2 = 100, c4 = 100,
             m1 = 100, m2 = 100, m4 = 100,
             r = 0.99999)
  set.seed(1234 + j)
  DE_res <- DEoptim(
    fn = loglik_R,
    lower = lower,
    upper = upper,
    data = sim.data,
    # DEoptim.control(NP = 150, itermax = 500, F = 0.9, CR = 0.8, strategy = 6,
    #                 trace = FALSE, paralleType = 1)
    DEoptim.control(NP = 150, itermax = 600, F = 0.9, CR = 0.8, strategy = 6,
                    trace = FALSE, parallelType = 1)
  )
  par.est.de <- DE_res$optim$bestmem
  res.prm <- rbind(res.prm.de, par.est.de)
  print(j)
}
end_time <- Sys.time()
end_time-start_time

colnames(res.prm.de) <- c("p1", "p2", "p4", "c1", "c2", "c4",
                       "m1", "m2", "m4", "r")

res.prm.de %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(mean = mean(Value),
            median = median(Value),
            percentile025 = quantile(Value, probs = 0.025),
            percentile975 = quantile(Value, probs = 0.975))

write.csv(res.prm.de, "destep1_start200_055_005_04_015_012_12_020_008_20_02.csv", row.names = FALSE)
# write.csv(res.prm.de, "destep1_start200_030_008_06_030_008_15_030_008_25_02.csv", row.names = FALSE)
# write.csv(res.prm.de, "destep1_start200_025_005_06_065_012_15_010_008_20_02.csv", row.names = FALSE)


## result from differential evolution
res.de <- read.csv("de_res/destep1_start200_055_005_04_015_012_12_020_008_20_02.csv")
# res.de <- read.csv("de_res/destep1_start200_025_005_06_065_012_15_010_008_20_02.csv")
# res.de <- read.csv("de_res/destep1_start200_030_008_06_030_008_15_030_008_25_02.csv")

res.prm <- data.frame()
start_time <- Sys.time()
for (j in 1:reps) {
  sim.data <- data.frame(t = sim.markov$time.steps[sim.markov$reps == j],
                         Y = sim.markov$ec.steps[sim.markov$reps == j],
                         M = sim.markov$sc.steps[sim.markov$reps == j] +
                           sim.markov$dc.steps[sim.markov$reps == j])
  data_list <- list(
    t = sim.data$t,
    M = sim.data$M,
    Y = sim.data$Y
  )
  de.est <- res.de[j,]
  parameters <- list(
    p1_raw = qlogis(as.numeric(de.est["p1"])),
    p2_raw = qlogis(as.numeric(de.est["p2"])),
    p4_raw = qlogis(as.numeric(de.est["p4"])),
    r_raw  = qlogis(as.numeric(de.est["r"])),
    logit_c1 = qlogis(as.numeric(de.est["c1"]/100)),
    logit_c2 = qlogis(as.numeric(de.est["c2"]/100)),
    logit_c4 = qlogis(as.numeric(de.est["c4"]/100)),
    logit_m1 = qlogis(as.numeric(de.est["m1"]/100)),
    logit_m2 = qlogis(as.numeric(de.est["m2"]/100)),
    logit_m4 = qlogis(as.numeric(de.est["m4"]/100))
  )
  obj <- MakeADFun(data_list, parameters, DLL="loglik", silent=TRUE)
  
  fit <- optim(obj$par, obj$fn, obj$gr,
               method = "BFGS")
  est <- fit$par
  p1 <- plogis(est["p1_raw"])
  p2 <- plogis(est["p2_raw"])
  p4 <- plogis(est["p4_raw"])
  r  <- plogis(est["r_raw"])
  c1 <- 100*plogis(est["logit_c1"])
  c2 <- 100*plogis(est["logit_c2"])
  c4 <- 100*plogis(est["logit_c4"])
  m1 <- 100*plogis(est["logit_m1"])
  m2 <- 100*plogis(est["logit_m2"])
  m4 <- 100*plogis(est["logit_m4"])
  
  par.est <- c(p1, p2, p4, c1, c2, c4, m1, m2, m4, r)
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
c(p1.true, p2.true, p4.true, c1.true, c2.true, c4.true, m1.true, m2.true, m4.true, r.true, s0.true)


res.prm %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = "", y = Value)) +
  geom_violin(aes(fill = Variable), color = "black") +
  geom_boxplot(fill = NA, width = 0.5)+
  facet_wrap(~ Variable, scale = "free_y", ncol = 3) +
  geom_jitter(width = 0.1, alpha = 0.3) + 
  theme_minimal() +
  labs(title = "Violin Plot for Each Parameter Estimates", x = "",
       fill= "Parameters", y = "Estimates")
c(p1.true, p2.true, p4.true, c1.true, c2.true, c4.true, m1.true, m2.true, m4.true, r.true, s0.true)

