library(tidyr)
library(dplyr)
library(ggplot2)
library(DEoptim)

log_likelihood <- function(params, data) {
  params <- as.numeric(params)
  names(params) <- c("p1", "p2", "p4", "c1", "c2", "c4", "m1", "m2", "m4")
  
  x <- as.numeric(data$X)
  z <- as.numeric(data$Z)
  t <- as.numeric(data$t)
  
  p1 <- params["p1"]; p2 <- params["p2"]; p4 <- params["p4"]
  c1 <- params["c1"]; c2 <- params["c2"]; c4 <- params["c4"]
  m1 <- params["m1"]; m2 <- params["m2"]; m4 <- params["m4"]
  
  A <- p1 / (1 + c1 * (t - m1)^2)
  B <- p2 / (1 + c2 * (t - m2)^2)
  D <- p4 / (1 + c4 * (t - m4)^2)
  # C <- 1-A-B
  
  epsilon <- 1e-10
  # A <- pmax(A, epsilon)
  # B <- pmax(B, epsilon)
  C <- pmax(1 - A - B-D, epsilon)
  
  ll <- sum((x == 1 & z == 0)*log(A + epsilon) + 
              (x==0 & z == 0)*log(B + epsilon) + 
              (x==0 & z == 1)*log(D + epsilon) +
              (x==-1 & z == 0)*log(C))
  
  return(as.numeric(-ll)) 
}

r_estimator <- function(data){
  Si <- data$Si
  deltat <- data$deltat
  n <- nrow(data)
  r <- n/(sum(Si*deltat))
  return(r)
}

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


p1.true <- 0.25
c1.true <- 0.005
m1.true <- 6


p2.true <- 0.65
c2.true <- 0.012
m2.true <- 15

p4.true <- 0.2
c4.true <- 0.008
m4.true <- 20

s0.true <- 200
r.true <- 0.2

reps <- 100 

sim.markov <- simulate_markov(
  reps = reps,
  p1 = p1.true, c1 = c1.true, m1 = m1.true,
  p2 = p2.true, c2 = c2.true, m2 = m2.true,
  p4 = p4.true, c4 = c4.true, m4 = m4.true,  
  s0 = s0.true,
  r = r.true,
  seed = 922125
)

res.prm <- data.frame()
start_time <- Sys.time()
for (j in 1:reps) {
  sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == j]),
                         Z = diff(sim.markov$dc.steps[sim.markov$reps == j]),
                         t = sim.markov$time.steps[sim.markov$reps == j][-1],
                         deltat = diff(sim.markov$time.steps[sim.markov$reps == j]),
                         Si = sim.markov$sc.steps[sim.markov$reps == j]
                         [-length(sim.markov$sc.steps[sim.markov$reps == j])])
  t_min <- min(sim.data$t)
  t_max <- max(sim.data$t, 100)
  lower <- c(p1 = 1e-5, p2 = 1e-5, p4 = 1e-5,
             c1 = 1e-5, c2 = 1e-5, c4 = 1e-5,
             m1 = t_min, m2 = t_min, m4 = t_min)
  upper <- c(p1 = 0.99999, p2 = 0.99999, p4 = 0.99999,
             c1 = 100, c2 = 100, c4 = 100,
             m1 = t_max, m2 = t_max, m4 = t_max)
  set.seed(1234 + j)
  DE_res <- DEoptim(
    fn = log_likelihood,
    lower = lower,
    upper = upper,
    data = sim.data,
    DEoptim.control(NP = 150, itermax = 600, F = 0.9, CR = 0.8, strategy = 6,
                    trace = FALSE)
  )
  par.est <- DE_res$optim$bestmem
  r.est <- r_estimator(sim.data)
  res.prm <- rbind(res.prm, c(par.est,r.est))
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
res.prm %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(mean = mean(Value),
            median = median(Value),
            percentile050 = quantile(Value, probs = 0.05),
            percentile950 = quantile(Value, probs = 0.95))

# write.csv(res.prm, "full_likelihood_start200_055_005_04_015_012_12_020_008_20_02.csv")

# write.csv(res.prm, "full_likelihood_start200_030_008_06_030_008_15_030_008_25_02.csv")

# write.csv(res.prm, "full_likelihood_start200_025_005_06_065_012_15_010_008_20_02.csv")

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