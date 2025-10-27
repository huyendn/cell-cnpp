library(DEoptim)
log_likelihood <- function(params, data) {
  params <- as.numeric(params)
  names(params) <- c("p1", "p2", "c1", "c2", "m1", "m2")
  
  x <- as.numeric(data$X)
  t <- as.numeric(data$t)
  
  p1 <- params["p1"]; p2 <- params["p2"]
  c1 <- params["c1"]; c2 <- params["c2"]
  m1 <- params["m1"]; m2 <- params["m2"]
  
  A <- p1 / (1 + c1 * (t - m1)^2)
  B <- p2 / (1 + c2 * (t - m2)^2)
  # C <- 1-A-B
  
  epsilon <- 1e-10
  # A <- pmax(A, epsilon)
  # B <- pmax(B, epsilon)
  C <- pmax(1 - A - B, epsilon)
  
  ll <- sum((x==1)*log(A + epsilon) + 
              (x==0)*log(B + epsilon) + 
              (x==-1)*log(C))
  
  return(as.numeric(ll)) 
}

log_likelihood_gradient <- function(params, data) {
  params <- as.numeric(params)
  names(params) <- c("p1", "p2", "c1", "c2", "m1", "m2")
  
  x <- as.numeric(data$X)
  t <- as.numeric(data$t)
  
  p1 <- params["p1"]; p2 <- params["p2"]
  c1 <- params["c1"]; c2 <- params["c2"]
  m1 <- params["m1"]; m2 <- params["m2"]
  
  A <- p1 / (1 + c1 * (t - m1)^2)
  B <- p2 / (1 + c2 * (t - m2)^2)
  epsilon <- 1e-10
  C <- pmax(1 - A - B, epsilon)
  
  dA_dp1 <- 1 / (1 + c1 * (t - m1)^2)
  dA_dc1 <- -p1 * (t - m1)^2 / (1 + c1 * (t - m1)^2)^2
  dA_dm1 <-  2 * p1 * c1 * (t - m1) / (1 + c1 * (t - m1)^2)^2
  
  dB_dp2 <- 1 / (1 + c2 * (t - m2)^2)
  dB_dc2 <- -p2 * (t - m2)^2 / (1 + c2 * (t - m2)^2)^2
  dB_dm2 <-  2 * p2 * c2 * (t - m2) / (1 + c2 * (t - m2)^2)^2
  
  dC_dp1 <- -dA_dp1
  dC_dp2 <- -dB_dp2
  dC_dc1 <- -dA_dc1
  dC_dc2 <- -dB_dc2
  dC_dm1 <- -dA_dm1
  dC_dm2 <- -dB_dm2
  
  grad <- numeric(6)
  names(grad) <- c("p1", "p2", "c1", "c2", "m1", "m2")
  
  grad["p1"] <- sum((x == 1)  * (dA_dp1 / A) + (x == -1) * (dC_dp1 / C))
  grad["p2"] <- sum((x == 0)  * (dB_dp2 / B) + (x == -1) * (dC_dp2 / C))
  grad["c1"] <- sum((x == 1)  * (dA_dc1 / A) + (x == -1) * (dC_dc1 / C))
  grad["c2"] <- sum((x == 0)  * (dB_dc2 / B) + (x == -1) * (dC_dc2 / C))
  grad["m1"] <- sum((x == 1)  * (dA_dm1 / A) + (x == -1) * (dC_dm1 / C))
  grad["m2"] <- sum((x == 0)  * (dB_dm2 / B) + (x == -1) * (dC_dm2 / C))
  
  return(as.numeric(grad))
}


run_optim <- function(start_params, data) {
  optim(
    par = start_params,
    fn = log_likelihood,
    gr = log_likelihood_gradient,
    data = data,
    method = "L-BFGS-B",
    lower = c(p1 = 1e-5, p2 = 1e-5, c1 = 1e-5, c2 = 1e-5, m1 = 1e-5, m2 = 1e-5),
    upper = c(p1 = 0.99999, p2 = 0.99999, c1 = Inf, c2 = Inf, m1 = Inf, m2 = Inf),
    control = list(fnscale = -1)
  )
}

r_estimator <- function(data){
  Si <- data$Si
  deltat <- data$deltat
  n <- nrow(data)
  r <- n/(sum(Si*deltat))
  return(r)
}

log_likelihood_deopt <- function(params, data) {
  p1 <- params[1]; p2 <- params[2]
  c1 <- params[3]; c2 <- params[4]
  m1 <- params[5]; m2 <- params[6]
  
  x <- as.numeric(data$X)
  t <- as.numeric(data$t)
  
  A <- pmax(p1 / (1 + c1 * (t - m1)^2), 1e-12)
  B <- pmax(p2 / (1 + c2 * (t - m2)^2), 1e-12)
  C <- pmax(1 - A - B, 1e-12)
  
  ll <- sum((x == 1) * log(A) +
              (x == 0) * log(B) +
              (x == -1) * log(C))
  return(-ll)  # DEoptim minimizes
}

# t_min <- min(sim.data$t)
# t_max <- max(sim.data$t)
# 
# lower <- c(p1 = 1e-5, p2 = 1e-5, c1 = 1e-5, c2 = 1e-5, m1 = t_min, m2 = t_min)
# upper <- c(p1 = 0.99999, p2 = 0.99999, c1 = 100, c2 = 100, m1 = t_max, m2 = t_max)
# 
# 
# DE_res <- DEoptim(
#   fn = log_likelihood_deopt,
#   lower = lower,
#   upper = upper,
#   data = sim.data,
#   DEoptim.control(NP = 60, itermax = 200, trace = FALSE)
# )
# 
# DE_res$optim$bestmem

res.prm <- data.frame()
start_time <- Sys.time()
for (j in 1:reps) {
  sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == j]),
                         t = sim.markov$time.steps[sim.markov$reps == j][-1],
                         deltat = diff(sim.markov$time.steps[sim.markov$reps == j]),
                         Si = sim.markov$sc.steps[sim.markov$reps == j]
                         [-length(sim.markov$sc.steps[sim.markov$reps == j])])
  t_min <- min(sim.data$t)
  t_max <- max(sim.data$t)
  lower <- c(p1 = 1e-5, p2 = 1e-5, c1 = 1e-5, c2 = 1e-5, m1 = t_min, m2 = t_min)
  upper <- c(p1 = 0.99999, p2 = 0.99999, c1 = 100, c2 = 100, m1 = t_max, m2 = t_max)
  set.seed(1234 + j)
  DE_res <- DEoptim(
    fn = log_likelihood_deopt,
    lower = lower,
    upper = upper,
    data = sim.data,
    DEoptim.control(NP = 60, itermax = 200, trace = FALSE)
  )
  par.est <- DE_res$optim$bestmem
  r.est <- r_estimator(sim.data)
  res.prm <- rbind(res.prm, c(par.est,r.est))
}
end_time <- Sys.time()
end_time-start_time

colnames(res.prm) <- c("p1", "p2", "c1", "c2", "m1", "m2", "r")

res.prm %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(mean = mean(Value),
            median = median(Value),
            percentile025 = quantile(Value, probs = 0.025),
            percentile975 = quantile(Value, probs = 0.975))
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


