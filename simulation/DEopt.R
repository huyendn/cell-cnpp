library(DEoptim)
library(ggplot2)
library(dplyr)
library(tidyr)

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
  
  return(as.numeric(-ll)) 
}

r_estimator <- function(data){
  Si <- data$Si
  deltat <- data$deltat
  n <- nrow(data)
  r <- n/(sum(Si*deltat))
  return(r)
}


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
    fn = log_likelihood,
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


