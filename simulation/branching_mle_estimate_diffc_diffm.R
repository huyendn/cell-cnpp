library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)

# True parameters
p1.true <- 0.55
c1.true <- 0.005
m1.true <- 4


p2.true <- 0.15
c2.true <- 0.01
m2.true <- 18

parm.true <- c(p1.true, c1.true, m1.true, p2.true, c2.true, m2.true)
s0.true <- 200
r.true <- 0.2

# Simulation set-up
min_p_t <- 0.0001
max_it <- 50000
seed0 <- 922125
t_max <- 100
reps <- 100

###
### plot probability functions
p1_fun <- function(t) p1.true / (1 + c1.true * (t - m1.true)^2)
p2_fun <- function(t) p2.true / (1 + c2.true * (t - m2.true)^2)
p3_fun <- function(t) 1 - p1_fun(t) - p2_fun(t)

t_vals <- seq(0, 60, length.out = 500)
data <- data.frame(
  t = t_vals,
  p1 = p1_fun(t_vals),
  p2 = p2_fun(t_vals),
  p3 = p3_fun(t_vals)
)

data_long <- pivot_longer(data, cols = c(p1, p2, p3), names_to = "Function", values_to = "Value")

ggplot(data_long, aes(x = t, y = Value, color = Function)) +
  geom_line(linewidth = 1.2) +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = c("skyblue3", "indianred3", "darkgreen"),
                     labels = c(expression(p[1](t)), expression(p[2](t)), expression(p[3](t)))) +
  labs(
    title = "Functions p1(t), p2(t), and p3(t)",
    x = "t",
    y = "Probability",
    color = "Function"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.text.align = 0
  )

### 
### simulate data
f_t1 <- function(t, c1, m1) {
  1 / (1 + c1 * (t - m1)^2)
}

f_t2 <- function(t, c2, m2) {
  1 / (1 + c2 * (t - m2)^2)
}

probs_t <- function(t, p1, p2, c1, c2, m1, m2) {
  f1t <- f_t1(t, c1, m1)
  f2t <- f_t2(t, c2, m2)
  
  p1t <- p1 * f1t
  p2t <- p2 * f2t
  p3t <- 1 - p1t - p2t
  
  if (any(p3t < 0)) warning("Some p3(t) values are negative — check parameters!")
  
  cbind(p1t, p2t, p3t)
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
                     c1 = c1.true, c2 = c2.true,
                     m1 = m1.true, m2 = m2.true)
    
    if ((sum(probs < min_p_t) > 0 | (itnum > max_it))) {
      cat(q_t, "\n")
      break
    }
    if (sum(probs) > 1) {
      probs <- probs/sum(probs)
    }
    event_type <- sample(1:3, 1, prob=probs)
    switch(event_type,
           "1" = {
             sc.current <- sc.current + 1
             ec.current <- ec.current
           },
           "2" = {
             sc.current <- sc.current
             ec.current <- ec.current + 1
           },
           "3" = {
             sc.current <- sc.current - 1
             ec.current <- ec.current + 2
           })
    sc.steps <- c(sc.steps, sc.current)
    ec.steps <- c(ec.steps, ec.current)
    time.steps <- c(time.steps, time.current)
  }
  sim.markov <- rbind(sim.markov,cbind(reps = rep(i, length(sc.steps)), time.steps, sc.steps, ec.steps))
}

###
### plot one replication

index.test <- sample(1:100, size = 1) # test the optimization on one replications
index.test
sim.markov %>% filter(reps == index.test) %>%
  ggplot() +
  geom_point(aes(x = time.steps, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
  geom_point(aes(x = time.steps, y = sc.steps, col = "stem cell"), alpha = 0.5)+
  labs(x = "time (in hours)", y = "cell counts", color = "",
       title = "Simulated cell counts")+
  scale_color_manual(values = c("differentiated cell" = "skyblue3",
                                "stem cell" = "indianred1",
                                "expected differentiated cell" = "navy",
                                "expected stem cell" = "brown4")) +
  theme_classic() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13))

### plot stopping time
sim.markov %>% filter(sc.steps == 0) %>%
  ggplot() +
  geom_histogram(aes(x = time.steps), 
                 binwidth = 1.5, 
                 color = "black",fill = "skyblue3") +
  labs(x = "time", 
       title = "Distribution of stopping time",
       subtitle = paste0("Number of replications: ",reps)) +
  theme_minimal()

sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == index.test]),
                       t = sim.markov$time.steps[sim.markov$reps == index.test][-1],
                       deltat = diff(sim.markov$time.steps[sim.markov$reps == index.test]),
                       Si = sim.markov$sc.steps[sim.markov$reps == index.test][-length(sim.markov$sc.steps[sim.markov$reps == index.test])])
head(sim.data)


###
###
###
log_likelihood <- function(params, data) {
  x <- data$X
  t <- data$t
  
  p1 <- params["p1"]
  p2 <- params["p2"]
  c1 <- params["c1"]
  c2 <- params["c2"]
  m1 <- params["m1"]
  m2 <- params["m2"]
  
  A <- p1 / (1 + c1 * (t - m1)^2)
  B <- p2 / (1 + c2 * (t - m2)^2)
  C <- 1 - A - B
  
  epsilon <- 1e-10 # small positive number to avoid log(0)
  ll <- sum((x == 1) * log(pmax(A, epsilon)) + (x == 0) * log(pmax(B, epsilon)) +
              (x == -1) * log(pmax(C, epsilon))
  )
  return(ll)
}


log_likelihood_gradient <- function(params, data) {
  x <- data$X
  t <- data$t
  
  p1 <- params["p1"]
  p2 <- params["p2"]
  c1 <- params["c1"]
  c2 <- params["c2"]
  m1 <- params["m1"]
  m2 <- params["m2"]
  
  A <- p1 / (1 + c1 * (t - m1)^2)
  B <- p2 / (1 + c2 * (t - m2)^2)
  C <- 1 - A - B
  
  dA_dp1 <- 1 / (1 + c1 * (t - m1)^2)
  dA_dc1 <- -p1 * (t - m1)^2 / (1 + c1 * (t - m1)^2)^2
  dA_dm1 <- 2 * p1 * c1 * (t - m1) / (1 + c1 * (t - m1)^2)^2
  
  dB_dp2 <- 1 / (1 + c2 * (t - m2)^2)
  dB_dc2 <- -p2 * (t - m2)^2 / (1 + c2 * (t - m2)^2)^2
  dB_dm2 <- 2 * p2 * c2 * (t - m2) / (1 + c2 * (t - m2)^2)^2
  
  dC_dp1 <- -dA_dp1
  dC_dc1 <- -dA_dc1
  dC_dm1 <- -dA_dm1
  dC_dp2 <- -dB_dp2
  dC_dc2 <- -dB_dc2
  dC_dm2 <- -dB_dm2
  
  grad <- numeric(6)
  names(grad) <- c("p1", "p2", "c1", "c2", "m1", "m2")
  
  grad["p1"] <- sum((x == 1) * (1 / A) * dA_dp1 + (x == -1) * (1 / C) * dC_dp1)
  grad["p2"] <- sum((x == 0) * (1 / B) * dB_dp2 + (x == -1) * (1 / C) * dC_dp2)
  grad["c1"] <- sum((x == 1) * (1 / A) * dA_dc1 + (x == -1) * (1 / C) * dC_dc1)
  grad["c2"] <- sum((x == 0) * (1 / B) * dB_dc2 + (x == -1) * (1 / C) * dC_dc2)
  grad["m1"] <- sum((x == 1) * (1 / A) * dA_dm1 + (x == -1) * (1 / C) * dC_dm1)
  grad["m2"] <- sum((x == 0) * (1 / B) * dB_dm2 + (x == -1) * (1 / C) * dC_dm2)
  
  return(grad)
}


run_optim <- function(start_params, data) {
  optim(
    par = start_params,
    fn = log_likelihood,
    gr = log_likelihood_gradient,
    data = data,
    method = "L-BFGS-B",
    lower = c(p1 = 1e-5, p2 = 1e-5, c1 = 1e-5, c2 = 1e-5, m1 = 1e-5, m2 = 1e-5),
    upper = c(p1 = 1, p2 = 1, c1 = Inf, c2 = Inf, m1 = -Inf, m2 = -Inf),
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


res.prm <- data.frame()

for (j in 1:reps) {
  sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == j]),
                         t = sim.markov$time.steps[sim.markov$reps == j][-1],
                         deltat = diff(sim.markov$time.steps[sim.markov$reps == j]),
                         Si = sim.markov$sc.steps[sim.markov$reps == j]
                         [-length(sim.markov$sc.steps[sim.markov$reps == j])])
  results <- lapply(starts, run_optim, data = sim.data)
  best <- results[[which.max(sapply(results, function(r) r$value))]]
  r.est <- r_estimator(sim.data)
  res.prm <- rbind(res.prm, c(best$par,r.est))
}

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
