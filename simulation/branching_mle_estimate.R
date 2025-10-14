library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)

# True parameters
p1.true <- 0.5
p2.true <- 0.2
c.true <- 0.005
m.true <- 4

s0.true <- 200
r.true <- 0.2

# Simulation set-up
min_p_t <- 0.0001
max_it <- 50000
seed0 <- 922125
t_max <- 100
reps <- 100

# Probability mass function
f_t <- function(t, c, m) {
  1 / (1 + c * (t - m)^2)
}

probs_t <- function(t, p1, p2, c, m) {
  ft <- f_t(t, c, m)
  p1t <- p1 * ft
  p2t <- p2 * ft
  p3t <- 1 - p1t - p2t
  cbind(p1t, p2t, p3t)
}
# probs_t(0.05,p1 = p1.true, p2 = p2.true, c = c.true, m = m.true)
# p1.true/(1+c.true*(0.05-m.true)^2)

# Simulate data
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
    probs <- probs_t(t = time.current, p1 = p1.true, p2 = p2.true, c = c.true, m = m.true)
    
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

sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == index.test]),
                       t = sim.markov$time.steps[sim.markov$reps == index.test][-1],
                       deltat = diff(sim.markov$time.steps[sim.markov$reps == index.test]),
                       Si = sim.markov$sc.steps[sim.markov$reps == index.test][-length(sim.markov$sc.steps[sim.markov$reps == index.test])])


### Optimization
neg_log_lik <- function(par, data) { #negative log-likelihood
  X <- data$X
  t <- data$t
  deltat <- data$deltat
  Si <- data$Si
  
  p1 <- par[1]
  p2 <- par[2]
  c  <- par[3]
  m  <- par[4]
  r  <- par[5]  
  
  if (p1 <= 0 || p2 <= 0 || c <= 0 || (p1 + p2) >= 1 || r <= 0) return(Inf)
  
  ft <- 1 / (1 + c * (t - m)^2)
  
  loglik <- numeric(length(X))
  loglik[X == 1]  <- log(p1 * ft[X == 1])
  loglik[X == 0]  <- log(p2 * ft[X == 0])
  loglik[X == -1] <- log(1 - (p1 + p2) * ft[X == -1])
  
  # interarrival time
  loglik_exp <- log(r * Si) - r * Si * deltat
  
  return(-sum(loglik) - sum(loglik_exp))
}

grad_log_lik <- function(par, data) { # gradient
  X <- data$X
  t <- data$t
  deltat <- data$deltat
  Si <- data$Si
  
  p1 <- par[1]
  p2 <- par[2]
  c  <- par[3]
  m  <- par[4]
  r  <- par[5]
  
  if (p1 <= 0 || p2 <= 0 || c <= 0 || (p1 + p2) >= 1 || r <= 0) {
    return(rep(NA, 5))
  }
  
  ft <- 1 / (1 + c * (t - m)^2)
  one_minus_p_ft <- 1 - (p1 + p2) * ft
  
  I1 <- as.numeric(X == 1)
  I0 <- as.numeric(X == 0)
  Im <- as.numeric(X == -1)
  
  denom_m <- one_minus_p_ft
  denom_m[denom_m == 0] <- 1e-10  # avoid division by zero
  
  # Derivatives of f(t)
  df_dc <- - (t - m)^2 * ft^2
  df_dm <- 2 * c * (t - m) * ft^2
  
  # pmf part
  dL_dp1 <- sum(I1 / p1 - Im * ft / denom_m)
  dL_dp2 <- sum(I0 / p2 - Im * ft / denom_m)
  
  temp_c <- ( (I1 + I0) / ft - Im * (p1 + p2) / denom_m )
  dL_dc <- sum(temp_c * df_dc)
  dL_dm <- sum(temp_c * df_dm)
  
  # exponential part
  dL_dr <- sum(1 / r - Si * deltat)
  
  -c(dL_dp1, dL_dp2, dL_dc, dL_dm, dL_dr)
}

ui <- matrix(c(
  1, 0, 0, 0, 0,   # p1 > 0
  0, 1, 0, 0, 0,   # p2 > 0
  0, 0, 1, 0, 0,   # c > 0
  -1, -1, 0, 0, 0, # p1 + p2 < 1
  0, 0, 0, 0, 1    # r > 0
), byrow = TRUE, nrow = 5)

ci <- c(0, 0, 0, -1, 0)

start <- c(p1 = 0.5, p2 = 0.4, c = 10, m = 10, r = 0.01)


result <- constrOptim(
  theta = start,
  f = neg_log_lik,
  grad = grad_log_lik,
  ui = ui,
  ci = ci,
  data = sim.data
)

result$par


###
###
### estimate on all replications
start.time <- Sys.time()
# start <- c(0.1, 0.1, 5, 5, 1)
start <- c(p1 = 0.2, p2 = 0.4, c = 10, m = 10, r = 1)

res.prm <- data.frame(p1 = numeric(), p2 = numeric(), c = numeric(),
                      m = numeric(), r = numeric())

for (j in 1:reps) {
  sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == j]),
                         t = sim.markov$time.steps[sim.markov$reps == j][-1],
                         deltat = diff(sim.markov$time.steps[sim.markov$reps == j]),
                         Si = sim.markov$sc.steps[sim.markov$reps == j]
                         [-length(sim.markov$sc.steps[sim.markov$reps == j])])
  result <- constrOptim(
    theta = start,
    f = neg_log_lik,
    grad = grad_log_lik,
    ui = ui,
    ci = ci,
    data = sim.data
  )
  res.prm <- rbind(res.prm, result$par)
}
colnames(res.prm) <- c("p1", "p2", "c", "m", "r")
dim(res.prm)

end.time <- Sys.time()

end.time-start.time

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
  facet_wrap(~ Variable, scale = "free_y", ncol = 5) +
  geom_jitter(width = 0.1, alpha = 0.3) + 
  theme_minimal() +
  labs(title = "Violin Plot for Each Parameter Estimates", x = "",
       fill= "Parameters", y = "Estimates")

# plot probablity mass function
p1_fun <- function(t) p1.true / (1 + c.true * (t - m.true)^2)
p2_fun <- function(t) p2.true / (1 + c.true * (t - m.true)^2)
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


####
####
####
#### including p4

# True parameters
p1.true <- 0.5
p2.true <- 0.2
p4.true <- 0.05
c.true <- 0.005
m.true <- 4

s0.true <- 200
r.true <- 0.2

# Simulation set-up
min_p_t <- 0.0001
max_it <- 50000
seed0 <- 922125
t_max <- 100
reps <- 100

# Probability mass function
f_t <- function(t, c, m) {
  1 / (1 + c * (t - m)^2)
}

probs_t <- function(t, p1, p2, p4, c, m) {
  ft <- f_t(t, c, m)
  p1t <- p1 * ft
  p2t <- p2 * ft
  p3t <- 1 - p1t - p2t - p4
  cbind(p1t, p2t, p3t, p4)
}
# probs_t(0.05,p1 = p1.true, p2 = p2.true, p4 =p4.true, c = c.true, m = m.true)
# p1.true/(1+c.true*(0.05-m.true)^2)

# Simulate data
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
                     p4 = p4.true, c = c.true, m = m.true)
    
    if ((sum(probs < min_p_t) > 0 | (itnum > max_it))) {
      cat(q_t, "\n")
      break
    }
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
  sim.markov <- rbind(sim.markov,cbind(reps = rep(i, length(sc.steps)), 
                                       time.steps, sc.steps, ec.steps, dc.steps))
}

index.test <- sample(1:100, size = 1) # test the optimization on one replications
index.test
sim.markov %>% filter(reps == index.test) %>%
  ggplot() +
  geom_point(aes(x = time.steps, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
  geom_point(aes(x = time.steps, y = sc.steps, col = "viable stem cell"), alpha = 0.5)+
  geom_point(aes(x = time.steps, y = dc.steps, col = "non-viable stem cell"), alpha = 0.5)+
  labs(x = "time (in hours)", y = "cell counts", color = "",
       title = "Simulated cell counts")+
  scale_color_manual(values = c("differentiated cell" = "skyblue3",
                                "viable stem cell" = "indianred1",
                                "non-viable stem cell" = "lightgreen",
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

sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == index.test]),
                       Z = diff(sim.markov$dc.steps[sim.markov$reps == index.test]),
                       t = sim.markov$time.steps[sim.markov$reps == index.test][-1],
                       deltat = diff(sim.markov$time.steps[sim.markov$reps == index.test]),
                       Si = sim.markov$sc.steps[sim.markov$reps == index.test][-length(sim.markov$sc.steps[sim.markov$reps == index.test])])


### Optimization
neg_log_lik <- function(par, data) {
  X <- data$X
  Z <- data$Z
  t <- data$t
  deltat <- data$deltat
  Si <- data$Si
  
  p1 <- par[1]
  p2 <- par[2]
  p4 <- par[3]
  c  <- par[4]
  m  <- par[5]
  r  <- par[6]
  
  ft <- 1 / (1 + c * (t - m)^2)
  
  # Parameter constraints
  if (p1 <= 0 || p2 <= 0 || p4 <= 0 || c <= 0 || (p1 + p2) >= 1 || (p1 + p2) * max(ft) + p4 >= 1 || r <= 0)
    return(Inf)
  
  loglik <- numeric(length(X))
  
  # indicator
  is_x1_z0 <- (X == 1 & Z == 0)
  is_x0_z0 <- (X == 0 & Z == 0)
  is_x0_z1 <- (X == 0 & Z == 1)
  is_xm1_z0 <- (X == -1 & Z == 0)
  
  loglik[is_x1_z0] <- log(p1 * ft[is_x1_z0])
  loglik[is_x0_z0] <- log(p2 * ft[is_x0_z0])
  loglik[is_x0_z1] <- log(p4)
  loglik[is_xm1_z0] <- log(1 - (p1 + p2) * ft[is_xm1_z0] - p4)
  
  loglik_exp <- log(r * Si) - r * Si * deltat
  
  return(-sum(loglik) - sum(loglik_exp))
}

grad_log_lik <- function(par, data) {
  X <- data$X
  Z <- data$Z
  t <- data$t
  deltat <- data$deltat
  Si <- data$Si
  
  p1 <- par[1]
  p2 <- par[2]
  p4 <- par[3]
  c  <- par[4]
  m  <- par[5]
  r  <- par[6]
  
  ft <- 1 / (1 + c * (t - m)^2)
  one_minus_p_ft_p4 <- 1 - (p1 + p2) * ft - p4
  
  if (p1 <= 0 || p2 <= 0 || p4 <= 0 || c <= 0 || (p1 + p2) >= 1 || any(one_minus_p_ft_p4 <= 0) || r <= 0)
    return(rep(NA, 6))
  
  I1 <- as.numeric(X == 1 & Z == 0)
  I2 <- as.numeric(X == 0 & Z == 0)
  I3 <- as.numeric(X == 0 & Z == 1)
  I4 <- as.numeric(X == -1 & Z == 0)
  
  df_dc <- - (t - m)^2 * ft^2
  df_dm <- 2 * c * (t - m) * ft^2
  
  denom <- one_minus_p_ft_p4
  denom[denom == 0] <- 1e-10
  
  dL_dp1 <- sum(I1 / p1 - I4 * ft / denom)
  dL_dp2 <- sum(I2 / p2 - I4 * ft / denom)
  dL_dp4 <- sum(I3 / p4 - I4 / denom)
  
  temp_c <- ( (I1 + I2) / ft - I4 * (p1 + p2) / denom )
  dL_dc <- sum(temp_c * df_dc)
  dL_dm <- sum(temp_c * df_dm)
  
  dL_dr <- sum(1 / r - Si * deltat)
  
  -c(dL_dp1, dL_dp2, dL_dp4, dL_dc, dL_dm, dL_dr)
}

ui <- matrix(c(
  1, 0, 0, 0, 0, 0,   # p1 > 0
  0, 1, 0, 0, 0, 0,   # p2 > 0
  0, 0, 1, 0, 0, 0,   # p4 > 0
  0, 0, 0, 1, 0, 0,   # c > 0
  -1, -1, -1, 0, 0, 0,  # p1 + p2 + p4 < 1
  0, 0, 0, 0, 0, 1    # r > 0
), byrow = TRUE, nrow = 6)

ci <- c(0, 0, 0, 0, -1, 0)

start <- c(p1 = 0.3, p2 = 0.3, p4 = 0.1, c = 10, m = 10, r = 0.01)
result <- constrOptim(
  theta = start,
  f = neg_log_lik,
  grad = grad_log_lik,
  ui = ui,
  ci = ci,
  data = sim.data
)

result$par


###
###
### estimate on all replications
start.time <- Sys.time()
# start <- c(0.1, 0.1, 0.1, 5, 5, 1)
start <- c(p1 = 0.2, p2 = 0.4, p4 = 0.1, c = 20, m = 20, r = 1)

res.prm <- data.frame(p1 = numeric(), p2 = numeric(), p4 = numeric(), c = numeric(),
                      m = numeric(), r = numeric())

for (j in 1:reps) {
  sim.data <- data.frame(X = diff(sim.markov$sc.steps[sim.markov$reps == j]),
                         Z = diff(sim.markov$dc.steps[sim.markov$reps == j]),
                         t = sim.markov$time.steps[sim.markov$reps == j][-1],
                         deltat = diff(sim.markov$time.steps[sim.markov$reps == j]),
                         Si = sim.markov$sc.steps[sim.markov$reps == j]
                         [-length(sim.markov$sc.steps[sim.markov$reps == j])])
  result <- constrOptim(
    theta = start,
    f = neg_log_lik,
    grad = grad_log_lik,
    ui = ui,
    ci = ci,
    data = sim.data
  )
  res.prm <- rbind(res.prm, result$par)
}
colnames(res.prm) <- c("p1", "p2", "p4", "c", "m", "r")
dim(res.prm)

end.time <- Sys.time()

end.time-start.time

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

# plot probablity mass function
p1_fun <- function(t) p1.true / (1 + c.true * (t - m.true)^2)
p2_fun <- function(t) p2.true / (1 + c.true * (t - m.true)^2)
p4_fun <- function(t) p4.true
p3_fun <- function(t) 1 - p1_fun(t) - p2_fun(t) - p4.true

t_vals <- seq(0, 60, length.out = 500)
data <- data.frame(
  t = t_vals,
  p1 = p1_fun(t_vals),
  p2 = p2_fun(t_vals),
  p3 = p3_fun(t_vals),
  p4 = p4_fun(t_vals)
)

data_long <- pivot_longer(data, cols = c(p1, p2, p3, p4), names_to = "Function", values_to = "Value")

ggplot(data_long, aes(x = t, y = Value, color = Function)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("skyblue3", "indianred3", "darkgreen", "black"),
                     labels = c(expression(p[1](t)), expression(p[2](t)), expression(p[3](t)),expression(p[4](t)))) +
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

