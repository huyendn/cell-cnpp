library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)

# True parameters
p1.true <- 0.2
c1.true <- 0.005
m1.true <- 4

p2.true <- 0.5
c2.true <- 0.005
m2.true <- 12

parm.true <- c(p1.true, c1.true, m1.true, p2.true, c2.true, m2.true)
s0.true <- 200
r.true <- 0.2

# Simulation set-up
min_p_t <- 0.0001
max_it <- 50000
seed0 <- 922125
t_max <- 100
reps <- 100


# plot functions
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


# Probability mass functions for each component
f_t1 <- function(t, c1, m1) {
  1 / (1 + c1 * (t - m1)^2)
}

f_t2 <- function(t, c2, m2) {
  1 / (1 + c2 * (t - m2)^2)
}

# Combined probabilities
probs_t <- function(t, p1, p2, c1, c2, m1, m2) {
  f1t <- f_t1(t, c1, m1)
  f2t <- f_t2(t, c2, m2)
  
  p1t <- p1 * f1t
  p2t <- p2 * f2t
  p3t <- 1 - p1t - p2t
  
  # Optionally check validity
  if (any(p3t < 0)) warning("Some p3(t) values are negative — check parameters!")
  
  cbind(p1t, p2t, p3t)
}


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
head(sim.data)

###
###
neg_log_likelihood <- function(params, data) {
  # Parameters
  p1 <- params[1]  # between 0 and 1
  c1 <- params[2]  # > 0
  m1 <- params[3]
  
  p2 <- params[4]  # between 0 and 1
  c2 <- params[5]  # > 0
  m2 <- params[6]
  
  # Extract data
  x <- data$X
  t <- data$t
  
  # Compute time-dependent probabilities
  f1 <- p1 / (1 + c1 * (t - m1)^2)
  f0 <- p2 / (1 + c2 * (t - m2)^2)
  f_neg1 <- 1 - f1 - f0
  
  # Check for invalid probabilities
  if (any(f1 < 0 | f1 > 1 | f0 < 0 | f0 > 1 | f_neg1 < 0 | f_neg1 > 1)) {
    return(Inf)
  }
  
  # Assign probabilities based on observed x
  epsilon <- 1e-12
  prob <- ifelse(x == 1, f1,
                 ifelse(x == 0, f0,
                        f_neg1))
  
  # Avoid log(0)
  prob <- pmax(prob, epsilon)
  
  # Negative log-likelihood
  nll <- -sum(log(prob))
  return(nll)
}

grad_neg_log_likelihood <- function(params, data) {
  # Parameters
  p1 <- params[1]; c1 <- params[2]; m1 <- params[3]
  p2 <- params[4]; c2 <- params[5]; m2 <- params[6]
  
  # Data
  t <- data$t
  x <- data$X
  
  # Precompute
  delta1 <- t - m1
  delta2 <- t - m2
  
  A <- 1 + c1 * delta1^2
  B <- 1 + c2 * delta2^2
  
  f1 <- p1 / A
  f0 <- p2 / B
  f_neg1 <- 1 - f1 - f0
  
  # Probabilities for each x
  prob <- ifelse(x == 1, f1,
                 ifelse(x == 0, f0,
                        f_neg1))
  epsilon <- 1e-12
  prob <- pmax(prob, epsilon)
  
  # Indicator vectors
  I1 <- as.numeric(x == 1)
  I0 <- as.numeric(x == 0)
  Ineg1 <- as.numeric(x == -1)
  
  # Compute partial derivatives
  df1_dp1 <- 1 / A
  df1_dc1 <- -p1 * delta1^2 / A^2
  df1_dm1 <- -2 * p1 * c1 * delta1 / A^2
  
  df0_dp2 <- 1 / B
  df0_dc2 <- -p2 * delta2^2 / B^2
  df0_dm2 <- -2 * p2 * c2 * delta2 / B^2
  
  dfneg1_dp1 <- -df1_dp1
  dfneg1_dc1 <- -df1_dc1
  dfneg1_dm1 <- -df1_dm1
  
  dfneg1_dp2 <- -df0_dp2
  dfneg1_dc2 <- -df0_dc2
  dfneg1_dm2 <- -df0_dm2
  
  # Select relevant derivative based on x
  dL_dp1 <- -(I1 * df1_dp1 + Ineg1 * dfneg1_dp1) / prob
  dL_dc1 <- -(I1 * df1_dc1 + Ineg1 * dfneg1_dc1) / prob
  dL_dm1 <- -(I1 * df1_dm1 + Ineg1 * dfneg1_dm1) / prob
  
  dL_dp2 <- -(I0 * df0_dp2 + Ineg1 * dfneg1_dp2) / prob
  dL_dc2 <- -(I0 * df0_dc2 + Ineg1 * dfneg1_dc2) / prob
  dL_dm2 <- -(I0 * df0_dm2 + Ineg1 * dfneg1_dm2) / prob
  
  # Sum over all observations
  grad <- c(
    sum(dL_dp1),
    sum(dL_dc1),
    sum(dL_dm1),
    sum(dL_dp2),
    sum(dL_dc2),
    sum(dL_dm2)
  )
  
  return(grad)
}

start <- c(p1 = 0.1, c1 = 0.1, m1 = 1, p2 = 0.4, c2 = 0.1, m2 = 1)
# Define the constraint matrix (ui) and rhs (ci)
ui <- rbind(
  c( 1,  0, 0,  0,  0, 0),  # p1 > 0.001
  c(-1,  0, 0,  0,  0, 0),  # p1 < 0.999
  c( 0,  1, 0,  0,  0, 0),  # c1 > 1e-6
  c( 0,  0, 0,  1,  0, 0),  # p2 > 0.001
  c( 0,  0, 0, -1,  0, 0),  # p2 < 0.999
  c( 0,  0, 0,  0,  1, 0),  # c2 > 1e-6
  c( 0,  0, 1,  0,  0, 0),  # m1 > 0
  c( 0,  0, 0,  0,  0, 1)  # m2 > 0
)

ci <- c(0.001, -0.999, 1e-6, 0.001, -0.999, 1e-6, 0, 0)


result <- constrOptim(
  theta = start,
  f = neg_log_likelihood,
  grad = grad_neg_log_likelihood,
  ui = ui,
  ci = ci,
  data = sim.data,
  control = list(reltol = 1e-8)
)

# Output results
cat("Optimized Parameters:\n")
print(result$par)
parm.true
cat("Final NLL:\n")
print(result$value)


###
###
### estimate on all replications
start.time <- Sys.time()
# start <- c(0.1, 0.1, 5, 5, 1)
start <- c(p1 = 0.1, c1 = 0.1, m1 = 1, p2 = 0.4, c2 = 0.1, m2 = 1)

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
    f = neg_log_likelihood,
    grad = grad_neg_log_likelihood,
    ui = ui,
    ci = ci,
    data = sim.data,
    control = list(reltol = 1e-8)
  )
  res.prm <- rbind(res.prm, result$par)
}
colnames(res.prm) <- c("p1","c1","m1", "p2", "c2", "m2")
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
