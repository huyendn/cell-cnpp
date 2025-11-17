library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)
library(DEoptim)


# True parameters
p1.true <- 0.55
c1.true <- 0.005
m1.true <- 4


p2.true <- 0.15
c2.true <- 0.012
m2.true <- 18

p4.true <- 0.2
c4.true <- 0.008
m4.true <- 28

s0.true <- 200
r.true <- 0.35

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
p4_fun <- function(t) p4.true / (1 + c4.true * (t - m4.true)^2)
p3_fun <- function(t) 1 - p1_fun(t) - p2_fun(t)

t_vals <- seq(0, 60, length.out = 500)
data <- data.frame(
  t = t_vals,
  p1 = p1_fun(t_vals),
  p2 = p2_fun(t_vals),
  p4 = p4_fun(t_vals),
  p3 = p3_fun(t_vals)
)

data_long <- pivot_longer(data, cols = c(p1, p2, p4, p3), names_to = "Function", values_to = "Value")

ggplot(data_long, aes(x = t, y = Value, color = Function)) +
  geom_line(linewidth = 1.2) +
  geom_hline(aes(yintercept = 1)) +
  scale_color_manual(values = c("skyblue3", "indianred3", "darkgreen", "orange3"),
                     labels = c(expression(p[1](t)), expression(p[2](t)), 
                                expression(p[3](t)), expression(p[4](t)))) +
  labs(
    title = "Functions p1(t), p2(t), p3(t), and p4(t)",
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

f_t4 <- function(t, c4, m4) {
  1 / (1 + c4 * (t - m4)^2)
}

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
    probs <- probs_t(t = time.current, 
                     p1 = p1.true, p2 = p2.true, p4 = p4.true, 
                     c1 = c1.true, c2 = c2.true, c4 = c4.true,
                     m1 = m1.true, m2 = m2.true, m4 = m4.true)
    
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
###
### plot one replication

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

###
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
  ggplot(aes(x = "", y = Value)) +
  geom_violin(aes(fill = Variable), color = "black") +
  geom_boxplot(fill = NA, width = 0.5)+
  facet_wrap(~ Variable, scale = "free_y", ncol = 3) +
  geom_jitter(width = 0.1, alpha = 0.3) + 
  theme_minimal() +
  labs(title = "Violin Plot for Each Parameter Estimates", x = "",
       fill= "Parameters", y = "Estimates")

###
### choosing a stop time
r_estimator_stop <- function(data, t_max){
  Si <- data$Si
  S_last <- Si[nrow(data)]
  deltat <- data$deltat
  t_last <- data$t[nrow(data)]
  n <- nrow(data)
  r <- n/(sum(Si*deltat) + (t_max - t_last)*S_last)
  return(r)
}

min(sim.markov %>% filter(sc.steps == 0) %>% pull(time.steps))

t.stop <- 40

sim.markov.stop <- sim.markov %>% filter(time.steps <= t.stop)

res.prm <- data.frame()
start_time <- Sys.time()
for (j in 1:reps) {
  sim.data <- data.frame(X = diff(sim.markov.stop$sc.steps[sim.markov.stop$reps == j]),
                         Z = diff(sim.markov.stop$dc.steps[sim.markov.stop$reps == j]),
                         t = sim.markov.stop$time.steps[sim.markov.stop$reps == j][-1],
                         deltat = diff(sim.markov.stop$time.steps[sim.markov.stop$reps == j]),
                         Si = sim.markov.stop$sc.steps[sim.markov.stop$reps == j]
                         [-length(sim.markov.stop$sc.steps[sim.markov.stop$reps == j])])
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
  r.est <- r_estimator_stop(data = sim.data, t_max = t.stop)
  res.prm <- rbind(res.prm, c(par.est,r.est))
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
  ggplot(aes(x = "", y = Value)) +
  geom_violin(aes(fill = Variable), color = "black") +
  geom_boxplot(fill = NA, width = 0.5)+
  facet_wrap(~ Variable, scale = "free_y", ncol = 3) +
  geom_jitter(width = 0.1, alpha = 0.3) + 
  theme_minimal() +
  labs(title = "Violin Plot for Each Parameter Estimates", x = "",
       fill= "Parameters", y = "Estimates")



library(statmod)
library(fitdistrplus)
library(goftest)

head(sim.markov)

sim.markov %>% filter(sc.steps == 0) %>%
  ggplot(aes(x = time.steps)) +
  geom_histogram(
    aes(y = after_stat(count)/sum(after_stat(count))),
    binwidth = 1.5, 
    color = "black",fill = "skyblue3") +
  scale_x_continuous(limits = c(20, 100))+
  scale_y_continuous(limits = c(0, 0.15)) +
  labs(x = "time", y = "relative frequency",
       title = "Distribution of stopping time",
       subtitle = paste0(reps, " replications with ", s0.true, " starting stem cells")) +
  theme_minimal()

stop.time <- sim.markov %>% filter(sc.steps == 0) %>% pull(time.steps)
min.stop.time <- 1/r.true*log(s0.true) - digamma(1)
stop.time.shift <- stop.time - min.stop.time
min.stop.time
fit <- fitdist(stop.time.shift, "invgauss", 
               start = list(mean = mean(stop.time.shift), shape = 1))
summary(fit)
plot(fit)


ks.test(stop.time.shift, "pinvgauss", 
        mean = fit$estimate["mean"], 
        shape = fit$estimate["shape"])

ad.test(stop.time.shift, null = "pinvgauss", 
        mean = fit$estimate["mean"], shape = fit$estimate["shape"])

fit_gamma <- fitdist(stop.time.shift, "gamma")
fit_lognorm <- fitdist(stop.time.shift, "lnorm")

gofstat(list(fit, fit_gamma, fit_lognorm))