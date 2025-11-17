library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
s0.true <- 200
r.true <- 0.2

p1.true <- 0.5
c1.true <- 0.005
m1.true <- 4

p2.true <- 0.2
c2.true <- 0.01
m2.true <- 12

p1.func <- function(t){
  p1.prob <- p1.true/(1+c1.true*(t-m1.true)^2)
  return(p1.prob)
}

p2.func <- function(t){
  p2.prob <- p2.true/(1+c2.true*(t-m2.true)^2)
  return(p2.prob)
}

p3.func <- function(t){
  p3.prob <- 1- p1.true/(1+c1.true*(t-m1.true)^2)-p2.true/(1+c2.true*(t-m2.true)^2)
  return(p3.prob)
}

int.p1.func <- function(t){
  (p1.true/sqrt(c1.true))*(atan(sqrt(c1.true)*(t-m1.true))+atan(sqrt(c1.true)*m1.true))
}

int.p2.func <- function(t){
  (p2.true/sqrt(c2.true))*(atan(sqrt(c2.true)*(t-m2.true))+atan(sqrt(c2.true)*m2.true))
}

int.p3.func <- function(t){
  t-int.p1.func(t) - int.p2.func(t)
}


# sc.mean <- function(t){
#   integrand <- function(x){
#     p1.func(x) - p3.func(x)
#   }
#   integral.out <- integrate(integrand, 0, t)$value
#   return(as.numeric(s0.true*exp(r.true*integral.out)))
# }

sc.mean <- function(t){
  s0.true*exp(r.true*(int.p1.func(t)-int.p3.func(t)))
}

sc.var.cnpp <- function(t){
  integrand <- function(x){
    (p1.func(x) + p3.func(x)) * r.true * sapply(x, sc.mean)
  }
  integral.out <- integrate(integrand, 0, t)$value
  return(integral.out)
}

sc.var.bp <- function(t) {
  A <- exp(2*r.true*(int.p1.func(t)-int.p3.func(t)))
  integrand_B <- function(x){
    (p1.func(x) + p3.func(x))*exp(-r.true*(int.p1.func(x)-int.p3.func(x)))
  }
  B <- integrate(integrand_B, lower = 0, upper = t)$value
  return(s0.true*r.true*A * B)
}

# sc.var.bp <- function(t) {
#   integrand_A <- function(u) {
#     2 * r.true * (p1.func(u) - p3.func(u))
#   }
#   A <- exp(integrate(integrand_A, lower = 0, upper = t)$value)
#   
#   integrand_B <- function(u_vec) {
#     sapply(u_vec, function(u) {
#       inner_exp_integrand <- function(v) {
#         p1.func(v) - p3.func(v)
#       }
#       inner_int <- integrate(inner_exp_integrand, lower = 0, upper = u)$value
#       
#       (p1.func(u) + p3.func(u)) * r.true * s0.true * exp(-r.true * inner_int)
#     })
#   }
#   
#   B <- integrate(integrand_B, lower = 0, upper = t)$value
#   # result: A(t) * B(t)
#   return(A * B)
# }

sc.autocorr.bp <- function(t, u) {
  if (t >= u) {
    return(sc.mean(t) / sc.mean(u) * sc.var.bp(u))
  } else {
    return(sc.mean(u) / sc.mean(t) * sc.var.bp(t))
  }
}

time.input <- seq(0, 60, by = 0.1)

mean.output <- sapply(time.input, sc.mean)
var.cpnn.output <- sapply(time.input, sc.var.cnpp)
var.bp.output <- sapply(time.input, sc.var.bp)

theoretical.dat <- data.frame(time = time.input,
                              mean = mean.output,
                              var.cpnn = var.cpnn.output,
                              var.bp = var.bp.output)

ggplot(theoretical.dat) +
  geom_line(aes(x = time, y = var.bp, col = "bp"), linetype = "dashed") +
  geom_line(aes(x = time, y = var.cpnn, col = "cpnn")) +
  labs(x = "time", y = "variance", title = "variance comparison", color = "") +
  scale_color_manual(values = c("cpnn" = "blue", "bp" = "red")) +
  theme_minimal() + theme(legend.position = "bottom")

ggplot(theoretical.dat) +
  geom_line(data = theoretical.dat, aes(x = time, y = mean), 
            col = "red", linewidth = 1.5) +
  geom_ribbon(data = theoretical.dat, 
              aes(x = time, ymax = mean + 2*sqrt(var.bp), 
                  ymin = mean-2*sqrt(var.bp),fill = "bp"),
              alpha = 0.2) +
  geom_ribbon(data = theoretical.dat, 
              aes(x = time, ymax = mean + 2*sqrt(var.cpnn), 
                  ymin = mean-2*sqrt(var.cpnn),fill = "cpnn"),
              alpha = 0.2)+
  labs(x = "time", y = "cell count", title = "two standard deviation from the mean", fill = "") +
  scale_fill_manual(values = c("cpnn" = "blue", "bp" = "red")) +
  theme_minimal() + theme(legend.position = "bottom")


# Simulation set-up
min_p_t <- 0.0001
max_it <- 50000
seed0 <- 922125
t_max <- 100
reps <- 100

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

sim.markov %>% filter(reps <= 50) %>% filter(time.steps <= 60) %>%
  ggplot() +
  geom_point(aes(x = time.steps, y = sc.steps), size = 0.2, alpha = 0.05)+
  geom_line(data = theoretical.dat, aes(x = time, y = mean), 
            col = "red", linewidth = 1.5) +
  geom_ribbon(data = theoretical.dat, 
              aes(x = time, ymax = mean + 2*sqrt(var.bp), ymin = mean-2*sqrt(var.bp)),
              fill = "red", alpha = 0.3) +
  # geom_ribbon(data = theoretical.dat, 
  #             aes(x = time, ymax = mean + 2*sqrt(var.cpnn), ymin = mean-2*sqrt(var.cpnn)),
  #             fill = "blue", alpha = 0.2) +
  labs(x = "time", y = "cell counts", title = "BP simulated cell counts",
       subtitle = "theoretical mean and bp variance") +
  theme_minimal() 

sim.markov %>% filter(reps <= 50) %>% filter(time.steps <= 60) %>%
  ggplot() +
  geom_point(aes(x = time.steps, y = sc.steps), size = 0.2, alpha = 0.05)+
  geom_line(data = theoretical.dat, aes(x = time, y = mean), 
            col = "red", linewidth = 1.5) +
  # geom_ribbon(data = theoretical.dat, 
  #             aes(x = time, ymax = mean + 2*sqrt(var.bp), ymin = mean-2*sqrt(var.bp)),
  #             fill = "red", alpha = 0.3) +
  geom_ribbon(data = theoretical.dat,
              aes(x = time, ymax = mean + 2*sqrt(var.cpnn), ymin = mean-2*sqrt(var.cpnn)),
              fill = "blue", alpha = 0.3) +
  labs(x = "time", y = "cell counts", title = "BP simulated cell counts",
       subtitle = "theoretical mean and cpnn variance") +
  theme_minimal() 

tms <- seq(0, 60, by = 5)

dat_list <- vector("list", reps)

for (i in seq_len(reps)) {
  sim_sub <- sim.markov %>% filter(reps == i)
  dat_i <- data.frame()
  
  for (j in 2:length(tms)) {
    time.before <- tms[j - 1]
    time.collected <- tms[j]
    
    # Find observations in the time window
    valid_idx <- which(sim_sub$time <= time.collected & sim_sub$time > time.before)
    
    if (length(valid_idx) > 0) {
      ind <- max(valid_idx)  # last observation in the window
      dat_i <- bind_rows(dat_i, cbind(sim_sub[ind, ], time.collected))
    }
  }
  
  dat_list[[i]] <- dat_i
}

dat <- bind_rows(dat_list)

ggplot() +
  # geom_point(data = dat, aes(x = time.collected, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
  geom_point(data = dat, aes(x = time.steps, y = sc.steps), alpha = 0.5)+
  # geom_point(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_point(data = theoretical.dat, aes(x = time, y = mean, col = "expected stem cell"))+
  # geom_line(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  # geom_line(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_ribbon(data = theoretical.dat, aes(x = time,
                                   ymin = mean-2*sqrt(var.bp),
                                   ymax = mean+2*sqrt(var.bp)), fill = "red", alpha = 0.25)+
  geom_ribbon(data = theoretical.dat, aes(x = time,
                                          ymin = mean-2*sqrt(var.cpnn),
                                          ymax = mean+2*sqrt(var.cpnn)), fill = "blue", alpha = 0.25)+
  # geom_ribbon(data = theo.dat, aes(x = time.steps, ymin = ec.steps-1.96*sqrt(ec.var),
  #                                  ymax = ec.steps+1.96*sqrt(ec.var)), color = "gray", alpha = 0.25)+
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

### check autocorrelation formula
tms.cor <- seq(5,25, by = 5)

# Create a matrix of E[X(t) X(u)]
autocorr.matrix <- outer(tms.cor, tms.cor, Vectorize(sc.autocorr.bp))

colnames(autocorr.matrix) <- tms.cor
rownames(autocorr.matrix) <- tms.cor

theosd_vec <- sqrt(diag(autocorr.matrix))
theo_corr <- autocorr.matrix / (theosd_vec %o% theosd_vec)

print(round(autocorr.matrix, 4))
print(round(theo_corr, 4))
# df <- melt(autocorr.matrix)
# colnames(df) <- c("t", "u", "E_XtXu")
# 
# ggplot(df, aes(x = t, y = u, fill = E_XtXu)) +
#   geom_tile() +
#   scale_fill_gradient2(
#     low = "#4575b4", mid = "white", high = "#d73027",
#     midpoint = 0,  # center of the scale
#     name = "E[X(t)X(u)]"
#   ) +
#   labs(title = "Theoretical E[X(t)X(u)]", x = "t", y = "u") +
#   theme_minimal()

dat_wide <- dat %>% filter(time.collected %in% tms.cor)%>%
  select(reps, time.collected, sc.steps) %>%
  pivot_wider(names_from = time.collected, values_from = sc.steps)

X <- as.matrix(dat_wide[,-1])  # remove 'reps' column
R <- nrow(X)

# Center data
X_centered <- sweep(X, 2, colMeans(X))
empirical_cov <- (t(X_centered) %*% X_centered) / R

sd_vec <- sqrt(diag(empirical_cov))
empirical_corr <- empirical_cov / (sd_vec %o% sd_vec)

colnames(empirical_corr) <- time_points
rownames(empirical_corr) <- time_points
print(round(empirical_corr, 3))
