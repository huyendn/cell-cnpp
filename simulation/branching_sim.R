library(ggplot2)
library(dplyr)
library(plotly)

q_t1 <- function(tms, prm) {
  return(2 - 1/(prm[1] + prm[2]*tms + prm[3]*tms^2))
}

## model parameters
s0.test <- 200
r.test <- 0.2
prm.test <-c(1.2, -0.106, 0.005)

# simulations settings
min_q_t = 0.0001
max_it = 50000
k = 3
seed0 = 922125
t_max = 100
reps <- 100

# expected differentiated cell function
ec_func <- function(t){
  integrand <- function(x){
    out <- q_t1(tms = x, prm = prm.test)*intensity_func(x)
  }
  out.ec <- integrate(integrand, lower = 0, upper = t)$value
  return(out.ec)
}

# expected stem cell function
sc_func <- function(t){
  out <- s0.test*exp(r.test*(t-int_q(t,prm.test)+int_q(0,prm.test)))
  return(out)
}

# variance differentiated cell function
var_ec_func <- function(t){
  integrand <- function(x){
    q_t <- q_t1(tms = x, prm = prm.test)
    out <- ((k-2)*(2-q_t)/k + 4*(2-k + (k-1)*q_t)/k)*intensity_func(x)
    return(out)
  }
  out.var <- integrate(integrand, lower = 0, upper = t)$value
  return(out.var)
}
# variance stem cell function
var_sc_func <- function(t){
  integrand <- function(x){
    q_t <- q_t1(tms = x, prm = prm.test)
    out <- ((2-q_t)/k + (2-k +(k-1)*q_t)/k)*intensity_func(x)
    return(out)
  }
  out.var.sc <- integrate(integrand, lower = 0, upper = t)$value
  return(out.var.sc)
}

# covariance 
cov_func <- function(t){
  integrand <- function(x){
    q_t <- q_t1(tms = x, prm = prm.test)
    out <- (-2*(2-k + (k-1)*q_t)/k)*intensity_func(x)
    return(out)
  }
  out.cov <- integrate(integrand, lower = 0, upper = t)$value
  return(out.cov)
}

set.seed(seed0)

sim.markov <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric())

for (i in 1:reps) {
  sc.current <- s0.test
  ec.current <- dc.current <- time.current <- 0
  sc.steps <- c(s0.test)
  ec.steps <- c(ec.current)
  dc.steps <- c(dc.current)
  time.steps <- c(time.current)
  itnum <- 0
  while (sc.current > 0) {
    itnum <- itnum + 1
    time.current <- time.current + rexp(1, rate = r.test*sc.current)
    q_t <- q_t1(time.current, prm.test)
    if ((q_t < min_q_t| (itnum > max_it))) {
      cat(q_t, "\n")
      break
    }
    probs <- c((2-q_t)/k, max(0, (k-2)*(2-q_t)/k), max(0,(2-k + (k-1)*q_t)/k))
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


tms <- seq(0, 80, by = 5)

dat.markov <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric(),
                         time.collected = numeric())

for (i in 1:reps) {
  subset <- sim.markov %>% filter(reps == i)
  for (j in 2:length(tms)) {
    time.before <- tms[j-1]
    time.collected <- tms[j]
    if (sum(subset$time <= time.collected  & subset$time > time.before) == 0){
      # dat[nrow(dat)+1, ] <- c(dat[nrow(dat), 1:4], time.collected)
      dat.markov <- dat.markov
    } else{
      ind <- max(which(subset$time <= time.collected  & subset$time > time.before))
      dat.markov <- rbind(dat.markov, cbind(subset[ind,], time.collected))
    }
    
  }
}

ggplot() +
  geom_point(data = dat.markov, aes(x = time.steps, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
  geom_point(data = dat.markov, aes(x = time.steps, y = sc.steps, col = "stem cell"), alpha = 0.5)+
  # geom_point(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  # geom_point(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  # geom_line(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  # geom_line(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  # geom_ribbon(data = theo.dat, aes(x = time.steps, 
  #                                  ymin = sc.steps-1.96*sqrt(sc.var), 
  #                                  ymax = sc.steps+1.96*sqrt(sc.var)), color = "gray", alpha = 0.25)+
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


###
pmf_t <- function(x, t, params) {
  a <- params[1]
  b <- params[2]
  c <- params[3]
  k <- params[4]
  
  # Compute q(t)
  q <- 2 - 1 / (a + b * t + c * t^2)
  
  # Compute raw (unnormalized) probabilities
  p1   <- (2 - q) / k
  p0   <- max(0, (k - 2) * (2 - q) / k)
  p_1  <- max(0, (2 - k + (k - 1) * q) / k)
  
  probs_raw <- c(p1, p0, p_1)
  support   <- c(1, 0, -1)
  
  # Normalize to sum to 1
  probs <- probs_raw / sum(probs_raw)
  
  # Return probability for this x
  prob <- probs[match(x, support)]
  if (is.na(prob)) prob <- 0
  return(prob)
}

pmf_vector <- function(x_vec, t_vec, params) {
  if (length(x_vec) != length(t_vec)) {
    stop("x_vec and t_vec must be the same length.")
  }
  
  # Vectorized call using mapply
  sapply(seq_along(x_vec), function(i) pmf_t(x_vec[i], t_vec[i], params))
}

# test function
x_vec <- c(1, 0, -1, 1, -1)
t_vec <- c(0, 2, 4, 6, 8)
params <- c(a = 1, b = 0.5, c = 0.1, k = 3)
pmf_vector(x_vec, t_vec, params)


get_interarrival_times <- function(t_vec) {
  return(diff(t_vec))
}


joint_prob <- function(x_vec, t_vec, S_vec, params, r) {
  n <- length(x_vec)
  
  if (length(t_vec) != (n + 1) || length(S_vec) != (n + 1)) {
    stop("t_vec and S_vec must be of length length(x_vec) + 1")
  }
  
  interarrival_times <- diff(t_vec)
  probs <- numeric(n)
  
  for (i in 1:n) {
    S_i <- S_vec[i]
    delta_t <- interarrival_times[i]
    x_i <- x_vec[i]
    t_i <- t_vec[i+1]  
    
    # 1. Transition probability P(x_i | t_i)
    px <- pmf_t(x_i, t_i, params)
    
    # 2. Interarrival time density P(Δt_i | S_i)
    if (S_i <= 0 || r <= 0) {
      pt <- 1e-12
    } else {
      pt <- dexp(delta_t, rate = r * S_i)
    }
    
    # 3. Joint probability
    probs[i] <- px * pt
  }
  return(probs)
}


log_likelihood_joint <- function(params_with_r, x_vec, t_vec, S_vec) {
  params <- params_with_r[1:4]
  r <- params_with_r[5]
  
  if (r <= 0) return(-Inf)
  
  probs <- joint_prob(x_vec, t_vec, S_vec, params, r)
  probs[probs < 1e-12] <- 1e-12
  
  return(-sum(log(probs)))  # for minimization
}

# Simulated data
# S_vec <- c(2, 1, 2, 1, 1.5)
# x_vec <- diff(S_vec)
# t_vec <- c(0, 0.5, 1.5, 2.7, 4.0) 

S_vec <- sim.markov$sc.steps[sim.markov$reps == 1]
x_vec <- diff(S_vec)
t_vec <- sim.markov$time.steps[sim.markov$reps == 1]

# Initial value: a, b, c, k, r
init_params <- c(1, 0.5, 0.5, 5, 1)

# Lower bounds: a > 0.5, b free, c ≥ 0, k > 0, r > 0
lower_bounds <- c(a = 0.5, b = -100, c = 0, k = 1e-6, r = 1e-6)

# Upper bounds (set large enough to not interfere)
upper_bounds <- c(a = 100, b = 100, c = 100, k = 10, r = 10)

# Optimization
fit <- optim(par = init_params,
             fn = log_likelihood_joint,
             x_vec = x_vec,
             t_vec = t_vec,
             S_vec = S_vec,
             method = "L-BFGS-B",
             lower = lower_bounds,
             upper = upper_bounds)


names(init_params) <- c("a", "b", "c", "k", "r")
names(lower_bounds) <- names(init_params)
names(upper_bounds) <- names(init_params)
fit$par

res.prm <- data.frame(a = numeric(), b = numeric(), c = numeric(),
                      k = numeric(), r = numeric())
res.fn.value <- c()

for (j in 1:reps) {
  S_vec <- sim.markov$sc.steps[sim.markov$reps == j]
  x_vec <- diff(S_vec)
  t_vec <- sim.markov$time.steps[sim.markov$reps == j]
  fit <- optim(par = init_params,
               fn = log_likelihood_joint,
               x_vec = x_vec,
               t_vec = t_vec,
               S_vec = S_vec,
               method = "L-BFGS-B",
               lower = lower_bounds,
               upper = upper_bounds)
  res.prm <- rbind(res.prm, fit$par)
  res.fn.value <- c(res.fn.value, fit$value)
}

