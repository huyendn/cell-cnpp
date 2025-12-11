library(tidyr)
library(dplyr)
library(ggplot2)

loglik_R <- function(params, data) {
  
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
  
  # hidden state set
  hidden_states <- 0
  alpha_prev <- 1
  
  # gradients of α wrt parameters
  alpha_prev_grad_p1 <- alpha_prev_grad_p2 <- alpha_prev_grad_p4 <- 0
  alpha_prev_grad_c1 <- alpha_prev_grad_c2 <- alpha_prev_grad_c4 <- 0
  alpha_prev_grad_m1 <- alpha_prev_grad_m2 <- alpha_prev_grad_m4 <- 0
  alpha_prev_grad_r  <- 0
  
  # log-likelihood gradient accumulators
  grad_p1 <- 0; grad_p2 <- 0; grad_p4 <- 0
  grad_c1 <- 0; grad_c2 <- 0; grad_c4 <- 0
  grad_m1 <- 0; grad_m2 <- 0; grad_m4 <- 0
  grad_r  <- 0
  
  total_loglik <- 0
  
  clamp <- function(x) max(min(x, 1 - 1e-12), 1e-12)
  
  for (step in 2:n_steps) {
    
    deltaT <- t[step] - t[step-1]
    deltaY <- Y[step] - Y[step-1]
    deltaM <- M[step] - M[step-1]
    
    # event probabilities
    p1_t <- clamp(p1 / (1 + c1*(t[step] - m1)^2))
    p2_t <- clamp(p2 / (1 + c2*(t[step] - m2)^2))
    p4_t <- clamp(p4 / (1 + c4*(t[step] - m4)^2))
    
    p_other <- 1 - p1_t - p2_t - p4_t
    p_other <- clamp(p_other)
    
    # exponential rates
    rates <- r * (M[step-1] - hidden_states)
    rates[!is.finite(rates)] <- 1e-8
    rates[rates < 1e-8] <- 1e-8
    
    exp_term <- dexp(deltaT, rates)
    if (any(!is.finite(exp_term))) return(safe_fail())
    
    # derivatives wrt parameters
    dh_dp1 <- 1/(1 + c1*(t[step] - m1)^2)
    dh_dp2 <- 1/(1 + c2*(t[step] - m2)^2)
    dh_dp4 <- 1/(1 + c4*(t[step] - m4)^2)
    
    dh_dc1 <- -p1*((t[step]-m1)^2)/((1 + c1*(t[step] - m1)^2)^2)
    dh_dc2 <- -p2*((t[step]-m2)^2)/((1 + c2*(t[step] - m2)^2)^2)
    dh_dc4 <- -p4*((t[step]-m4)^2)/((1 + c4*(t[step] - m4)^2)^2)
    
    dh_dm1 <- p1 * (-2 * c1 * (t[step] - m1)) / (1 + c1*(t[step] - m1)^2)^2
    dh_dm2 <- p2 * (-2 * c2 * (t[step] - m2)) / (1 + c2*(t[step] - m2)^2)^2
    dh_dm4 <- p4 * (-2 * c4 * (t[step] - m4)) / (1 + c4*(t[step] - m4)^2)^2
    
    # correct dh/dr
    lambda <- rates
    dh_dr <- (M[step-1] - hidden_states) * exp(-lambda * deltaT) * (1 - lambda * deltaT)
    
    ### ================================
    ### EVENT BRANCHES
    ### ================================
    
    if (deltaY == 1 && deltaM == 0) {  
      # p2 event
      alpha_new <- alpha_prev * exp_term * p2_t
      
      grad_term_p1 <- alpha_prev_grad_p1 * exp_term * p2_t
      grad_term_p2 <- alpha_prev_grad_p2 * exp_term * p2_t + alpha_prev * dh_dp2
      grad_term_p4 <- alpha_prev_grad_p4 * exp_term * p2_t
      
      grad_term_c1 <- alpha_prev_grad_c1 * exp_term * p2_t
      grad_term_c2 <- alpha_prev_grad_c2 * exp_term * p2_t + alpha_prev * dh_dc2
      grad_term_c4 <- alpha_prev_grad_c4 * exp_term * p2_t
      
      grad_term_m1 <- alpha_prev_grad_m1 * exp_term * p2_t
      grad_term_m2 <- alpha_prev_grad_m2 * exp_term * p2_t + alpha_prev * dh_dm2
      grad_term_m4 <- alpha_prev_grad_m4 * exp_term * p2_t
      
      grad_term_r <- alpha_prev_grad_r * exp_term * p2_t + alpha_prev * dh_dr
    }
    
    else if (deltaY == 2 && deltaM == -1) {  
      # p3 event
      alpha_new <- alpha_prev * exp_term * p_other
      
      grad_term_p1 <- alpha_prev_grad_p1 * exp_term * p_other + alpha_prev * (-dh_dp1)
      grad_term_p2 <- alpha_prev_grad_p2 * exp_term * p_other + alpha_prev * (-dh_dp2)
      grad_term_p4 <- alpha_prev_grad_p4 * exp_term * p_other + alpha_prev * (-dh_dp4)
      
      grad_term_c1 <- alpha_prev_grad_c1 * exp_term * p_other + alpha_prev * (-dh_dc1)
      grad_term_c2 <- alpha_prev_grad_c2 * exp_term * p_other + alpha_prev * (-dh_dc2)
      grad_term_c4 <- alpha_prev_grad_c4 * exp_term * p_other + alpha_prev * (-dh_dc4)
      
      grad_term_m1 <- alpha_prev_grad_m1 * exp_term * p_other + alpha_prev * (-dh_dm1)
      grad_term_m2 <- alpha_prev_grad_m2 * exp_term * p_other + alpha_prev * (-dh_dm2)
      grad_term_m4 <- alpha_prev_grad_m4 * exp_term * p_other + alpha_prev * (-dh_dm4)
      
      grad_term_r <- alpha_prev_grad_r * exp_term * p_other + alpha_prev * dh_dr
      
      # filter states
      valid <- (hidden_states <= M[step])
      if (!any(valid)) return(safe_fail())
      
      alpha_new     <- alpha_new[valid]
      hidden_states <- hidden_states[valid]
      
      grad_term_p1 <- grad_term_p1[valid]
      grad_term_p2 <- grad_term_p2[valid]
      grad_term_p4 <- grad_term_p4[valid]
      
      grad_term_c1 <- grad_term_c1[valid]
      grad_term_c2 <- grad_term_c2[valid]
      grad_term_c4 <- grad_term_c4[valid]
      
      grad_term_m1 <- grad_term_m1[valid]
      grad_term_m2 <- grad_term_m2[valid]
      grad_term_m4 <- grad_term_m4[valid]
      
      grad_term_r <- grad_term_r[valid]
    }
    
    else if (deltaY == 0 && deltaM == 1) {   
      # p1/p4 event
      diag_part  <- alpha_prev * exp_term * p1_t
      shift_part <- alpha_prev * exp_term * p4_t
      
      L <- length(diag_part)
      if (L == 0) return(safe_fail())
      
      # combine α
      alpha_new <- numeric(L + L)
      alpha_new[1:L] <- diag_part
      alpha_new[(L+1):(2*L)] <- shift_part
      
      # update hidden_states (stay + shift)
      hidden_states <- c(hidden_states, hidden_states + 1)
      
      # filter
      valid <- (hidden_states <= M[step])
      if (!any(valid)) return(safe_fail())
      
      alpha_new     <- alpha_new[valid]
      hidden_states <- hidden_states[valid]
      
      # gradients same pattern
      diag_grad_term_p1 <- alpha_prev_grad_p1 * exp_term * p1_t + alpha_prev * dh_dp1
      diag_grad_term_p2 <- alpha_prev_grad_p2 * exp_term * p1_t
      diag_grad_term_p4 <- alpha_prev_grad_p4 * exp_term * p1_t
      
      diag_grad_term_c1 <- alpha_prev_grad_c1 * exp_term * p1_t + alpha_prev * dh_dc1
      diag_grad_term_c2 <- alpha_prev_grad_c2 * exp_term * p1_t
      diag_grad_term_c4 <- alpha_prev_grad_c4 * exp_term * p1_t
      
      diag_grad_term_m1 <- alpha_prev_grad_m1 * exp_term * p1_t + alpha_prev * dh_dm1
      diag_grad_term_m2 <- alpha_prev_grad_m2 * exp_term * p1_t
      diag_grad_term_m4 <- alpha_prev_grad_m4 * exp_term * p1_t
      
      diag_grad_term_r <- alpha_prev_grad_r * exp_term * p1_t + alpha_prev * dh_dr
      
      shift_grad_term_p1 <- alpha_prev_grad_p1 * exp_term * p4_t
      shift_grad_term_p2 <- alpha_prev_grad_p2 * exp_term * p4_t
      shift_grad_term_p4 <- alpha_prev_grad_p4 * exp_term * p4_t + alpha_prev * dh_dp4
      
      shift_grad_term_c1 <- alpha_prev_grad_c1 * exp_term * p4_t
      shift_grad_term_c2 <- alpha_prev_grad_c2 * exp_term * p4_t
      shift_grad_term_c4 <- alpha_prev_grad_c4 * exp_term * p4_t + alpha_prev * dh_dc4
      
      shift_grad_term_m1 <- alpha_prev_grad_m1 * exp_term * p4_t
      shift_grad_term_m2 <- alpha_prev_grad_m2 * exp_term * p4_t
      shift_grad_term_m4 <- alpha_prev_grad_m4 * exp_term * p4_t + alpha_prev * dh_dm4
      
      shift_grad_term_r <- alpha_prev_grad_r * exp_term * p4_t + alpha_prev * dh_dr
      
      grad_term_p1 <- c(diag_grad_term_p1, shift_grad_term_p1)[valid]
      grad_term_p2 <- c(diag_grad_term_p2, shift_grad_term_p2)[valid]
      grad_term_p4 <- c(diag_grad_term_p4, shift_grad_term_p4)[valid]
      
      grad_term_c1 <- c(diag_grad_term_c1, shift_grad_term_c1)[valid]
      grad_term_c2 <- c(diag_grad_term_c2, shift_grad_term_c2)[valid]
      grad_term_c4 <- c(diag_grad_term_c4, shift_grad_term_c4)[valid]
      
      grad_term_m1 <- c(diag_grad_term_m1, shift_grad_term_m1)[valid]
      grad_term_m2 <- c(diag_grad_term_m2, shift_grad_term_m2)[valid]
      grad_term_m4 <- c(diag_grad_term_m4, shift_grad_term_m4)[valid]
      
      grad_term_r <- c(diag_grad_term_r, shift_grad_term_r)[valid]
    }
    
    else return(safe_fail())
    
    ### ================================
    ### NORMALIZATION AND GRADIENT
    ### ================================
    s <- sum(alpha_new)
    if (!is.finite(s) || s <= 0) return(safe_fail())
    
    total_loglik <- total_loglik + log(s)
    
    # normalized α
    alpha_norm <- alpha_new / s
    
    # gradients of normalized α: α' = (u'/s) - α*(sum(u')/s)
    sum_grad_p1 <- sum(grad_term_p1)
    sum_grad_p2 <- sum(grad_term_p2)
    sum_grad_p4 <- sum(grad_term_p4)
    sum_grad_c1 <- sum(grad_term_c1)
    sum_grad_c2 <- sum(grad_term_c2)
    sum_grad_c4 <- sum(grad_term_c4)
    sum_grad_m1 <- sum(grad_term_m1)
    sum_grad_m2 <- sum(grad_term_m2)
    sum_grad_m4 <- sum(grad_term_m4)
    sum_grad_r  <- sum(grad_term_r)
    
    alpha_prev_grad_p1 <- (grad_term_p1 / s) - alpha_norm * (sum_grad_p1 / s)
    alpha_prev_grad_p2 <- (grad_term_p2 / s) - alpha_norm * (sum_grad_p2 / s)
    alpha_prev_grad_p4 <- (grad_term_p4 / s) - alpha_norm * (sum_grad_p4 / s)
    
    alpha_prev_grad_c1 <- (grad_term_c1 / s) - alpha_norm * (sum_grad_c1 / s)
    alpha_prev_grad_c2 <- (grad_term_c2 / s) - alpha_norm * (sum_grad_c2 / s)
    alpha_prev_grad_c4 <- (grad_term_c4 / s) - alpha_norm * (sum_grad_c4 / s)
    
    alpha_prev_grad_m1 <- (grad_term_m1 / s) - alpha_norm * (sum_grad_m1 / s)
    alpha_prev_grad_m2 <- (grad_term_m2 / s) - alpha_norm * (sum_grad_m2 / s)
    alpha_prev_grad_m4 <- (grad_term_m4 / s) - alpha_norm * (sum_grad_m4 / s)
    
    alpha_prev_grad_r <- (grad_term_r / s) - alpha_norm * (sum_grad_r / s)
    
    # accumulate log-likelihood gradient
    grad_p1 <- grad_p1 + sum_grad_p1 / s
    grad_p2 <- grad_p2 + sum_grad_p2 / s
    grad_p4 <- grad_p4 + sum_grad_p4 / s
    
    grad_c1 <- grad_c1 + sum_grad_c1 / s
    grad_c2 <- grad_c2 + sum_grad_c2 / s
    grad_c4 <- grad_c4 + sum_grad_c4 / s
    
    grad_m1 <- grad_m1 + sum_grad_m1 / s
    grad_m2 <- grad_m2 + sum_grad_m2 / s
    grad_m4 <- grad_m4 + sum_grad_m4 / s
    
    grad_r  <- grad_r  + sum_grad_r  / s
    
    # update α
    alpha_prev <- alpha_norm
  }
  
  grad <- c(grad_p1, grad_p2, grad_p4,
            grad_c1, grad_c2, grad_c4,
            grad_m1, grad_m2, grad_m4,
            grad_r)
  
  return(list(nll = -total_loglik, grad = -grad))
}

fn_opt <- function(par, data) {
  out <- loglik_R(par, data)
  return(out$nll)
}

gr_opt <- function(par, data) {
  out <- loglik_R(par, data)
  return(out$grad)
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

res.de<- read.csv("destep1_start200_025_008_06_065_012_15_020_008_20_02.csv")

res.prm <- data.frame()
start_time <- Sys.time()
for (j in 1:reps) {
  sim.data <- data.frame(t = sim.markov$time.steps[sim.markov$reps == j],
                         Y = sim.markov$ec.steps[sim.markov$reps == j],
                         M = sim.markov$sc.steps[sim.markov$reps == j] +
                           sim.markov$dc.steps[sim.markov$reps == j])
  de.est <- res.de[j,]
  init_par <- c(p1 = as.numeric(de.est["p1"]),
                p2 = as.numeric(de.est["p2"]),
                p4 = as.numeric(de.est["p4"]),
                c1 = as.numeric(de.est["c1"]),
                c2 = as.numeric(de.est["c2"]),
                c4 = as.numeric(de.est["c4"]),
                m1 = as.numeric(de.est["m1"]),
                m2 = as.numeric(de.est["m2"]),
                m4 = as.numeric(de.est["m4"]),
                r = as.numeric(de.est["r"])
  )
  
  lower <- c(p1 = 1e-5, p2 = 1e-5, p4 = 1e-5,
             c1 = 1e-5, c2 = 1e-5, c4 = 1e-5,
             m1 = 1e-5, m2 = 1e-5, m4 = 1e-5,
             r = 1e-5)
  upper <- c(p1 = 0.99999, p2 = 0.99999, p4 = 0.99999,
             c1 = 100, c2 = 100, c4 = 100,
             m1 = 100, m2 = 100, m4 = 100,
             r = 0.99999)
  
  fit <- optim(
    par = init_par,
    fn  = fn_opt,
    gr  = gr_opt,
    data = sim.data,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      maxit = 2000,
      factr = 1e7,
      trace = 1
    )
  )
  
  par.est <- fit$par
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

# 
# res.prm %>% 
#   pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
#   ggplot(aes(x = "", y = Value)) +
#   geom_violin(aes(fill = Variable), color = "black") +
#   geom_boxplot(fill = NA, width = 0.5)+
#   facet_wrap(~ Variable, scale = "free_y", ncol = 3) +
#   geom_jitter(width = 0.1, alpha = 0.3) + 
#   theme_minimal() +
#   labs(title = "Violin Plot for Each Parameter Estimates", x = "",
#        fill= "Parameters", y = "Estimates")
