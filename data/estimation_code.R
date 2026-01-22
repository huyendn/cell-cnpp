library(tidyr)
library(dplyr)
library(ggplot2)
library(DEoptim)
library(patchwork)

### functions for estimation
log_likelihood_noduds_de <- function(params, data) {
  params <- as.numeric(params)
  names(params) <- c("u", "v", "c1", "c2", "m1", "m2")
  
  # Reparameterization
  u <- params["u"]
  v <- params["v"]
  
  p1 <- u
  p2 <- (1 - u) * v
  
  c1 <- params["c1"]; c2 <- params["c2"]
  m1 <- params["m1"]; m2 <- params["m2"]
  
  x <- as.numeric(data$X)
  t <- as.numeric(data$t)
  
  A <- p1 / (1 + c1 * (t - m1)^2)
  B <- p2 / (1 + c2 * (t - m2)^2)
  
  epsilon <- 1e-10
  C <- pmax(1 - A - B, epsilon)
  
  ll <- sum(
    (x ==  1) * log(A + epsilon) +
      (x ==  0) * log(B + epsilon) +
      (x == -1) * log(C)
  )
  
  return(-ll)
}

r_estimator <- function(data){
  Si <- data$Si
  deltat <- data$deltat
  n <- nrow(data)
  r <- n/(sum(Si*deltat))
  return(r)
}

### functions for prediction
int.p.func <- function(t, p, c, m) {
  (p / sqrt(c)) * (atan(sqrt(c) * (t - m)) + atan(sqrt(c) * m))
}

make_expected_noduds_data <- function(time.input,
                                      p1.true, c1.true, m1.true,
                                      p2.true, c2.true, m2.true,
                                      # p4.true, c4.true, m4.true,
                                      s0.true, r.true) {
  p1.func <- function(t){
    p1.prob <- p1.true/(1+c1.true*(t-m1.true)^2)
    return(p1.prob)
  }
  
  p2.func <- function(t){
    p2.prob <- p2.true/(1+c2.true*(t-m2.true)^2)
    return(p2.prob)
  }
  
  p3.func <- function(t){
    # p3.prob <- 1- p1.true/(1+c1.true*(t-m1.true)^2)-p2.true/(1+c2.true*(t-m2.true)^2) - p4.true/(1+c4.true*(t-m4.true)^2)
    p3.prob <-1- p1.true/(1+c1.true*(t-m1.true)^2)-p2.true/(1+c2.true*(t-m2.true)^2)
    return(p3.prob)
  }
  
  int.p1 <- function(t) int.p.func(t, p1.true, c1.true, m1.true)
  int.p2 <- function(t) int.p.func(t, p2.true, c2.true, m2.true)
  # int.p4 <- function(t) int.p.func(t, p4.true, c4.true, m4.true)
  
  int.p3 <- function(t) {
    # t - int.p1(t) - int.p2(t) - int.p4(t)
    t - int.p1(t) - int.p2(t)
  }
  
  sc.mean <- function(t) {
    s0.true * exp(r.true * (int.p1(t) - int.p3(t)))
  }
  
  sc.var.bp <- function(t) {
    A <- exp(2 * r.true * (int.p1(t) - int.p3(t)))

    integrand_B <- function(x) {
      (p1.func(x) + p3.func(x)) *
        exp(-r.true * (int.p1(x) - int.p3(x)))
    }

    B <- integrate(integrand_B, lower = 0, upper = t)$value
    s0.true * r.true * A * B
  }
  ec.mean <- function(t){
    integrand <- function(x){
      (p2.func(x) + 2 * p3.func(x))* exp(r.true * (int.p1(x) - int.p3(x)))
    }
    out <- integrate(integrand, lower = 0, upper = t)$value
    return(out*s0.true*r.true)
  }
  
  data.frame(
    time = time.input,
    sc.mean = sapply(time.input, sc.mean),
    sc.var.bp = sapply(time.input, sc.var.bp),
    ec.mean = sapply(time.input, ec.mean)
  )
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

### for plotting
colors1 <- c("observed MkP+ErP" = "lightpink", "observed MEP" = "skyblue",
             "estimated MkP+ErP" = "indianred", "estimated MEP" = "darkblue")


###
###
###
mep_tpo_count <- read.csv("mep_tpo_count.csv")
mep_tpo_count <- mep_tpo_count[-c(1:4),]

mep_tpo_count$time_shift <- mep_tpo_count$time - min(mep_tpo_count$time)

### TPO
mep_tpo_fitdata <- data.frame(X = diff(mep_tpo_count$mep),
                              t = mep_tpo_count$time_shift[-1],
                              deltat = diff(mep_tpo_count$time_shift),
                              Si = mep_tpo_count$mep
                              [-length(mep_tpo_count$mep)])
head(mep_tpo_fitdata)
t_min <- min(mep_tpo_fitdata$t, 0)
t_max <- max(mep_tpo_fitdata$t, 100)
lower <- c(u = 1e-5, v = 1e-5, 
           c1 = 1e-5, c2 = 1e-5, 
           m1 = t_min, m2 = t_min)
upper <- c(u = 0.99999, v = 0.99999,
           c1 = 100, c2 = 100, 
           m1 = t_max, m2 = t_max)
set.seed(123)
DE_res <- DEoptim(
  fn = log_likelihood_noduds_de,
  lower = lower,
  upper = upper,
  data = mep_tpo_fitdata,
  DEoptim.control(NP = 150, itermax = 600, F = 0.9, CR = 0.8, strategy = 6,
                  trace = FALSE)
)
u_hat <- DE_res$optim$bestmem["u"]
v_hat <- DE_res$optim$bestmem["v"]

p1_hat <- u_hat
p2_hat <- (1 - u_hat) * v_hat

p1_hat; p2_hat

par.est <- DE_res$optim$bestmem
r.est <- r_estimator(mep_tpo_fitdata)
tpo.noduds.est <- c(p1_hat, p2_hat, par.est[3:6], r.est)
names(tpo.noduds.est) <- c("p1", "p2", "c1", "c2", "m1", "m2", "r")
tpo.noduds.est
tpo.noduds.est["p1"]; tpo.noduds.est["c1"]; tpo.noduds.est["m1"]
# 0.4560184 0.001138745 19.0687 
tpo.noduds.est["p2"]; tpo.noduds.est["c2"]; tpo.noduds.est["m2"]
# 0.5439761  0.0001757318  21.74811
tpo.noduds.est["r"]
# 0.04528165

mep_tpo_pred <- make_expected_noduds_data(time.input = seq(min(mep_tpo_count$time_shift), max(mep_tpo_count$time_shift+10), by = 0.1),
                                          p1.true = tpo.noduds.est["p1"],
                                          c1.true = tpo.noduds.est["c1"],
                                          m1.true = tpo.noduds.est["m1"],
                                          p2.true = tpo.noduds.est["p2"],
                                          c2.true = tpo.noduds.est["c2"],
                                          m2.true = tpo.noduds.est["m2"],
                                          s0.true = mep_tpo_count$mep[1],
                                          r.true = tpo.noduds.est["r"])

tpo_est_plot <- ggplot() +
  geom_point(data = mep_tpo_count, 
             aes(x = time_shift + min(time), y = daughter, col = "observed MkP+ErP"),
             pch = 19, size = 2) +
  geom_point(data = mep_tpo_count, 
             aes(x = time_shift+ min(time), y = mep, col = "observed MEP"),
             pch = 19, size = 2) +
  geom_line(data = mep_tpo_pred, 
            aes(x = time + min(mep_tpo_count$time), y = ec.mean, col = "estimated MkP+ErP"),
            linetype = 4, linewidth = 1) +
  geom_line(data = mep_tpo_pred, 
            aes(x = time + min(mep_tpo_count$time), y = sc.mean, col = "estimated MEP"),
            linetype = 4, linewidth = 1) +
  scale_y_continuous(limits = c(-10, 250)) +
  scale_x_continuous(limits = c(20, 200)) +
  labs(x = "Time (in hours)", y = "Cell counts", color = "",
       title = "Lacking TPO condition",
       # subtitle = "Observed and Estimated Cell Counts"
  )+
  scale_color_manual(values = colors1) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.text.align = 0,
    legend.position = "bottom"
  )

##
##
mep_ctrl_count <- read.csv("mep_ctrl_count.csv")
mep_ctrl_count <- mep_ctrl_count[-c(1:14),]

mep_ctrl_count$time_shift <- mep_ctrl_count$time - min(mep_ctrl_count$time)

mep_ctrl_fitdata <- data.frame(X = diff(mep_ctrl_count$mep),
                               t = mep_ctrl_count$time_shift[-1],
                               deltat = diff(mep_ctrl_count$time_shift),
                               Si = mep_ctrl_count$mep
                               [-length(mep_ctrl_count$mep)])

head(mep_ctrl_fitdata)
t_min <- min(mep_ctrl_fitdata$t, 0)
t_max <- max(mep_ctrl_fitdata$t, 100)
lower <- c(u = 1e-5, v = 1e-5, 
           c1 = 1e-5, c2 = 1e-5, 
           m1 = t_min, m2 = t_min)
upper <- c(u = 0.99999, v = 0.99999,
           c1 = 1, c2 = 1, 
           m1 = t_max, m2 = t_max)
set.seed(124)
DE_res <- DEoptim(
  fn = log_likelihood_noduds_de,
  lower = lower,
  upper = upper,
  data = mep_ctrl_fitdata,
  DEoptim.control(NP = 150, itermax = 600, F = 0.9, CR = 0.8, strategy = 6,
                  trace = FALSE
  )
)

u_hat <- DE_res$optim$bestmem["u"]
v_hat <- DE_res$optim$bestmem["v"]

p1_hat <- u_hat
p2_hat <- (1 - u_hat) * v_hat

p1_hat; p2_hat

par.est <- DE_res$optim$bestmem
r.est <- r_estimator(mep_ctrl_fitdata)
par.est
r.est


ctrl.est <- c(p1_hat, p2_hat, par.est[3:6], r.est)
names(ctrl.est) <- c("p1", "p2", "c1", "c2", "m1", "m2", "r")

mep_ctrl_pred <- make_expected_noduds_data(time.input = seq(min(mep_ctrl_count$time_shift), max(mep_ctrl_count$time_shift+10), by = 0.1),
                                          p1.true = ctrl.est["p1"],
                                          c1.true = ctrl.est["c1"],
                                          m1.true = ctrl.est["m1"],
                                          p2.true = ctrl.est["p2"],
                                          c2.true = ctrl.est["c2"],
                                          m2.true = ctrl.est["m2"],
                                          s0.true = mep_ctrl_count$mep[1],
                                          r.true = ctrl.est["r"])

ctrl_est_plot <- ggplot() +
  geom_point(data = mep_ctrl_count, 
             aes(x = time_shift + min(time), y = daughter, col = "observed MkP+ErP"),
             pch = 19, size = 2) +
  geom_point(data = mep_ctrl_count, 
             aes(x = time_shift+ min(time), y = mep, col = "observed MEP"),
             pch = 19, size = 2) +
  geom_line(data = mep_ctrl_pred, 
            aes(x = time + min(mep_ctrl_count$time), y = ec.mean + mep_ctrl_count$daughter[1], col = "estimated MkP+ErP"),
            linetype = 4, linewidth = 1) +
  geom_line(data = mep_ctrl_pred, 
            aes(x = time + min(mep_ctrl_count$time), y = sc.mean, col = "estimated MEP"),
            linetype = 4, linewidth = 1) +
  scale_y_continuous(limits = c(-10, 250)) +
  scale_x_continuous(limits = c(20, 200)) +
  labs(x = "Time (in hours)", y = "Cell counts", color = "",
       title = "Control condition",
       # subtitle = "Observed and Expected Cell Counts"
  )+
  scale_color_manual(values = colors1) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.text.align = 0,
    legend.position = "bottom"
  )

(ctrl_est_plot | tpo_est_plot) +  plot_annotation(
  title = "Observed and Estimated Cell Counts"
) + plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )
