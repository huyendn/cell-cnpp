library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(patchwork)
library(reshape2)

library(statmod)
library(fitdistrplus)
library(goftest)

# probability plot
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

configuration <- list(
  set1 = c("p1" = 0.55, "p2" = 0.15, "p4" = 0.2,
           "c1" = 0.005, "c2" = 0.012, "c4" = 0.008,
           "m1" = 4, "m2" = 12, "m4" = 20),
  set2 = c("p1" = 0.3, "p2" = 0.3, "p4" = 0.3,
           "c1" = 0.008, "c2" = 0.008, "c4" = 0.008,
           "m1" = 6, "m2" = 15, "m4" = 25),
  set3 = c("p1" = 0.25, "p2" = 0.65, "p4" = 0.1,
           "c1" = 0.005, "c2" = 0.012, "c4" = 0.008,
           "m1" = 6, "m2" = 15, "m4" = 20)
)

make_plot <- function(index, configuration) {
  
  par <- configuration[[index]]
  
  p1.true <- par["p1"];  c1.true <- par["c1"];  m1.true <- par["m1"]
  p2.true <- par["p2"];  c2.true <- par["c2"];  m2.true <- par["m2"]
  p4.true <- par["p4"];  c4.true <- par["c4"];  m4.true <- par["m4"]
  
  p1_fun <- function(t) p1.true / (1 + c1.true * (t - m1.true)^2)
  p2_fun <- function(t) p2.true / (1 + c2.true * (t - m2.true)^2)
  p4_fun <- function(t) p4.true / (1 + c4.true * (t - m4.true)^2)
  p3_fun <- function(t) 1 - p1_fun(t) - p2_fun(t) - p4_fun(t)
  
  t_vals <- seq(0, 60, length.out = 500)
  
  data <- data.frame(
    t = t_vals,
    p1 = p1_fun(t_vals),
    p2 = p2_fun(t_vals),
    p4 = p4_fun(t_vals),
    p3 = p3_fun(t_vals)
  )
  
  data_long <- pivot_longer(
    data, 
    cols = c(p1, p2, p4, p3), 
    names_to = "Function", 
    values_to = "Value"
  )
  
  ggplot(data_long, aes(x = t, y = Value, color = Function)) +
    geom_line(linewidth = 1.2) +
    geom_hline(aes(yintercept = 1)) +
    scale_color_manual(
      values = c("skyblue3", "indianred3", "darkgreen", "orange3"),
      labels = c(
        expression(p[1](t)),
        expression(p[2](t)),
        expression(p[3](t)),
        expression(p[4](t))
      )
    ) +
    labs(
      title = paste0("Configuration ", index),
      x = "Time",
      y = "Probability",
      color = "Function"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      legend.text.align = 0
    )
}

# Generate the three plots
p1 <- make_plot(1, configuration)
p2 <- make_plot(2, configuration)
p3 <- make_plot(3, configuration)

# Patchwork combination
(p1 | p2 | p3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

###
###
###
### prediction plot
compute_ratio_for_set <- function(sim.markov, res.prm, param_set_id) {
  
  time.seq <- seq(5, 40, by = 5)
  
  sim.aggregate <- sim.markov %>%
    group_by(reps) %>%
    arrange(time.steps) %>%
    summarise(
      tibble(
        target.time = time.seq,
        idx = map_dbl(time.seq, ~max(which(time.steps <= .x))),
        time.steps = time.steps[idx],
        sc.steps  = sc.steps[idx],
        ec.steps  = ec.steps[idx],
        dc.steps  = dc.steps[idx]
      )
    ) %>% 
    unnest(cols = c(target.time, idx, time.steps, sc.steps, ec.steps, dc.steps))
  
  reps <- 100
  s0.true <- 200
  
  pred <- data.frame()
  for (i in 1:reps) {
    sim.aggregate.subset <- sim.aggregate %>% filter(reps == i)
    time.input <- sim.aggregate.subset$time.steps
    
    r.true <- res.prm[i, "r"]
    p1.true <- res.prm[i, "p1"]; c1.true <- res.prm[i, "c1"]; m1.true <- res.prm[i, "m1"]
    p2.true <- res.prm[i, "p2"]; c2.true <- res.prm[i, "c2"]; m2.true <- res.prm[i, "m2"]
    p4.true <- res.prm[i, "p4"]; c4.true <- res.prm[i, "c4"]; m4.true <- res.prm[i, "m4"]
    
    p1.func <- function(t) p1.true/(1+c1.true*(t-m1.true)^2)
    p2.func <- function(t) p2.true/(1+c2.true*(t-m2.true)^2)
    p4.func <- function(t) p4.true/(1+c4.true*(t-m4.true)^2)
    p3.func <- function(t) 1 - p1.func(t) - p2.func(t) - p4.func(t)
    
    int.p1.func <- function(t) (p1.true/sqrt(c1.true))*(atan(sqrt(c1.true)*(t-m1.true))+atan(sqrt(c1.true)*m1.true))
    int.p2.func <- function(t) (p2.true/sqrt(c2.true))*(atan(sqrt(c2.true)*(t-m2.true))+atan(sqrt(c2.true)*m2.true))
    int.p4.func <- function(t) (p4.true/sqrt(c4.true))*(atan(sqrt(c4.true)*(t-m4.true))+atan(sqrt(c4.true)*m4.true))
    int.p3.func <- function(t) t - int.p1.func(t) - int.p2.func(t) - int.p4.func(t)
    
    sc.mean <- function(t) s0.true*exp(r.true*(int.p1.func(t)-int.p3.func(t)))
    
    ec.mean <- function(t) {
      integrate(function(x) (p2.func(x)+2*p3.func(x))*r.true*sc.mean(x), 0, t)$value
    }
    
    dc.mean <- function(t) {
      integrate(function(x) p4.func(x)*r.true*sc.mean(x), 0, t)$value
    }
    
    pred.subset <- data.frame(
      reps = sim.aggregate.subset$reps,
      target.time = sim.aggregate.subset$target.time,
      time.steps = sim.aggregate.subset$time.steps,
      sc.steps = sapply(time.input, sc.mean),
      ec.steps = sapply(time.input, ec.mean),
      dc.steps = sapply(time.input, dc.mean)
    )
    pred <- rbind(pred, pred.subset)
  }
  
  pred$m.steps <- pred$sc.steps + pred$dc.steps
  sim.aggregate$m.steps <- sim.aggregate$sc.steps + sim.aggregate$dc.steps
  
  ratio <- data.frame(
    reps = sim.aggregate$reps,
    target.time = sim.aggregate$target.time,
    time.steps = sim.aggregate$time.steps,
    SC = (sim.aggregate$sc.steps+0.5)/(pred$sc.steps+0.5),
    FC = (sim.aggregate$ec.steps+0.5)/(pred$ec.steps+0.5),
    DC = (sim.aggregate$dc.steps+0.5)/(pred$dc.steps+0.5),
    SCandDC = (sim.aggregate$m.steps+0.5)/(pred$m.steps+0.5),
    param_set = param_set_id
  )
  
  return(ratio)
}

###
sim.markov.1 <- read.csv("est_res/sim_start200_055_005_04_015_012_12_020_008_20_02.csv")
res.prm.1 <- read.csv("est_res/est_start200_055_005_04_015_012_12_020_008_20_02.csv")

sim.markov.2 <- read.csv("est_res/sim_start200_030_008_06_030_008_15_030_008_25_02.csv")
res.prm.2 <- read.csv("est_res/est_start200_030_008_06_030_008_15_030_008_25_02.csv")

sim.markov.3 <- read.csv("est_res/sim_start200_025_005_06_065_012_15_010_008_20_02.csv")
res.prm.3 <- read.csv("est_res/est_start200_025_005_06_065_012_15_010_008_20_02.csv")

ratio1 <- compute_ratio_for_set(sim.markov.1, res.prm.1, 1)
ratio2 <- compute_ratio_for_set(sim.markov.2, res.prm.2, 2)
ratio3 <- compute_ratio_for_set(sim.markov.3, res.prm.3, 3)

ratio_all <- bind_rows(ratio1, ratio2, ratio3)

ratio_long <- ratio_all %>%
  pivot_longer(
    cols = c(SC, FC, DC, SCandDC),
    names_to = "type",
    values_to = "value"
  ) %>% 
  mutate(
    type = factor(type, levels = c("SC", "DC", "SCandDC", "FC")),
    param_set = factor(param_set)
  )
facet_labels <- c(
  "SC" = "Viable Stem Cells",
  "DC" = "Nonviable Stem Cells",
  "FC" = "Differentiated Cells",
  "SCandDC" = "Total Stem Cells"
)

ratio_long %>%
  filter(target.time > 15 & target.time < 40) %>%
  ggplot(aes(x = factor(target.time), y = value, fill = param_set)) +
  
  geom_boxplot(outlier.size = 1, position = position_dodge(width = 0.8)) +
  geom_hline(aes(yintercept = 1)) +
  
  scale_fill_manual(values = c("skyblue3", "indianred3", "darkgreen")) +
  labs(
    y = "Observed / Predicted Ratios",
    x = "Time",
    fill = "Parameter Configuration Set"
  ) +
  facet_wrap(~ type, labeller = labeller(type = facet_labels)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    # legend.title = element_text(face = "bold"),
    legend.text.align = 0,
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "bottom"
  )

###
###
### variance plot

par <- configuration[[1]]
par
s0.true <- 200
r.true <- 0.2
p1.true <- par["p1"]; p2.true <- par["p2"]; p4.true <- par["p4"]
c1.true <- par["c1"]; c2.true <- par["c2"]; c4.true <- par["c4"]
m1.true <- par["m1"]; m2.true <- par["m2"]; m4.true <- par["m4"]


p1.func <- function(t){
  p1.prob <- p1.true/(1+c1.true*(t-m1.true)^2)
  return(p1.prob)
}

p2.func <- function(t){
  p2.prob <- p2.true/(1+c2.true*(t-m2.true)^2)
  return(p2.prob)
}

p3.func <- function(t){
  p3.prob <- 1- p1.true/(1+c1.true*(t-m1.true)^2)-p2.true/(1+c2.true*(t-m2.true)^2) - p4.true/(1+c4.true*(t-m4.true)^2)
  return(p3.prob)
}

p4.func <- function(t){
  p4.prob <- p4.true/(1+c4.true*(t-m4.true)^2)
  return(p4.prob)
}

int.p1.func <- function(t){
  (p1.true/sqrt(c1.true))*(atan(sqrt(c1.true)*(t-m1.true))+atan(sqrt(c1.true)*m1.true))
}

int.p2.func <- function(t){
  (p2.true/sqrt(c2.true))*(atan(sqrt(c2.true)*(t-m2.true))+atan(sqrt(c2.true)*m2.true))
}

int.p4.func <- function(t){
  (p4.true/sqrt(c4.true))*(atan(sqrt(c4.true)*(t-m4.true))+atan(sqrt(c4.true)*m4.true))
}

int.p3.func <- function(t){
  t-int.p1.func(t) - int.p2.func(t) - int.p4.func(t)
}

sc.mean <- function(t){
  s0.true*exp(r.true*(int.p1.func(t)-int.p3.func(t)))
}

sc.var.bp <- function(t) {
  A <- exp(2*r.true*(int.p1.func(t)-int.p3.func(t)))
  integrand_B <- function(x){
    (p1.func(x) + p3.func(x))*exp(-r.true*(int.p1.func(x)-int.p3.func(x)))
  }
  B <- integrate(integrand_B, lower = 0, upper = t)$value
  return(s0.true*r.true*A * B)
}

sc.autocorr.bp <- function(t, u) {
  if (t >= u) {
    return(sc.mean(t) / sc.mean(u) * sc.var.bp(u))
  } else {
    return(sc.mean(u) / sc.mean(t) * sc.var.bp(t))
  }
}

time.input <- seq(0, 60, by = 0.1)

mean.output <- sapply(time.input, sc.mean)
var.bp.output <- sapply(time.input, sc.var.bp)

theoretical.dat <- data.frame(time = time.input,
                              mean = mean.output,
                              var.bp = var.bp.output)

sim.markov.1 %>% filter(reps <= 50) %>% filter(time.steps <= 60) %>%
  ggplot() +
  geom_point(aes(x = time.steps, y = sc.steps), size = 0.02, alpha = 0.05)+
  geom_line(data = theoretical.dat, aes(x = time, y = mean), 
            col = "skyblue3", linewidth = 1.5) +
  geom_ribbon(data = theoretical.dat, 
              aes(x = time, ymax = mean + 2*sqrt(var.bp), ymin = mean-2*sqrt(var.bp)),
              fill = "skyblue3", alpha = 0.5) +
  # geom_ribbon(data = theoretical.dat, 
  #             aes(x = time, ymax = mean + 2*sqrt(var.cpnn), ymin = mean-2*sqrt(var.cpnn)),
  #             fill = "blue", alpha = 0.2) +
  labs(x = "Time", y = "Cell counts", title = "Theoretical Mean and Variance of Stem Cells") +
  theme_minimal()
        

mean_var_plot <- sim.markov.1 %>% 
  filter(reps <= 50, time.steps <= 60) %>%
  ggplot() +
  
  # Empirical points
  geom_point(
    aes(x = time.steps, y = sc.steps),
    size = 0.02, alpha = 0.05
  ) +
  
  # Theoretical mean line
  geom_line(
    data = theoretical.dat,
    aes(x = time, y = mean),
    color = "skyblue3", linewidth = 1.5
  ) +
  
  # Theoretical 2 SD confidence band
  geom_ribbon(
    data = theoretical.dat,
    aes(x = time,
        ymin = mean - 2 * sqrt(var.bp),
        ymax = mean + 2 * sqrt(var.bp)),
    fill = "skyblue3", alpha = 0.3
  ) +
  
  # Labels
  labs(
    x = "Time",
    y = "Cell counts",
    title = "Mean and Variance of Stem Cell Counts"
  ) +
  
  # Minimal theme but customized to match boxplot formatting
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # centered bold title
    axis.title = element_text(face = "bold"),               # bold axis titles
    strip.text = element_text(face = "bold", size = 16),    # larger + bold facet labels (if any)
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # box outline
    legend.position = "bottom"
  )

### autocovariance
tms.cor <- seq(5,25, by = 5)

# theoretical
# Create a matrix of E[X(t) X(u)]
theo_cov <- outer(tms.cor, tms.cor, Vectorize(sc.autocorr.bp))

colnames(theo_cov) <- tms.cor
rownames(theo_cov) <- tms.cor

theosd_vec <- sqrt(diag(theo_cov))
theo_corr <- theo_cov / (theosd_vec %o% theosd_vec)

print(round(theo_cov, 4))
print(round(theo_corr, 4))

sim.aggregate <- sim.markov.1 %>% filter(reps <= 50) %>%
  group_by(reps) %>%
  arrange(time.steps) %>%
  summarise(
    tibble(
      target.time = tms.cor,
      idx = map_dbl(tms.cor, ~max(which(time.steps <= .x))),
      time.steps = time.steps[idx],
      sc.steps  = sc.steps[idx],
      ec.steps  = ec.steps[idx],
      dc.steps  = dc.steps[idx]
    )
  ) %>% 
  unnest(cols = c(target.time, idx, time.steps, sc.steps, ec.steps, dc.steps))

head(sim.aggregate)

dat_wide <- sim.aggregate %>% select(reps, target.time, sc.steps)%>%
  pivot_wider(names_from = target.time, values_from = sc.steps)

X <- as.matrix(dat_wide[,-1])  # remove 'reps' column
R <- nrow(X)


# Center data
X_centered <- sweep(X, 2, colMeans(X))
empirical_cov <- (t(X_centered) %*% X_centered) / R

sd_vec <- sqrt(diag(empirical_cov))
empirical_corr <- empirical_cov / (sd_vec %o% sd_vec)

print(round(empirical_corr, 4)); print(round(theo_corr, 4))
print(round(empirical_cov, 4)); print(round(theo_cov, 4))


# Convert matrices to long format
empirical_long <- melt(empirical_corr)
theoretical_long <- melt(theo_corr)

# Empirical autocorrelation heatmap with skyblue3 gradient
emp_plot <- ggplot(empirical_long, aes(x = Time1, y = Time2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "skyblue3"  
  ) +
  labs(title = "Empirical Autocorrelation Heatmap",
       x = "Time",
       y = "Time",
       fill = "Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Theoretical autocorrelation heatmap with skyblue3 gradient
theo_plot <-ggplot(theoretical_long, aes(x = Time1, y = Time2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "skyblue3"
  ) +
  labs(title = "Theoretical Autocorrelation Heatmap",
       x = "Time",
       y = "Time",
       fill = "Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

theo_plot + emp_plot

### 
###
### stopping time
stop_time_plot <- sim.markov.1 %>% filter(sc.steps == 0) %>%
  ggplot(aes(x = time.steps)) +
  geom_histogram(
    aes(y = after_stat(count)/sum(after_stat(count))),
    binwidth = 1.5, 
    color = "black",fill = "skyblue3") +
  geom_vline(aes(xintercept = 1/0.2*log(200) - digamma(1)), linetype = "dashed") +
  scale_x_continuous(limits = c(25, 115))+
  # scale_y_continuous(limits = c(0, 0.15)) +
  labs(x = "Time", y = "Relative frequency",
       title = "Stopping Time"
       ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # centered bold title
    axis.title = element_text(face = "bold"),               # bold axis titles
    strip.text = element_text(face = "bold", size = 16),    # larger + bold facet labels (if any)
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # box outline
    legend.position = "bottom"
  )


stop.time <- sim.markov.1 %>% filter(sc.steps == 0) %>% pull(time.steps)
min.stop.time <- 1/0.2*log(200) - digamma(1)
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


