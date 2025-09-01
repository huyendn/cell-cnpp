library(nhppp)
library(ggplot2)
library(dplyr)

# proliferation function q_t = p2_t + 2p3_t
q_t1 <- function(tms, prm) {
  return(2 - 1/(prm[1] + prm[2]*tms + prm[3]*tms^2))
}

# integrate proliferation function
int_q <- function(t, prm) {
  2*t - (2*atan((prm[2] + 2*prm[3]*t)/sqrt(4*prm[1]*prm[3] - prm[2]^2))/
    sqrt(4*prm[1]*prm[3] - prm[2]^2))
}

## model parameters
s0 <- 200
r.test <- 0.2
prm <-c(1.2, -0.106, 0.005)

# expected differentiated cell function
ec_func <- function(t){
  integrand <- function(t){
    out <- q_t1(tms = t, prm = prm)*intensity_func(t)
  }
  out.ec <- integrate(integrand, lower = 0, upper = t)$value
  return(out.ec)
}

# expected stem cell function
sc_func <- function(t){
  out <- s0*exp(r.test*(t-int_q(t,prm)+int_q(0,prm)))
  return(out)
}

# variance differentiated cell function

# variance stem cell function


# intensity function 
intensity_func <- function(t){
  out <- r.test*s0*exp(r.test*(t-int_q(t,prm)+int_q(0,prm)))
}

# slope and intercept of the log-linear for thinning algorithm
# need to check if this function is always greater than the intensity function
slope <- -0.01
intercept <- 6

# simulation function, 
cnpp_sim <- function(s0, prm, slope = slope, intercept = slope,
                     min_q_t = 0.0001, max_it = 50000, k = 3, seed0 = 0,
                     t_max = 100) {
  set.seed(seed0)
  # draw event time from the nhppp package using thinning algorithm
  time.draw <- draw_intensity(lambda = intensity_func,
                              line_majorizer_intercept = intercept,
                              line_majorizer_slope = slope,
                              line_majorizer_is_loglinear = TRUE,
                              t_min = 0, t_max = t_max)
  sc.current <- s0
  ec.current <- dc.current <- time.current <- 0
  sc.steps <- c(s0)
  ec.steps <- c(ec.current)
  # dc.steps <- c(dc.current)
  time.steps <- c(time.current)
  itnum <- 0
  index <- 0
  while (sc.current > 0 & index < length(time.draw)) {
    itnum <- itnum + 1
    index <- index + 1
    time.current <- time.draw[index]
    # calculate proliferation function at event time
    q_t <- q_t1(time.current, prm)
    if ((q_t < min_q_t| (itnum > max_it))) {
      cat(q_t, "\n")
      break
    }
    # determine the probability of each type of division
    probs <- c((2-q_t)/k, max(0, (k-2)*(2-q_t)/k), max(0,(2-k + (k-1)*q_t)/k))
    if (sum(probs) > 1) {
      probs <- probs/sum(probs)
    }
    # sample the type of division
    event_type <- sample(3, 1, prob=probs)
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
    reps.steps <- rep(i, length(sc.steps))
  }
  sim.dat <- data.frame(time = time.steps, sc = sc.steps, ec = ec.steps)
  return(sim.dat)
}

reps <- 30
sim <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric())
for (i in 1:reps) {
  test.sim <- cnpp_sim(s0 = s0, prm, seed0 = seed0+i)
  sim <- rbind(sim,cbind(reps = rep(i, nrow(test.sim)), test.sim))
}

tms <- seq(0, max(sim$time), by = 5)

dat <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric(),
                  time.collected = numeric())

for (i in 1:reps) {
  subset <- sim %>% filter(reps == i)
  for (j in 2:length(tms)) {
    time.before <- tms[j-1]
    time.collected <- tms[j]
    if (sum(subset$time <= time.collected  & subset$time > time.before) == 0){
      dat <- dat
    } else{
      ind <- max(which(subset$time <= time.collected  & subset$time > time.before))
      dat <- rbind(dat, cbind(subset[ind,], time.collected))
    }
    
  }
}


ggplot() +
  geom_point(data = dat,aes(x = time, y = ec, col = "differentiated cell"), alpha = 0.5)+
  geom_point(data = dat,aes(x = time, y = sc, col = "stem cell"), alpha = 0.5)+
  labs(x = "time (in hours)", y = "cell counts", color = "",
       title = "Simulated cell counts")+
  scale_color_manual(values = c("differentiated cell" = "skyblue3",
                                "stem cell" = "indianred3")) +
  theme_classic() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14), 
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13))

sim %>% 
ggplot() +
  geom_point(aes(x = time, y = ec, col = "differentiated cell"), alpha = 0.5)+
  geom_point(aes(x = time, y = sc, col = "stem cell"), alpha = 0.5)+
  labs(x = "time (in hours)", y = "cell counts", color = "",
       title = "Simulated cell counts")+
  scale_color_manual(values = c("differentiated cell" = "skyblue3",
                                "stem cell" = "indianred3")) +
  theme_classic() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14), 
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13))

exp.mean.sc <- exp.mean.ec <- c()
for (i in 2:length(tms)) {
  exp.mean.sc <- c(exp.mean.sc,sc_func(tms[i]))
  exp.mean.ec <- c(exp.mean.ec,ec_func(tms[i]))
}

knitr::kable(dat %>% group_by(time.collected) %>%
  summarise(n = n(), samp.mean.sc = mean(sc), samp.mean.ec = mean(ec)) %>%
  mutate(theo.mean.sc = exp.mean.sc,
         theo.mean.ec = exp.mean.ec))
