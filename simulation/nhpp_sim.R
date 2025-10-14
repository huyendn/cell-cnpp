library(nhppp)
library(ggplot2)
library(dplyr)
library(plotly)
# library(car)
# library(rgl)

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
s0.test <- 200
r.test <- 0.2
prm.test <-c(1.2, -0.106, 0.005)

# intensity function 
intensity_func <- function(t){
  out <- r.test*s0.test*exp(r.test*(t-int_q(t,prm.test)+int_q(0,prm.test)))
}

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

# slope and intercept of the log-linear for thinning algorithm
# need to check if this function is always greater than the intensity function
slope <- -0.008
intercept <- 5

# slope <- -0.008
# intercept <- 3.5
### check intensity function is always smaller than the majorization function
majorization_func <- function(t, a, b){
  out <- exp(a + b*t)
  return(out)
}
input.time <- seq(0, 100)
output.intensity <- output.majorization <- c()
for (i in 1:length(input.time)) {
  output.intensity <- c(output.intensity, intensity_func(input.time[i]))
  output.majorization <- c(output.majorization, 
                           majorization_func(input.time[i], a = intercept, b = slope))
}
dat.majorization <- data.frame(time = input.time,
                               intensity = output.intensity,
                               majorization = output.majorization)
ggplot(data = dat.majorization) +
  geom_line(aes(x = time, y = intensity), color = "indianred")+
  geom_line(aes(x = time, y = majorization), color = "skyblue3")

# simulations settings
min_q_t = 0.0001
max_it = 50000
k = 3
seed0 = 922125
t_max = 100
reps <- 100



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
sim <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric())
for (i in 1:reps) {
  # draw event time from the nhppp package using thinning algorithm
  time.draw <- draw_intensity(lambda = intensity_func,
                              line_majorizer_intercept = intercept,
                              line_majorizer_slope = slope,
                              line_majorizer_is_loglinear = TRUE,
                              t_min = 0, t_max = t_max)
  sc.current <- s0.test
  ec.current <- dc.current <- time.current <- 0
  sc.steps <- c(s0.test)
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
    q_t <- q_t1(time.current, prm.test)
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
  sim <- rbind(sim,cbind(reps = rep(i, length(sc.steps)), time.steps, sc.steps, ec.steps))
}
# View(sim)

tms <- seq(0, 80, by = 5)

dat <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric(),
                  time.collected = numeric())

for (i in 1:reps) {
  subset <- sim %>% filter(reps == i)
  for (j in 2:length(tms)) {
    time.before <- tms[j-1]
    time.collected <- tms[j]
    if (sum(subset$time <= time.collected  & subset$time > time.before) == 0){
      # dat[nrow(dat)+1, ] <- c(dat[nrow(dat), 1:4], time.collected)
      dat <- dat
    } else{
      ind <- max(which(subset$time <= time.collected  & subset$time > time.before))
      dat <- rbind(dat, cbind(subset[ind,], time.collected))
    }
    
  }
}


# sim %>% 
#   ggplot() +
#   geom_point(aes(x = time.steps, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
#   geom_point(aes(x = time.steps, y = sc.steps, col = "stem cell"), alpha = 0.5)+
#   labs(x = "time (in hours)", y = "cell counts", color = "",
#        title = "Simulated cell counts")+
#   scale_color_manual(values = c("differentiated cell" = "skyblue3",
#                                 "stem cell" = "indianred3")) +
#   theme_classic() + 
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14), 
#         axis.text.x = element_text(size = 14),
#         axis.title.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         legend.text = element_text(size = 13))

exp.mean.sc <- exp.mean.ec <-var.sc <- var.ec <- cov.sc.ec <- c()

for (i in 2:length(tms)) {
  exp.mean.sc <- c(exp.mean.sc,sc_func(tms[i]))
  exp.mean.ec <- c(exp.mean.ec,ec_func(tms[i]))
  var.sc <- c(var.sc, var_sc_func(tms[i]))
  var.ec <- c(var.ec, var_ec_func(tms[i]))
  cov.sc.ec <- c(cov.sc.ec, cov_func(tms[i]))
}

theo.dat <- data.frame(time.steps = tms[-1],
                       sc.steps = exp.mean.sc, 
                       ec.steps = exp.mean.ec,
                       sc.var = var.sc,
                       ec.var = var.ec, 
                       cov = cov.sc.ec)
ggplot() +
  geom_point(data = dat, aes(x = time.steps, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
  geom_point(data = dat, aes(x = time.steps, y = sc.steps, col = "stem cell"), alpha = 0.5)+
  geom_point(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_point(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_ribbon(data = theo.dat, aes(x = time.steps, 
                                   ymin = sc.steps-1.96*sqrt(sc.var), 
                                   ymax = sc.steps+1.96*sqrt(sc.var)), color = "gray", alpha = 0.25)+
  geom_ribbon(data = theo.dat, aes(x = time.steps, ymin = ec.steps-1.96*sqrt(ec.var), 
                                   ymax = ec.steps+1.96*sqrt(ec.var)), color = "gray", alpha = 0.25)+
  labs(x = "time (in hours)", y = "cell counts", color = "",
       title = "Simulated cell counts",
       subtitle = "Compound nonhomogeneous Poisson process")+
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



ggplot() +
  geom_point(data = dat, aes(x = time.collected, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
  geom_point(data = dat, aes(x = time.collected, y = sc.steps, col = "stem cell"), alpha = 0.5)+
  geom_point(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_point(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_ribbon(data = theo.dat, aes(x = time.steps, 
                                   ymin = sc.steps-1.96*sqrt(sc.var), 
                                   ymax = sc.steps+1.96*sqrt(sc.var)), color = "gray", alpha = 0.25)+
  geom_ribbon(data = theo.dat, aes(x = time.steps, ymin = ec.steps-1.96*sqrt(ec.var), 
                                   ymax = ec.steps+1.96*sqrt(ec.var)), color = "gray", alpha = 0.25)+
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

knitr::kable(dat %>% group_by(time.collected) %>%
               summarise(n = n(), 
                         samp.mean.sc = mean(sc.steps), 
                         samp.mean.ec = mean(ec.steps),
                         samp.var.sc = var(sc.steps),
                         samp.var.ec = var(ec.steps),
                         samp.cov = cov(sc.steps, ec.steps)) %>%
               mutate(theo.mean.sc = theo.dat$sc.steps[theo.dat$time.steps <= max(dat$time.collected)],
                      theo.mean.ec = theo.dat$ec.steps[theo.dat$time.steps <= max(dat$time.collected)],
                      theo.var.sc = theo.dat$sc.var[theo.dat$time.steps <= max(dat$time.collected)],
                      theo.var.ec = theo.dat$ec.var[theo.dat$time.steps <= max(dat$time.collected)],
                      theo.cov = theo.dat$cov[theo.dat$time.steps <= max(dat$time.collected)]))
# sample.vs.exp <- dat %>% group_by(time.collected) %>%
#   summarise(n = n(),
#             samp.mean.sc = mean(sc.steps),
#             samp.mean.ec = mean(ec.steps),
#             samp.var.sc = var(sc.steps),
#             samp.var.ec = var(ec.steps),
#             samp.cov = cov(sc.steps, ec.steps)) %>%
#   mutate(theo.mean.sc = theo.dat$sc.steps,
#          theo.mean.ec = theo.dat$ec.steps,
#          theo.var.sc = theo.dat$sc.var,
#          theo.var.ec = theo.dat$ec.var,
#          theo.cov = theo.dat$cov)
# 
# write.csv(sample.vs.exp, file = "sim_example.csv")



confidence_ellipse <- function(mu, Sigma, level = 0.95, npoints = 100) {
  # mu: mean vector (length 2)
  # Sigma: 2x2 covariance matrix
  # level: confidence level (default 95%)
  # npoints: number of points to draw the ellipse
  
  if (length(mu) != 2 || !all(dim(Sigma) == c(2, 2))) {
    stop("mu must be a length-2 vector and Sigma must be a 2x2 matrix.")
  }
  
  # Chi-squared quantile for the confidence level
  radius <- sqrt(qchisq(level, df = 2))
  
  # Eigen decomposition
  eig <- eigen(Sigma)
  angles <- seq(0, 2 * pi, length.out = npoints)
  
  # Unit circle scaled by eigenvalues and rotated by eigenvectors
  ellipse <- t(eig$vectors %*% diag(sqrt(eig$values)) %*% 
                 rbind(cos(angles), sin(angles))) * radius
  
  # Shift by mean
  ellipse <- sweep(ellipse, 2, mu, "+")
  
  colnames(ellipse) <- c("x", "y")
  return(as.data.frame(ellipse))
}

ellipse_dat <- data.frame(time.collected = numeric(),
                          x = numeric(), y = numeric())


for (j in 1:length(tms)) {
  xc <- theo.dat$sc.steps[j]
  yc <- theo.dat$ec.steps[j]
  cov.mat <- matrix(c(theo.dat$sc.var[j], theo.dat$cov[j],
                      theo.dat$cov[j], theo.dat$ec.var[j]),
                    nrow=2, ncol = 2)
  test.ellipse <- confidence_ellipse(mu = c(xc, yc), Sigma = cov.mat)
  ellipse_dat <- rbind(ellipse_dat, 
                       cbind(rep(tms[j+1], nrow = test.ellipse), test.ellipse))
}

colnames(ellipse_dat) <- c("time.collected", "sc.dim", "ec.dim")

# subset_dat <- dat %>% filter(time.collected <= 40)
# subset_ellipse_dat <- ellipse_dat%>% filter(time.collected <= 40)


plot_ly(dat, x = ~time.collected, y = ~sc.steps, z = ~ec.steps,  
        type = "scatter3d", mode = "markers", marker = list(size = 2),
        color = I("navyblue")) %>%
  add_trace(x = ellipse_dat$time.collected, 
            y = ellipse_dat$sc.dim, z = ellipse_dat$ec.dim, 
            opacity = 0.5, color = I("skyblue3"),
            marker = list(size = 4)) %>%
  layout(scene = list(xaxis = list(title = "time"),
                      yaxis = list(title = "stem cell"),
                      zaxis = list(title = "differentiated cell")))

###

tms.plot <- seq(5, 40, by = 0.2)

exp.mean.sc <- exp.mean.ec <-var.sc <- var.ec <- cov.sc.ec <- c()

for (i in 1:length(tms.plot)) {
  exp.mean.sc <- c(exp.mean.sc,sc_func(tms.plot[i]))
  exp.mean.ec <- c(exp.mean.ec,ec_func(tms.plot[i]))
  var.sc <- c(var.sc, var_sc_func(tms.plot[i]))
  var.ec <- c(var.ec, var_ec_func(tms.plot[i]))
  cov.sc.ec <- c(cov.sc.ec, cov_func(tms.plot[i]))
}

theo.dat.plot <- data.frame(time.steps = tms.plot,
                            sc.steps = exp.mean.sc,
                            ec.steps = exp.mean.ec,
                            sc.var = var.sc,
                            ec.var = var.ec,
                            cov = cov.sc.ec)

ellipse.dat.plot <- data.frame(time.collected = numeric(),
                               x = numeric(), y = numeric())


for (j in 1:length(tms.plot)) {
  xc <- theo.dat.plot$sc.steps[j]
  yc <- theo.dat.plot$ec.steps[j]
  cov.mat <- matrix(c(theo.dat.plot$sc.var[j], theo.dat.plot$cov[j],
                      theo.dat.plot$cov[j], theo.dat.plot$ec.var[j]),
                    nrow=2, ncol = 2)
  test.ellipse <- confidence_ellipse(mu = c(xc, yc), Sigma = cov.mat)
  ellipse.dat.plot <- rbind(ellipse.dat.plot, 
                            cbind(rep(tms.plot[j], nrow = test.ellipse), test.ellipse))
}

colnames(ellipse.dat.plot) <- c("time.collected", "sc.dim", "ec.dim")



plot_ly(ellipse.dat.plot, x = ~time.collected, y = ~sc.dim, z = ~ec.dim,  
        type = "scatter3d", mode = "markers",color = I("skyblue3"))%>%
  layout(scene = list(xaxis = list(title = "time"),
                      yaxis = list(title = "stem cell"),
                      zaxis = list(title = "differentiated cell")))



### simulate data from the multi-type branching process
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
  geom_point(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_point(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_ribbon(data = theo.dat, aes(x = time.steps, 
                                   ymin = sc.steps-1.96*sqrt(sc.var), 
                                   ymax = sc.steps+1.96*sqrt(sc.var)), color = "gray", alpha = 0.25)+
  geom_ribbon(data = theo.dat, aes(x = time.steps, ymin = ec.steps-1.96*sqrt(ec.var), 
                                   ymax = ec.steps+1.96*sqrt(ec.var)), color = "gray", alpha = 0.25)+
  labs(x = "time (in hours)", y = "cell counts", color = "",
       title = "Simulated cell counts",
       subtitle = "Branching process")+
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



ggplot() +
  geom_point(data = dat.markov, aes(x = time.collected, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
  geom_point(data = dat.markov, aes(x = time.collected, y = sc.steps, col = "stem cell"), alpha = 0.5)+
  geom_point(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_point(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
  geom_line(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
  geom_ribbon(data = theo.dat, aes(x = time.steps, 
                                   ymin = sc.steps-1.96*sqrt(sc.var), 
                                   ymax = sc.steps+1.96*sqrt(sc.var)), color = "gray", alpha = 0.25)+
  geom_ribbon(data = theo.dat, aes(x = time.steps, ymin = ec.steps-1.96*sqrt(ec.var), 
                                   ymax = ec.steps+1.96*sqrt(ec.var)), color = "gray", alpha = 0.25)+
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

knitr::kable(dat.markov %>% group_by(time.collected) %>%
               summarise(n = n(), 
                         samp.mean.sc = mean(sc.steps), 
                         samp.mean.ec = mean(ec.steps),
                         samp.var.sc = var(sc.steps),
                         samp.var.ec = var(ec.steps),
                         samp.cov = cov(sc.steps, ec.steps)) %>%
               mutate(theo.mean.sc = theo.dat$sc.steps,
                      theo.mean.ec = theo.dat$ec.steps,
                      theo.var.sc = theo.dat$sc.var,
                      theo.var.ec = theo.dat$ec.var,
                      theo.cov = theo.dat$cov))

# sample.vs.exp.markov <- dat.markov %>% group_by(time.collected) %>%
# summarise(n = n(), 
#           samp.mean.sc = mean(sc.steps), 
#           samp.mean.ec = mean(ec.steps),
#           samp.var.sc = var(sc.steps),
#           samp.var.ec = var(ec.steps),
#           samp.cov = cov(sc.steps, ec.steps)) %>%
#   mutate(theo.mean.sc = theo.dat$sc.steps,
#          theo.mean.ec = theo.dat$ec.steps,
#          theo.var.sc = theo.dat$sc.var,
#          theo.var.ec = theo.dat$ec.var,
#          theo.cov = theo.dat$cov)
# 
# write.csv(sample.vs.exp.markov, file = "sim_markov_example.csv")

plot_ly(dat.markov, x = ~time.collected, y = ~sc.steps, z = ~ec.steps,  
        type = "scatter3d", mode = "markers", marker = list(size = 2),
        color = I("navyblue")) %>%
  add_trace(x = ellipse_dat$time.collected, 
            y = ellipse_dat$sc.dim, z = ellipse_dat$ec.dim, 
            opacity = 0.5, color = I("skyblue3"),
            marker = list(size = 4)) %>%
  layout(scene = list(xaxis = list(title = "time"),
                      yaxis = list(title = "stem cell"),
                      zaxis = list(title = "differentiated cell")))
upper <- 80
lower <- 0

df.time <- data.frame(time = c(sim$time.steps[sim$time.steps >lower & sim$time.steps < upper],
                               sim.markov$time.steps[sim.markov$time.steps >lower & sim.markov$time.steps < upper]),
                      group = c(rep("cpnn", length(sim$time.steps[sim$time.steps >lower & sim$time.steps < upper])),
                                rep("branching",length(sim.markov$time.steps[sim.markov$time.steps >lower & sim.markov$time.steps < upper]))))

ggplot(df.time) +
  geom_histogram(aes(x = time, fill = group),binwidth = 1, color = "black") +
  facet_wrap(~group, ncol = 1) +  
  labs(title = "Simulated Event Times", x = "Simulated event time", y = "Count") +
  theme_minimal()

ggplot(df.time) +
  geom_histogram(aes(x = time, fill = group),binwidth = 0.2,
                 position = "identity", alpha = 0.5, color = "black") +
  labs(title = "Simulated Event Times", x = "Simulated event time", y = "Count") +
  theme_minimal()

index <- 15
df.cell <- data.frame(stem = c(dat$sc.steps[dat$time.collected == tms[index]],
                               dat.markov$sc.steps[dat.markov$time.collected == tms[index]]),
                      diff = c(dat$ec.steps[dat$time.collected == tms[index]],
                               dat.markov$ec.steps[dat.markov$time.collected == tms[index]]),
                      group = c(rep("cnpp", length(dat$sc.steps[dat$time.collected == tms[index]])),
                                rep("branching", length(dat.markov$sc.steps[dat.markov$time.collected == tms[index]]))))
ggplot(df.cell) +
  geom_point(aes(x = stem, y = diff, color = group))+
  labs(x = "stem cell", y = "differentiated cell",
       title = paste0("timepoint ",tms[index]))

dim(sim)
dim(sim.markov)

### stop at fixed time
# tau <- 50
# set.seed(seed0)
# sim.tau <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric())
# for (i in 1:reps) {
#   # draw event time from the nhppp package using thinning algorithm
#   time.draw <- draw_intensity(lambda = intensity_func,
#                               line_majorizer_intercept = intercept,
#                               line_majorizer_slope = slope,
#                               line_majorizer_is_loglinear = TRUE,
#                               t_min = 0, t_max = t_max)
#   sc.current <- s0.test
#   ec.current <- dc.current <- time.current <- 0
#   sc.steps <- c(s0.test)
#   ec.steps <- c(ec.current)
#   # dc.steps <- c(dc.current)
#   time.steps <- c(time.current)
#   itnum <- 0
#   index <- 0
#   while (index < length(time.draw)) {
#     itnum <- itnum + 1
#     index <- index + 1
#     time.current <- time.draw[index]
#     # calculate proliferation function at event time
#     q_t <- q_t1(time.current, prm.test)
#     if ((q_t < min_q_t| (itnum > max_it))) {
#       cat(q_t, "\n")
#       break
#     }
#     # determine the probability of each type of division
#     probs <- c((2-q_t)/k, max(0, (k-2)*(2-q_t)/k), max(0,(2-k + (k-1)*q_t)/k))
#     if (sum(probs) > 1) {
#       probs <- probs/sum(probs)
#     }
#     # sample the type of division
#     event_type <- sample(1:3, 1, prob=probs)
#     switch(event_type,
#            "1" = {
#              sc.current <- sc.current + 1
#              ec.current <- ec.current
#            },
#            "2" = {
#              sc.current <- sc.current
#              ec.current <- ec.current + 1
#            },
#            "3" = {
#              sc.current <- sc.current - 1
#              ec.current <- ec.current + 2
#            })
#     sc.steps <- c(sc.steps, sc.current)
#     ec.steps <- c(ec.steps, ec.current)
#     time.steps <- c(time.steps, time.current)
#   }
#   sim.tau <- rbind(sim.tau,cbind(reps = rep(i, length(sc.steps)), time.steps, sc.steps, ec.steps))
# }
# 
# dat.tau <- data.frame(reps = numeric(), time = numeric(), sc = numeric(), ec = numeric(),
#                   time.collected = numeric())
# 
# for (i in 1:reps) {
#   subset <- sim.tau %>% filter(reps == i)
#   for (j in 2:length(tms)) {
#     time.before <- tms[j-1]
#     time.collected <- tms[j]
#     if (sum(subset$time <= time.collected  & subset$time > time.before) == 0){
#       # dat[nrow(dat)+1, ] <- c(dat[nrow(dat), 1:4], time.collected)
#       dat.tau <- dat.tau
#     } else{
#       ind <- max(which(subset$time <= time.collected  & subset$time > time.before))
#       dat.tau <- rbind(dat.tau, cbind(subset[ind,], time.collected))
#     }
#     
#   }
# }
# 
# 
# ggplot() +
#   geom_point(data = dat.tau, aes(x = time.steps, y = ec.steps, col = "differentiated cell"), alpha = 0.5)+
#   geom_point(data = dat.tau, aes(x = time.steps, y = sc.steps, col = "stem cell"), alpha = 0.5)+
#   geom_point(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
#   geom_point(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
#   geom_line(data = theo.dat, aes(x = time.steps, y = ec.steps, col = "expected differentiated cell"))+
#   geom_line(data = theo.dat, aes(x = time.steps, y = sc.steps, col = "expected stem cell"))+
#   geom_ribbon(data = theo.dat, aes(x = time.steps, 
#                                    ymin = sc.steps-1.96*sqrt(sc.var), 
#                                    ymax = sc.steps+1.96*sqrt(sc.var)), color = "gray", alpha = 0.25)+
#   geom_ribbon(data = theo.dat, aes(x = time.steps, ymin = ec.steps-1.96*sqrt(ec.var), 
#                                    ymax = ec.steps+1.96*sqrt(ec.var)), color = "gray", alpha = 0.25)+
#   labs(x = "time (in hours)", y = "cell counts", color = "",
#        title = "Simulated cell counts")+
#   scale_color_manual(values = c("differentiated cell" = "skyblue3",
#                                 "stem cell" = "indianred1",
#                                 "expected differentiated cell" = "navy",
#                                 "expected stem cell" = "brown4")) +
#   theme_classic() +
#   theme(legend.position = "bottom",
#         plot.title = element_text(size = 14),
#         axis.text.x = element_text(size = 14),
#         axis.title.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title.y = element_text(size = 14),
#         legend.text = element_text(size = 13))
# 
# plot_ly(dat.tau, x = ~time.collected, y = ~sc.steps, z = ~ec.steps,  
#         type = "scatter3d", mode = "markers", marker = list(size = 2),
#         color = I("navyblue")) %>%
#   add_trace(x = ellipse_dat$time.collected, 
#             y = ellipse_dat$sc.dim, z = ellipse_dat$ec.dim, 
#             opacity = 0.5, color = I("skyblue3"),
#             marker = list(size = 4)) %>%
#   layout(scene = list(xaxis = list(title = "time"),
#                       yaxis = list(title = "stem cell"),
#                       zaxis = list(title = "differentiated cell")))
# 
# sim.tau %>% group_by(reps) %>% summarize(n = n(), max.time = max(time.steps), min.sc = min(sc.steps))
# 
# prob_func <- function(time.jump, jump.size, prm){
#   q_t <- q_t1(tms = time.jump, prm = prm)
#   probs <- c((2-q_t)/k, max(0, (k-2)*(2-q_t)/k), max(0,(2-k + (k-1)*q_t)/k))
#   if (sum(probs) > 1) {
#     probs <- probs/sum(probs)
#   }
#   if(sum(jump.size == c(1,0)) > 1){
#     out <- probs[1]
#   } else if(sum(jump.size == c(0,1))>1) {
#     out <- probs[2]
#   } else if(sum(jump.size == c(-1,2))>1){
#     out <- probs[3]
#   } else {
#     out <- 0
#   }
# }
# 
# 
# sim.tau.time.jump <- sim.tau %>% filter(reps == 1) %>%
#   pull(time.steps)
# sim.tau.sc.steps <- sim.tau %>% filter(reps == 1) %>%
#   pull(sc.steps)
# sim.tau.ec.steps <- sim.tau %>% filter(reps == 1) %>%
#   pull(ec.steps)
# diff.sc <- diff(sim.tau.sc.steps)
# diff.ec <- diff(sim.tau.ec.steps)
# 
# sim.tau.jump.size<- lapply(seq_along(diff.sc), function(i) c(diff.sc[i], diff.ec[i]))
# print(sim.tau.jump.size)
# length(sim.tau.jump.size)
# 
# prob.vec <- c()
# for (i in 1:length(sim.tau.jump.size)) {
#   probs <- prob_func(time.jump = sim.tau.time.jump[i+1], jump.size = sim.tau.jump.size[[i]], prm = prm.test)
#   prob.vec <- c(prob.vec, probs)
# }
# 
# log_likelihood <- function(stop.time = 10, prm, rate, time.jump, jump, k =3){
#   new_intensity_func <- function(t){
#     out <- rate*s0.test*exp(rate*(t-int_q(t,prm = prm)+int_q(0,prm = prm)))
#     return(out)
#   }
#   cummulative.intensity <- integrate(new_intensity_func, 0, stop.time)$value
#   log.likelihood <- (-1)*cummulative.intensity
#   for (i in length(time.jump)) {
#     intesity.time <- new_intensity_func(time.jump[i])
#     q_t <- q_t1(tms = time.jump[i], prm = prm)
#     probs <- c((2-q_t)/k, max(0, (k-2)*(2-q_t)/k), max(0,(2-k + (k-1)*q_t)/k))
#     if (sum(probs) > 1) {
#       probs <- probs/sum(probs)
#     }
#     jump.size <- jump[[i]]
#     if(sum(jump.size == c(1,0)) > 1){
#       out.prob <- probs[1]
#     } else if(sum(jump.size == c(0,1))>1) {
#       out.prob <- probs[2]
#     } else if(sum(jump.size == c(-1,2))>1){
#       out.prob <- probs[3]
#     } else {
#       out.prob <- 0
#     }
#     log.likelihood <- log.likelihood + log(intesity.time) + log(out.prob)
#   }
#   return(log.likelihood)
# }
# 
# 
# log_likelihood(stop.time = tau, prm = prm.test, rate = r.test, 
#                time.jump = sim.tau.time.jump,
#                jump = sim.tau.jump.size)
# 
# sim.tau.time.jump <- diff.sc <- diff.ec <- c()
# 
# for (i in 1:reps) {
#   sim.subset <- sim.tau %>% filter(reps == i) %>% filter(time.steps <= tau)
#   sim.subset.time.jump <- sim.subset %>% pull(time.steps)
#   sim.subset.time.jump <- sim.subset.time.jump[-1]
#   sim.subset.sc.steps <- sim.subset %>% pull(sc.steps)
#   sim.subset.ec.steps <- sim.subset %>% pull(ec.steps)
#   diff.subset.sc <- diff(sim.subset.sc.steps)
#   diff.subset.ec <- diff(sim.subset.ec.steps)
#   
#   sim.tau.time.jump <- c(sim.tau.time.jump, sim.subset.time.jump)
#   diff.sc <- c(diff.sc, diff.subset.sc)
#   diff.ec <- c(diff.ec, diff.subset.ec)
# }
# sim.tau.jump.size<- lapply(seq_along(diff.sc), function(i) c(diff.sc[i], diff.ec[i]))
# length(sim.tau.time.jump)
# length(sim.tau.jump.size)
# length(diff.sc)
# 
# r.vals <- seq(0.1, 1, length.out = 10)
# a1.vals <- seq(1.1,2, length.out = 10)
# a2.vals <- seq(-0.120, 0.006, length.out = 10)
# a3.vals <- seq(0.004, 0.013, length.out = 10)
# k.vals <- seq(0.5, 5, length.out = 10)
# 
# grid3D <- expand.grid(y = a1.vals, z = a2.vals, w = a3.vals)
# grid4D <- expand.grid(x = r.vals, y = a1.vals, z = a2.vals, w = a3.vals)
# grid5D <- expand.grid(x = r.vals, y = a1.vals, z = a2.vals, w = a3.vals, v = k.vals)
# 
# log.likelihood <- c()
# for (i in 1:nrow(grid5D)) {
#   r.input <- grid5D$x[i]
#   prm.input <- c(grid5D$y[i], grid5D$z[i], grid5D$w[i])
#   k.input <- grid5D$v[i]
#   output.lik <- log_likelihood(stop.time = 50, prm = prm.input, rate = r.input, 
#                                time.jump = sim.tau.time.jump,
#                                jump = sim.tau.jump.size,
#                                k = k.input)
#   log.likelihood <- c(log.likelihood, output.lik)
# }
# which.max(log.likelihood)
# grid5D[which.max(log.likelihood),]
# 
# 
# log.likelihood <- c()
# for (i in 1:nrow(grid4D)) {
#   r.input <- grid4D$x[i]
#   prm.input <- c(grid4D$y[i], grid4D$z[i], grid4D$w[i])
#   output.lik <- log_likelihood(stop.time = 50, prm = prm.input, rate = r.input, 
#                                time.jump = sim.tau.time.jump,
#                                jump = sim.tau.jump.size,
#                                k = 3)
#   log.likelihood <- c(log.likelihood, output.lik)
# }
# which.max(log.likelihood)
# grid4D[which.max(log.likelihood),]
# max(log.likelihood)
# 
# 
# log.likelihood <- c()
# for (i in 1:nrow(grid3D)) {
#   prm.input <- c(grid4D$y[i], grid4D$z[i], grid4D$w[i])
#   output.lik <- log_likelihood(stop.time = 50, prm = prm.input, rate = r.test, 
#                                time.jump = sim.tau.time.jump,
#                                jump = sim.tau.jump.size,
#                                k = 3)
#   log.likelihood <- c(log.likelihood, output.lik)
# }
# which.max(log.likelihood)
# grid3D[which.max(log.likelihood),]
# max(log.likelihood)
