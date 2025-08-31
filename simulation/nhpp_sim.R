library(nhppp)
library(ggplot2)
### Method 1

q_t1 <- function(tms, prm) {
  return(2 - 1/(prm[1] + prm[2]*tms + prm[3]*tms^2))
}


int_q <- function(t, prm) {
  2*t - (2*atan((prm[2] + 2*prm[3]*t)/sqrt(4*prm[1]*prm[3] - prm[2]^2))/
    sqrt(4*prm[1]*prm[3] - prm[2]^2))
}

s0 <- 200
sc.current <- s0
ec.current <- dc.current <- time.current <- 0
r.test <- 0.2
sc.steps <- c(s0)
ec.steps <- c(ec.current)
dc.steps <- c(dc.current)
time.steps <- c(time.current)
prm <-c(1.2, -0.106, 0.005)
min_q_t <- 0.0001
max_it = 50000
k <- 3
itnum <- 0
seed0 <- 712982

set.seed(seed0)
while (sc.current > 0) {
  itnum <- itnum + 1
  time.current <- ppp_next_n(n=1, rate= r.test*sc.current,
                             t_min = time.current)
  q_t <- q_t1(time.current, prm)
  if ((q_t < min_q_t| (itnum > max_it))) {
    cat(q_t, "\n")
    break
  }
  probs <- c((2-q_t)/k, max(0, (k-2)*(2-q_t)/k), max(0,(2-k + (k-1)*q_t)/k))
  if (sum(probs) > 1) {
    probs <- probs/sum(probs)
  }
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
}

sc.m1 <- sc.steps
ec.m1 <- ec.steps
time.m1 <- time.steps
plot(time.m1, sc.m1)
plot(time.m1, ec.m1)


### Method 2

s0 <- 200
sc.current <- s0
ec.current <- dc.current <- time.current <- 0
r.test <- 0.2
sc.steps <- c(s0)
ec.steps <- c(ec.current)
dc.steps <- c(dc.current)
time.steps <- c(time.current)
prm <-c(1.2, -0.106, 0.005)
min_q_t <- 0.0001
max_it = 50000
k <- 3
itnum <- 0
seed0 <- 712982

slope <- -0.01
intercept <- 6

intensity_func <- function(t){
  out <- r.test*s0*exp(r.test*(t-int_q(t,prm)+int_q(0,prm)))
}

ec_func <- function(t){
  integrand <- function(t){
    out <- q_t1(tms = t, prm = prm)*intensity_func(t)
  }
  out.ec <- integrate(integrand, lower = 0, upper = t)$value
  return(out.ec)
}
sc_func <- function(t){
  out <- s0*exp(r.test*(t-int_q(t,prm)+int_q(0,prm)))
  return(out)
}

set.seed(seed0)
time.draw <- draw_intensity(lambda = intensity_func,
                            line_majorizer_intercept = intercept,
                            line_majorizer_slope = slope,
                            line_majorizer_is_loglinear = TRUE,
                            t_min = 0, t_max = 70)
while (sc.current > 0) {
  itnum <- itnum + 1
  time.current <- time.draw[itnum]
  q_t <- q_t1(time.current, prm)
  if ((q_t < min_q_t| (itnum > max_it))) {
    cat(q_t, "\n")
    break
  }
  probs <- c((2-q_t)/k, max(0, (k-2)*(2-q_t)/k), max(0,(2-k + (k-1)*q_t)/k))
  if (sum(probs) > 1) {
    probs <- probs/sum(probs)
  }
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
}

sc.m2 <- sc.steps
ec.m2 <- ec.steps
time.m2 <- time.steps
plot(time.m2, sc.m2)
plot(time.m2, ec.m2)


intensity.steps.m1 <- c()
intensity.steps.m2 <- c()
ec.mean <- c()
sc.mean <- c()
majorizer <- c()
for (i in 1:length(time.m1)) {
  intensity.steps.m1 <- c(intensity.steps.m1, intensity_func(time.m1[i]))
  majorizer <- c(majorizer, exp(intercept+slope*time.m1[i]))
}
for (i in 1:length(time.draw)) {
  intensity.steps.m2 <- c(intensity.steps.m2, intensity_func(time.draw[i]))
  ec.mean <- c(ec.mean, ec_func(time.draw[i]))
  sc.mean <- c(sc.mean, sc_func(time.draw[i]))
}

plot(time.m1, intensity.steps.m1, type = "l", 
     xlab = "time", ylab = "intensity")
plot(time.draw, intensity.steps.m2, type = "l", 
     xlab = "time", ylab = "intensity")

plot(time.m1, intensity_steps)
plot(time.m1, majorizer)
plot(time.m1, majorizer - intensity_steps)
min(majorizer - intensity_steps)



data.m1 <- data.frame(time = time.m1,
                      sc = sc.m1,
                      ec = ec.m1)
data.m2 <- data.frame(time = time.m2,
                      sc = sc.m2,
                      ec = ec.m2,
                      ec.mean = ec.mean[1:length(time.m2)],
                      sc.mean = sc.mean[1:length(time.m2)])
color_values <- c("simulated EC" = "skyblue3", 
                  "simulated SC" = "indianred2",
                  "expected counts" = "black")
ggplot() +
  # geom_point(data = data.m1, aes(x = time, y = ec, color = "EC by method 1"))+
  geom_point(data = data.m2, aes(x = time, y = ec, color = "simulated EC"))+
  # geom_point(data = data.m1, aes(x = time, y = sc, color = "SC by method 1"))+
  geom_point(data = data.m2, aes(x = time, y = sc, color = "simulated SC"))+
  geom_line(data = data.m2, aes(x = time, y = ec.mean, color = "expected counts"), lty = "dashed") +
  geom_line(data = data.m2, aes(x = time, y = sc.mean, color = "expected counts"), lty = "dashed") +
  scale_color_manual(values = color_values) +
  labs(x = "time (in hours)", y = "cell counts", color = "",
       title = "Simulated cell counts")+
  theme_classic() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14), 
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13))



hist(time.m1)
hist(time.m2)
hist(time.draw)
