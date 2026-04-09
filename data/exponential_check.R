library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(goftest)
library(fitdistrplus)

## check if the interarrival time is exponential with constant rate
mep_tpo_count <- read.csv("mep_tpo_count.csv")
mep_tpo_count <- mep_tpo_count[-c(1:4),]
head(mep_tpo_count)
tpo_interarrival <- diff(mep_tpo_count$time)
tpo_scale_interarrival <- tpo_interarrival*mep_tpo_count$mep[-nrow(mep_tpo_count)]
tpo_scale_interarrival <- tpo_scale_interarrival[tpo_scale_interarrival != 0]

lambda_hat <- 1 / mean(tpo_scale_interarrival)


# convert to data frame
df.arrival <- data.frame(arr.time = tpo_scale_interarrival,
                         scale.arr.time = lambda_hat*tpo_scale_interarrival)

# plot
hist.arr.plot <- ggplot(df.arrival, aes(x = scale.arr.time)) +
  geom_histogram(breaks = c(h$breaks,6.5),
                 aes(y = after_stat(density)),
                 # bins = 11,
                 fill = "gray",
                 color = "black") +
  stat_function(fun = dexp,
                args = list(rate = 1),
                color = "indianred3",
                linewidth = 1) +
  scale_x_continuous(limits = c(0,6.5)) +
  labs(title = "Lacking TPO condition",
       x = "Scaled interarrival times",
       y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text.align = 0,
    legend.position = "bottom"
  )
hist.arr.plot
##
## formal test
ad.test(tpo_scale_interarrival, null = "pexp", rate = lambda_hat)
gofstat(fitdist(tpo_scale_interarrival, "exp"))

##
##


###
###
## control dataset
mep_ctrl_count <- read.csv("mep_ctrl_count.csv")
mep_ctrl_count <- mep_ctrl_count[-c(1:14),]

ctrl_interarrival <- diff(mep_ctrl_count$time)
ctrl_scale_interarrival <- ctrl_interarrival*mep_ctrl_count$mep[-nrow(mep_ctrl_count)]
ctrl_scale_interarrival <- ctrl_scale_interarrival[ctrl_scale_interarrival != 0]

lambda_hat_ctrl <- 1 / mean(ctrl_scale_interarrival)


# convert to data frame
df.arrival.ctrl <- data.frame(arr.time = ctrl_scale_interarrival,
                              scale.arr.time = lambda_hat_ctrl*ctrl_scale_interarrival)

# plot
h_ctrl <- hist(df.arrival.ctrl$scale.arr.time, plot = FALSE, 
               breaks = c(h$breaks,6.5))
hist(df.arrival.ctrl$scale.arr.time, breaks = c(h$breaks,6.5))
hist.arr.ctrl.plot <- ggplot(df.arrival.ctrl, aes(x = scale.arr.time)) +
  geom_histogram(aes(y = after_stat(density)),
                 breaks = h_ctrl$breaks,
                 fill = "gray",
                 color = "black") +
  stat_function(fun = dexp,
                args = list(rate = 1),
                color = "indianred3",
                linewidth = 1) +
  labs(title = "Control condition",
       x = "Scaled interarrival times",
       y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text.align = 0,
    legend.position = "bottom"
  )
hist.arr.ctrl.plot

## formal test
ad.test(ctrl_scale_interarrival, null = "pexp", rate = lambda_hat_ctrl)
gofstat(fitdist(ctrl_scale_interarrival, "exp"))


(hist.arr.ctrl.plot | hist.arr.plot) +  plot_annotation(
  title = "Scaled Arrival Times of MEP Divisions")+ 
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )