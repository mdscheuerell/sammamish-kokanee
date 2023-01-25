#### setup ####

library(here)
library(readr)
library(dplyr)
library(MARSS)

## read data
dat_counts <- read_csv(here("data", "samm_kokanee_counts.csv"))

## log-transform counts
log_counts <- dat_counts %>%
  select(count) %>%
  log() %>%
  unlist()


#### model fitting ####

control_list <- list(maxit = 5000, allow.degen = TRUE)
  
## 1) unbiased random walk

mod_list <- list(
  B = matrix(1),
  U = matrix(0),
  Q = matrix("q"),
  Z = matrix(1),
  A = matrix(0),
  R = matrix("r")
)

mod_fit_unbiased <- MARSS(y = matrix(log_counts, nrow = 1),
                          model = mod_list, control = control_list)


## 2) biased random walk

mod_list <- list(
  B = matrix(1),
  U = matrix("u"),
  Q = matrix("q"),
  Z = matrix(1),
  A = matrix(0),
  R = matrix("r")
)

mod_fit_biased <- MARSS(y = matrix(log_counts, nrow = 1),
                        model = mod_list, control = control_list)

## boostrapped CI's
mod_fit_biased_CI <- MARSSparamCIs(mod_fit_biased, method = "hessian", alpha = 0.05, nboot = 1000)


## 3) stationary AR(1)

mod_list <- list(
  B = matrix("b"),
  U = matrix(0),
  Q = matrix("q"),
  Z = matrix(1),
  A = matrix(0),
  R = matrix("r")
)

mod_fit_ar1 <- MARSS(y = matrix(log_counts, nrow = 1),
                     model = mod_list, control = control_list)

